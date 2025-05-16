# --- Download Functions (download.R) ---

#' Download Helper using wget or axel
#'
#' @param url The URL to download.
#' @param dest_file The destination file path.
#' @param parallel_n Number of parallel connections for axel (0 means use wget).
#' @param software_paths List containing paths to 'wget' and 'axel'.
#' @keywords internal
execute_download <- function(url, dest_file, parallel_n = 0, software_paths) {
  fs::dir_create(fs::path_dir(dest_file)) # Ensure directory exists

  if (parallel_n > 0 && !is.null(software_paths$axel)) {
    # Use axel
    run_command(software_paths$axel,
      args = c("-n", as.character(parallel_n), "-o", dest_file, "-a", "-c", url),
      spinner = TRUE
    )
  } else if (!is.null(software_paths$wget)) {
    # Use wget
    # Check wget version for progress bar style - complex, maybe just use quiet?
    # Let's use quiet + show-progress if possible, otherwise default.
    # Detecting version accurately is tricky. Assume modern wget supports --show-progress.
    # Older versions might need --progress=bar
    # For simplicity, let's just use -c (continue) and -O (output file). Progress might show by default.
    run_command(software_paths$wget,
      args = c("-c", url, "-O", dest_file, "--quiet", "--show-progress"), # Best guess for modern wget
      spinner = FALSE
    ) # Wget has its own progress
  } else {
    cli::cli_abort("Neither {.val wget} nor {.val axel} path is configured. Cannot download.")
  }
}

#' Download Helper using Aspera
#'
#' @param remote_path The Aspera remote path (e.g., era-fasp@...).
#' @param db Source database ("ENA" or "GSA").
#' @param software_paths List containing path to 'ascp'.
#' @param aspera_key_path Path to the Aspera key file.
#' @param output_dir Local directory to download into.
#' @keywords internal
execute_aspera <- function(remote_path, db = "ENA", software_paths, aspera_key_path, output_dir = ".") {
  if (is.null(software_paths$ascp)) {
    cli::cli_abort("Aspera download requested, but {.val ascp} path is not configured.")
  }
  if (is.null(aspera_key_path) || !fs::file_exists(aspera_key_path)) {
    cli::cli_abort("Aspera key file not found or not specified: {.path {aspera_key_path}}")
  }

  # Standard Aspera options from the script
  ascp_args <- c(
    "-P", "33001", # Port
    "-i", aspera_key_path, # Key file
    "-QT", # Quality of Transport (encryption/resume)
    "-l", "1000m", # Max rate (adjust?)
    "-k1", # Resume level (check documentation)
    "-d", # Create directory (usually needed)
    remote_path, # Source path
    output_dir # Destination directory
  )

  run_command(software_paths$ascp, args = ascp_args, spinner = TRUE)
}

#' Get Aspera Key Path
#'
#' Finds the appropriate Aspera key file based on standard locations or options.
#'
#' @param db "ENA" or "GSA".
#' @param software_paths List containing path to 'ascp'.
#' @return Path to the key file or NULL.
#' @keywords internal
get_aspera_key_path <- function(db = "ENA", software_paths) {
  key_path <- NULL
  option_name <- ifelse(db == "ENA", "iseq.aspera_key_ena", "iseq.aspera_key_gsa")
  option_path <- getOption(option_name)

  if (!is.null(option_path) && fs::file_exists(option_path)) {
    cli::cli_inform("Using Aspera key from R option {.val {option_name}}: {.path {option_path}}")
    return(option_path)
  }

  if (db == "ENA") {
    if (!is.null(software_paths$ascp)) {
      ascp_dir <- fs::path_dir(software_paths$ascp)
      # Common locations relative to ascp executable
      potential_paths <- c(
        fs::path(ascp_dir, "..", "etc", "aspera", "aspera_bypass_rsa.pem"),
        fs::path(ascp_dir, "..", "etc", "aspera_tokenauth_id_rsa")
      )
      for (p in potential_paths) {
        if (fs::file_exists(p)) {
          key_path <- fs::path_real(p) # Get absolute path
          cli::cli_inform("Found standard ENA Aspera key: {.path {key_path}}")
          break
        }
      }
    }
    if (is.null(key_path)) {
      cli::cli_alert_warning(
        "Could not find standard ENA Aspera key. Set {.code options({option_name} = ...)} if needed."
      )
    }
  } else { # GSA
    gsa_key_file <- ".asperaGSA.openssh" # Relative to working dir or output dir? Assume output_dir
    # Should be in the output_dir where the script runs
    potential_path <- fs::path(getwd(), gsa_key_file) # Or use output_dir argument if passed down

    if (fs::file_exists(potential_path)) {
      key_path <- potential_path
      cli::cli_inform("Found GSA Aspera key: {.path {key_path}}")
    } else {
      # Try downloading it as per bash script
      cli::cli_alert_info("GSA Aspera key {.file {gsa_key_file}} not found, attempting download...")
      gsa_key_url <- "https://ngdc.cncb.ac.cn/gsa/file/downFile?fileName=download/aspera01.openssh"
      tryCatch(
        {
          resp <- httr::GET(gsa_key_url, httr::write_disk(potential_path, overwrite = TRUE), httr::timeout(60))
          # Check if download was successful (may not return 200 OK for file downloads)
          if (fs::file_exists(potential_path) && fs::file_size(potential_path) > 0) {
            cli::cli_alert_success("Downloaded GSA Aspera key to {.path {potential_path}}")
            key_path <- potential_path
          } else {
            stop("Download failed or resulted in empty file.")
          }
        },
        error = function(e) {
          cli::cli_alert_warning("Could not download GSA Aspera key: {e$message}")
          if (fs::file_exists(potential_path)) fs::file_delete(potential_path)
        }
      )
    }
    if (is.null(key_path)) {
      cli::cli_alert_warning(
        "Could not find or download GSA Aspera key. Set {.code options({option_name} = ...)} if needed."
      )
    }
  }
  return(key_path)
}


#' Download a Single ENA/SRA Run
#'
#' Handles logic for downloading FASTQ or SRA, using FTP/HTTPS/Aspera.
#'
#' @param run_id SRR/ERR/DRR identifier.
#' @param run_info A single row tibble/data.frame from the metadata file for this run.
#' @param download_fastq_gz Should FASTQ.gz be preferred?
#' @param use_aspera Should Aspera be used if available?
#' @param database_preference "ena" or "sra".
#' @param parallel_n Parallel connections for axel.
#' @param output_dir Destination directory.
#' @param software_paths Paths to tools.
#' @return A list containing paths to downloaded files (or the SRA file path).
#'   Keys: `files`, `type` ("fastq_gz", "sra").
#' @keywords internal
download_ena_sra_run <- function(run_id, run_info, download_fastq_gz = TRUE, use_aspera = FALSE,
                                 database_preference = "ena", parallel_n = 0, output_dir = ".", software_paths) {
  downloaded_files <- list(files = character(), type = NULL)
  aspera_key_path <- NULL
  if (use_aspera) {
    aspera_key_path <- get_aspera_key_path("ENA", software_paths)
    if (is.null(aspera_key_path)) {
      cli::cli_alert_warning("Aspera requested but key not found/configured for ENA, falling back to FTP/HTTPS.")
      use_aspera <- FALSE
    }
  }

  # --- Try FASTQ.gz first if requested and available from ENA ---
  # ENA metadata TSV has columns like 'fastq_ftp', 'fastq_bytes', 'fastq_md5', 'fastq_aspera'
  # These columns can contain multiple semicolon-separated values for paired-end data.
  fastq_ftp_links <- NULL
  fastq_aspera_links <- NULL
  if (download_fastq_gz && database_preference == "ena" && "fastq_ftp" %in% names(run_info)) {
    fastq_ftp_links <- str_split(run_info$fastq_ftp, ";")[[1]] %>%
      discard(is.na) %>%
      discard(. == "")
    if (use_aspera && "fastq_aspera" %in% names(run_info)) {
      fastq_aspera_links <- str_split(run_info$fastq_aspera, ";")[[1]] %>%
        discard(is.na) %>%
        discard(. == "")
      # Ensure ftp and aspera links correspond if possible (e.g., same number)
      if (length(fastq_aspera_links) != length(fastq_ftp_links)) {
        cli::cli_alert_warning(
          "Mismatch between number of FTP ({length(fastq_ftp_links)}) and Aspera ({length(fastq_aspera_links)}) FASTQ links for {.val {run_id}}. Preferring FTP.")
        fastq_aspera_links <- NULL # Fallback to FTP
      }
    } else {
      fastq_aspera_links <- NULL # No aspera links or not requested
    }
  }

  if (length(fastq_ftp_links) > 0) {
    cli::cli_alert_info("Attempting direct FASTQ.gz download for {.val {run_id}} from ENA...")
    download_mode <- ifelse(use_aspera && !is.null(fastq_aspera_links), "Aspera", "HTTPS")
    cli::cli_inform("Download mode: {.val {download_mode}}")

    success <- TRUE
    temp_files <- character()

    links_to_use <- if (download_mode == "Aspera") fastq_aspera_links else fastq_ftp_links

    for (i in seq_along(links_to_use)) {
      link <- links_to_use[i]
      # Construct destination path from link basename
      dest_file <- fs::path(output_dir, fs::path_file(link)) # Assumes link ends with filename
      # Use FTP link for basename even if downloading via Aspera
      if (download_mode == "Aspera") dest_file <- fs::path(output_dir, fs::path_file(fastq_ftp_links[i]))

      temp_files <- c(temp_files, dest_file)
      cli::cli_inform("Downloading: {.url {link}} to {.path {dest_file}}")

      tryCatch(
        {
          if (download_mode == "Aspera") {
            # Aspera link format is usually like: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR...
            execute_aspera(
              remote_path = glue::glue("era-fasp@{link}"), db = "ENA", software_paths = software_paths,
              aspera_key_path = aspera_key_path, output_dir = output_dir
            )
          } else {
            # FTP link format: ftp.sra.ebi.ac.uk/vol1/fastq/SRR...
            execute_download(
              url = paste0("https://", link), dest_file = dest_file,
              parallel_n = parallel_n, software_paths = software_paths
            )
          }
          # Basic check after download
          if (!fs::file_exists(dest_file) || fs::file_size(dest_file) == 0) {
            stop("Download resulted in missing or empty file.")
          }
        },
        error = function(e) {
          cli::cli_alert_danger("Failed to download {.url {link}}: {e$message}")
          success <<- FALSE # Set success flag to FALSE
          if (fs::file_exists(dest_file)) fs::file_delete(dest_file) # Clean up failed part
        }
      )
      if (!success) break # Stop if one file failed
    }

    if (success) {
      downloaded_files$files <- temp_files
      downloaded_files$type <- "fastq_gz"
      cli::cli_alert_success("Successfully downloaded FASTQ.gz file(s) for {.val {run_id}}.")
      return(downloaded_files) # Success, return downloaded FASTQ paths
    } else {
      cli::cli_alert_warning("Direct FASTQ.gz download failed for {.val {run_id}}, attempting SRA download...")
      # Clean up any partially downloaded fastq files
      fs::file_delete(temp_files[fs::file_exists(temp_files)])
    }
  }


  # --- Fallback or Preference: Download SRA file ---
  sra_file_path <- fs::path(output_dir, run_id) # SRA files often named just by run_id
  downloaded_sra <- FALSE

  # Try ENA SRA link first if database_preference is 'ena'
  ena_sra_ftp_link <- NULL
  ena_sra_aspera_link <- NULL
  if (database_preference == "ena" && "submitted_ftp" %in% names(run_info)) { # Check ENA metadata format
    # Example: ftp.sra.ebi.ac.uk/vol1/srr/SRR600/004/SRR6000174
    # Find the link that looks like an SRA path, not lite version
    sra_links <- str_split(run_info$submitted_ftp, ";")[[1]] %>%
      discard(is.na) %>%
      discard(. == "") %>%
      str_subset("/srr/", negate = FALSE) %>% # Ensure it's an SRA path
      str_subset("lite.sra", negate = TRUE) # Exclude lite version if present
    if (length(sra_links) > 0) ena_sra_ftp_link <- sra_links[1] # Take the first match

    if (use_aspera && !is.null(ena_sra_ftp_link) && "submitted_aspera" %in% names(run_info)) {
      # Try to find corresponding aspera link
      aspera_sra_links <- str_split(run_info$submitted_aspera, ";")[[1]] %>%
        discard(is.na) %>%
        discard(. == "") %>%
        str_subset(":/vol1/srr/", negate = FALSE) %>% # Aspera path format
        str_subset("lite.sra", negate = TRUE)
      # Simple assumption: if one ftp link found, use first aspera link
      if (length(aspera_sra_links) > 0) ena_sra_aspera_link <- aspera_sra_links[1]
    }
  }

  if (!is.null(ena_sra_ftp_link)) {
    download_mode <- ifelse(use_aspera && !is.null(ena_sra_aspera_link), "Aspera", "FTP")
    cli::cli_alert_info("Attempting SRA download for {.val {run_id}} from ENA ({.val {download_mode}})...")
    link_to_use <- if (download_mode == "Aspera") ena_sra_aspera_link else ena_sra_ftp_link
    dest_file <- sra_file_path # Standard SRA file name

    tryCatch(
      {
        if (download_mode == "Aspera") {
          execute_aspera(
            remote_path = link_to_use, db = "ENA", software_paths = software_paths,
            aspera_key_path = aspera_key_path, output_dir = output_dir
          )
          # Aspera might download to current dir, ensure it's named correctly
          downloaded_basename <- fs::path_file(ena_sra_ftp_link) # Get basename from FTP link
          potential_downloaded_file <- fs::path(output_dir, downloaded_basename)
          # if(fs::file_exists(potential_downloaded_file) && !fs::path_equivalent(potential_downloaded_file, dest_file)) {
          if (fs::file_exists(potential_downloaded_file) && (potential_downloaded_file != dest_file)) {
            fs::file_move(potential_downloaded_file, dest_file)
          }
        } else {
          execute_download(
            url = paste0("https://", link_to_use), dest_file = dest_file,
            parallel_n = parallel_n, software_paths = software_paths
          )
        }
        if (fs::file_exists(dest_file) && fs::file_size(dest_file) > 0) {
          downloaded_sra <- TRUE
        } else {
          stop("Download resulted in missing or empty file.")
        }
      },
      error = function(e) {
        cli::cli_alert_warning("Failed to download SRA from ENA ({.val {download_mode}}) for {.val {run_id}}: {e$message}")
        if (fs::file_exists(dest_file)) fs::file_delete(dest_file) # Clean up
      }
    )
  }

  # If ENA failed or wasn't tried, use SRA (srapath -> https)
  if (!downloaded_sra) {
    cli::cli_alert_info("Attempting SRA download for {.val {run_id}} from SRA/NCBI (HTTPS)...")
    if (is.null(software_paths$srapath)) { # Note: srapath often part of sra-tools
      software_paths$srapath <- get_software_path("srapath", "sra-tools", required = TRUE) # Find srapath
      if (is.null(software_paths$srapath)) cli::cli_abort("Cannot find srapath.")
    }

    sra_https_link <- NULL
    tryCatch(
      {
        # srapath returns the https link(s)
        res <- run_command(software_paths$srapath, args = run_id, error_on_status = FALSE)
        if (res$status == 0 && nzchar(res$stdout)) {
          # Take the first link if multiple are returned
          sra_https_link <- str_split(str_trim(res$stdout), "\\n")[[1]][1]
        } else {
          # Check stderr for common 'not found' messages
          if (str_detect(res$stderr, "not found|resolve accession")) {
            stop("srapath could not find accession.")
          } else {
            stop(paste("srapath failed.", res$stderr))
          }
        }
      },
      error = function(e) {
        cli::cli_alert_warning("srapath failed for {.val {run_id}}: {e$message}")
      }
    )

    if (!is.null(sra_https_link)) {
      dest_file <- sra_file_path
      tryCatch(
        {
          cli::cli_inform("Downloading via HTTPS: {.url {sra_https_link}}")
          execute_download(
            url = sra_https_link, dest_file = dest_file,
            parallel_n = parallel_n, software_paths = software_paths
          )
          if (fs::file_exists(dest_file) && fs::file_size(dest_file) > 0) {
            downloaded_sra <- TRUE
          } else {
            stop("Download resulted in missing or empty file.")
          }
        },
        error = function(e) {
          cli::cli_alert_danger("Failed to download SRA from SRA/NCBI (HTTPS) for {.val {run_id}}: {e$message}")
          if (fs::file_exists(dest_file)) fs::file_delete(dest_file) # Clean up
        }
      )
    }
  }

  if (downloaded_sra) {
    downloaded_files$files <- sra_file_path
    downloaded_files$type <- "sra"
    cli::cli_alert_success("Successfully downloaded SRA file for {.val {run_id}} to {.path {sra_file_path}}.")
    return(downloaded_files)
  } else {
    cli::cli_abort("Failed to download any sequence data for run {.val {run_id}} from all sources.")
  }
}


#' Download a Single GSA Run
#'
#' Handles logic for downloading GSA files (FASTQ, BAM, etc.) using FTP/HTTPS/Aspera/Huawei.
#'
#' @param run_id CRR identifier.
#' @param run_info A single row tibble/data.frame from the metadata file for this run.
#' @param use_aspera Should Aspera be used if available?
#' @param parallel_n Parallel connections for axel.
#' @param output_dir Destination directory.
#' @param software_paths Paths to tools.
#' @return A list containing paths to downloaded files. Key: `files`. Type is implicit (various).
#' @keywords internal
download_gsa_run <- function(run_id, run_info, use_aspera = FALSE, parallel_n = 0,
                             output_dir = ".", software_paths) {
  downloaded_files <- list(files = character())
  aspera_key_path <- NULL
  if (use_aspera) {
    aspera_key_path <- get_aspera_key_path("GSA", software_paths)
    if (is.null(aspera_key_path)) {
      cli::cli_alert_warning("Aspera requested but key not found/configured for GSA, falling back to FTP/HTTPS.")
      use_aspera <- FALSE
    }
  }

  # GSA metadata CSV needs parsing to find download links.
  # The bash script scrapes the website. Let's try using the 'ftpFile' column if it exists in CSV.
  # GSA CSV Column names might be 'Run Accession', 'ftpFile', 'md5', 'fileSize', 'BioProject Accession' etc.
  # 'ftpFile' might contain multiple '|' separated links relative to a base FTP path.
  # This is less reliable than scraping but avoids complex scraping logic initially.

  # Alternative: Replicate scraping logic (more robust if website is stable)
  cra_id <- run_info$CRA_ID # Assuming this column exists and holds CRA
  if (is.null(cra_id)) {
    cli::cli_abort("Could not determine CRA ID from GSA metadata for run {.val {run_id}}.")
  }
  # https://download.cncb.ac.cn/gsa/CRA001160/CRR034520/
  run_page_url <- paste0("https://ngdc.cncb.ac.cn/gsa/browse/", cra_id, "/", run_id)
  cli::cli_alert_info("Scraping GSA run page for download links: {.url {run_page_url}}")

  all_links <- character()
  tryCatch(
    {
      resp <- httr::GET(run_page_url, httr::timeout(60))
      httr::stop_for_status(resp, task = paste("scrape GSA page for", run_id))
      content <- httr::content(resp, "text", encoding = "UTF-8")
      # Extract links matching patterns from bash script (adjust regex if needed)
      # Links like https://download.cncb.ac.cn/... or ftp://download.big.ac.cn/... or huaweicloud... ending in common suffixes
      links <- stringr::str_extract_all(content, '(https?://[^"]*\\.(gz|bam|tar|bz2|fastq|fq))|(ftp://[^"]*\\.(gz|bam|tar|bz2|fastq|fq))')[[1]]
      # Deduplicate and clean
      all_links <- unique(links) %>% discard(~ . == "")
      if (length(all_links) == 0) stop("No download links found on page.")
      cli::cli_alert_success("Found {length(all_links)} potential download link(s) for {.val {run_id}}.")
    },
    error = function(e) {
      cli::cli_abort("Failed to scrape GSA download links for {.val {run_id}}: {e$message}")
    }
  )

  # Prioritize links: Huawei > Aspera (if requested) > FTP > HTTPS
  huawei_links <- all_links[stringr::str_detect(all_links, "huaweicloud")]
  ftp_links <- all_links[stringr::str_detect(all_links, "^ftp://download\\.big\\.ac\\.cn")]
  https_links <- all_links[stringr::str_detect(all_links, "^https://download\\.cncb\\.ac\\.cn")]

  links_to_try <- list()
  download_mode <- NULL

  # TODO: Add default download mode, download_mode
  if (length(huawei_links) > 0) {
    links_to_try <- huawei_links
    download_mode <- "Huawei Cloud (HTTPS)"
    if (use_aspera) cli::cli_alert_info("Huawei Cloud links available, using them instead of Aspera for {.val {run_id}}.")
  } else if (use_aspera && length(ftp_links) > 0 && !is.null(aspera_key_path)) {
    # Convert FTP links to Aspera paths
    # ftp://download.big.ac.cn/ -> aspera01@download.cncb.ac.cn:
    links_to_try <- stringr::str_replace(ftp_links, "^ftp://download\\.big\\.ac\\.cn/", "aspera01@download.cncb.ac.cn:")
    download_mode <- "Aspera"
  } else if (length(ftp_links) > 0) {
    links_to_try <- ftp_links
    download_mode <- "FTP"
  } else if (length(https_links) > 0) {
    links_to_try <- https_links
    download_mode <- "HTTPS"
  } else {
    cli::cli_abort("No suitable download links found for {.val {run_id}} after filtering.")
  }

  cli::cli_alert_info("Attempting download for {.val {run_id}} using {.val {download_mode}} ({length(links_to_try)} file(s))...")

  success <- TRUE
  downloaded_run_files <- character()

  for (link in links_to_try) {
    # Determine destination filename
    dest_file <- fs::path(output_dir, fs::path_file(link))
    # If using Aspera, the source link is already the Aspera path, but basename should come from FTP/HTTPS equivalent
    if (download_mode == "Aspera") {
      # Find corresponding ftp link to get basename (this mapping assumes order is preserved)
      ftp_link_for_basename <- ftp_links[match(link, links_to_try)] # Match Aspera link back to FTP link list
      if (!is.na(ftp_link_for_basename)) {
        dest_file <- fs::path(output_dir, fs::path_file(ftp_link_for_basename))
      } else {
        # Fallback if matching failed - might create incorrect filename
        cli::cli_warn("Could not map Aspera link back to FTP link for basename, using Aspera path for filename.")
        # Try to infer from aspera path (e.g., after the colon)
        inferred_name <- tail(str_split(link, ":")[[1]], 1) %>% fs::path_file()
        if (nchar(inferred_name) > 0) dest_file <- fs::path(output_dir, inferred_name)
      }
    }

    cli::cli_inform("Downloading: {.url {link}} to {.path {dest_file}}")
    downloaded_run_files <- c(downloaded_run_files, dest_file)

    tryCatch(
      {
        if (download_mode == "Aspera") {
          execute_aspera(
            remote_path = link, db = "GSA", software_paths = software_paths,
            aspera_key_path = aspera_key_path, output_dir = output_dir
          )
        } else { # FTP, HTTPS, Huawei Cloud
          execute_download(
            url = link, dest_file = dest_file,
            parallel_n = parallel_n, software_paths = software_paths
          )
        }
        if (!fs::file_exists(dest_file) || fs::file_size(dest_file) == 0) {
          stop("Download resulted in missing or empty file.")
        }
      },
      error = function(e) {
        cli::cli_alert_danger("Failed to download {.url {link}}: {e$message}")
        success <<- FALSE
        if (fs::file_exists(dest_file)) fs::file_delete(dest_file) # Clean up failed part
      }
    )
    if (!success) break
  }

  if (success) {
    downloaded_files$files <- downloaded_run_files
    cli::cli_alert_success("Successfully downloaded file(s) for GSA run {.val {run_id}}.")
    return(downloaded_files)
  } else {
    cli::cli_alert_danger("Download failed for one or more files for GSA run {.val {run_id}}.")
    # Clean up any successful parts of this run? Or leave them? Let's leave them for potential resume.
    # fs::file_delete(downloaded_run_files[fs::file_exists(downloaded_run_files)])
    return(NULL) # Indicate failure
  }
}
