# --- Main Function (iseq_download.R) ---

#' Download Sequencing Data and Metadata
#'
#' Downloads data for Runs associated with a given accession (Project, Study,
#' Sample, Experiment, Run, GEO Series/Sample) from ENA/SRA or GSA databases.
#'
#' @param input Character vector of accession IDs, or a path to a file containing
#'   one accession per line.
#' @param output_dir Path to the output directory. Will be created if it doesn't exist.
#'   Defaults to the current working directory.
#' @param metadata_only Logical. If TRUE, only download metadata files, skip sequence data. (Default: FALSE)
#' @param gzip Logical. If TRUE, prefer downloading gzipped FASTQ files directly.
#'   If not available or if SRA files are downloaded, SRA will be converted and
#'   the resulting FASTQ files will be compressed using `pigz`. Requires `pigz`. (Default: FALSE)
#' @param fastq Logical. If TRUE, convert downloaded SRA files to FASTQ format using
#'   `fasterq-dump`. Requires `sra-tools`. (Default: FALSE)
#' @param merge Character. Merge FASTQ files based on level: "ex" (Experiment),
#'   "sa" (Sample), "st" (Study), or "none". Requires `sra-tools` if conversion needed. (Default: "none")
#' @param threads Integer. Number of threads for `fasterq-dump` and `pigz`. (Default: 8)
#' @param database Character. Preferred database for SRA/Run download: "ena" or "sra".
#'   Note: Some data might only be available from one source. "auto" is not implemented here yet, defaults to ENA preference. (Default: "ena")
#' @param parallel Integer. Number of parallel connections for download using `axel`.
#'   Requires `axel`. If 0, uses `wget`. (Default: 0)
#' @param aspera Logical. Use Aspera (`ascp`) for download if available (ENA/GSA FTP sources only).
#'   Requires `aspera-cli` and configured keys (see Details). (Default: FALSE)
#' @param retries Integer. Number of times to retry download/validation on failure. (Default: 2)
#' @param remove_sra Logical. If TRUE and `fastq` or `gzip` is TRUE, remove the original SRA file
#'   after successful conversion/compression. (Default: TRUE)
#'
#' @details
#'   **Software Requirements:** This function relies on external command-line tools.
#'   Ensure the following are installed and accessible in your system's PATH:
#'   - `wget` (required for basic download)
#'   - `md5sum` (from `coreutils`, for validation)
#'   Optional, based on arguments:
#'   - `axel` (if `parallel > 0`)
#'   - `pigz` (if `gzip = TRUE`)
#'   - `sra-tools` (>= 2.11.0 recommended, provides `fasterq-dump`, `vdb-validate`, `srapath`) (if `fastq=TRUE` or SRA download occurs)
#'   - `aspera-cli` (provides `ascp`) (if `aspera = TRUE`)
#'
#'   You can provide paths to these tools via R's `options()` if they are not in the PATH, e.g.,
#'   `options(iseq.wget = "/path/to/wget", iseq.pigz = "/path/to/pigz")`.
#'
#'   **Aspera Keys:** For Aspera downloads:
#'   - ENA: The key file (`aspera_bypass_rsa.pem` or `aspera_tokenauth_id_rsa`) is usually found
#'     within the `aspera-cli` installation (e.g., `.../etc/`). The function attempts to find it.
#'     Alternatively, set `options(iseq.aspera_key_ena = "/path/to/ena_key.pem")`.
#'   - GSA: Requires a specific key file (`.asperaGSA.openssh`). The function will attempt to download it
#'     from the GSA website to the `output_dir` if not found. Alternatively, place it there manually
#'     or set `options(iseq.aspera_key_gsa = "/path/to/gsa_key.openssh")`.
#'
#'   **Metadata:** Metadata files (`<accession>.metadata.tsv` for ENA/SRA, `<accession>.metadata.csv`
#'   and potentially `<CRA>.metadata.xlsx` for GSA) are saved in the `output_dir`.
#'
#'   **Logging:** Download success and failures for individual runs are logged to `success.log` and `fail.log`
#'   in the `output_dir`. Existing entries in `success.log` will cause runs to be skipped unless the log file is removed/edited.
#'
#'   **Merging:** Merging concatenates FASTQ files. For GSA data, merging assumes FASTQ(.gz) files exist and follows naming conventions inferred during download; BAM or other file types are not merged.
#'
#' @return Invisibly returns a list containing paths to successfully processed top-level accessions' metadata files.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage (ensure required software is installed)
#'
#' # Set paths if not in PATH (example)
#' options(iseq.pigz = "/opt/homebrew/bin/pigz")
#' options(iseq.fasterq_dump = "/path/to/fasterq-dump")
#'
#' # Download metadata only for a project
#' iseq_download("PRJNA12345", metadata_only = TRUE, output_dir = "project_meta")
#'
#' # Download SRA, convert to FASTQ, compress, use 4 threads
#' iseq_download("SRR67890", fastq = TRUE, gzip = TRUE, threads = 4, output_dir = "srr_data")
#'
#' # Download project data using Aspera, 10 parallel connections (if Aspera fails),
#' # prefer gzipped FASTQ, merge by sample
#' iseq_download("PRJEB9876",
#'   output_dir = "project_data",
#'   gzip = TRUE, aspera = TRUE, parallel = 10,
#'   merge = "sa", threads = 8
#' )
#'
#' # Download from a file containing multiple accessions
#' writeLines(c("SRR1000", "SRR1001"), "accessions.txt")
#' iseq_download("accessions.txt", output_dir = "multi_srr", fastq = TRUE)
#' unlink("accessions.txt") # Clean up example file
#' }
# --- Helper Functions ---

#' Validate and initialize parameters
#' @keywords internal
validate_and_init_params <- function(input, output_dir, merge, database, threads, parallel, retries) {
  # Validate numeric parameters
  if (!is.numeric(threads) || threads < 1) {
    cli::cli_abort("{.arg threads} must be a positive integer.")
  }
  threads <- as.integer(threads)

  if (!is.numeric(parallel) || parallel < 0) {
    cli::cli_abort("{.arg parallel} must be a non-negative integer.")
  }
  parallel <- as.integer(parallel)

  if (!is.numeric(retries) || retries < 0) {
    cli::cli_abort("{.arg retries} must be a non-negative integer.")
  }
  retries <- as.integer(retries)

  # Process input: file or vector
  if (length(input) == 1 && fs::is_file(input)) {
    accessions <- readr::read_lines(input) %>%
      stringr::str_trim() %>%
      discard(~ . == "")
    cli::cli_alert_info("Read {length(accessions)} accessions from file: {.path {input}}")
  } else {
    accessions <- as.character(input) %>%
      stringr::str_trim() %>%
      discard(~ . == "")
  }
  if (length(accessions) == 0) {
    cli::cli_abort("No valid accessions provided in {.arg input}.")
  }

  # Setup output directory
  if (!fs::dir_exists(output_dir)) {
    cli::cli_alert_info("Creating output directory: {.path {output_dir}}")
    fs::dir_create(output_dir, recurse = TRUE)
  }
  output_dir <- fs::path_real(output_dir)

  return(list(
    accessions = accessions,
    output_dir = output_dir,
    merge = merge,
    database = database,
    threads = threads,
    parallel = parallel,
    retries = retries
  ))
}

#' Check and initialize required software paths
#' @keywords internal
init_software_paths <- function(gzip, fastq, merge, parallel, aspera) {
  software_paths <- list()

  # Required tools
  software_paths$wget <- get_software_path("wget", required = (parallel == 0))
  software_paths$md5sum <- get_software_path("md5sum", "coreutils", required = TRUE)

  # Optional tools based on parameters
  if (parallel > 0) software_paths$axel <- get_software_path("axel", required = TRUE)
  if (gzip) software_paths$pigz <- get_software_path("pigz", required = TRUE)
  if (fastq || merge != "none" || gzip) {
    software_paths$fasterq_dump <- get_software_path("fasterq_dump", "sra-tools>=2.11.0", required = TRUE)
    software_paths$vdb_validate <- get_software_path("vdb_validate", "sra-tools>=2.11.0", required = TRUE)
  }
  if (aspera) {
    software_paths$ascp <- get_software_path("ascp", "aspera-cli=4.14.0", required = TRUE)
  }

  cli::cli_alert_success("Required software checks passed.")
  return(software_paths)
}

#' Process metadata for a single accession
#' @keywords internal
process_metadata <- function(accession, output_dir = NULL) {
  tryCatch(
    {
      metadata_file <- get_metadata(accession, output_dir)
      metadata_source <- attr(metadata_file, "source")
      # cli::cli_alert_success("Metadata ready: {.path {metadata_file}} (Source: {metadata_source})")

      # Parse metadata to get run list
      if (metadata_source == "GSA") {
        run_col_gsa <- "Run"
        meta_df <- readr::read_csv(metadata_file, show_col_types = FALSE)
        if (!run_col_gsa %in% names(meta_df)) {
          cli::cli_abort("Cannot find run column {.val {run_col_gsa}} in GSA metadata.")
        }
        run_list_df <- meta_df %>% distinct(across(all_of(run_col_gsa)), .keep_all = TRUE)
        run_list_df$CRA_ID <- accession
      } else { # ENA/SRA
        run_col_ena <- "run_accession"
        meta_df <- readr::read_tsv(metadata_file, show_col_types = FALSE)
        if (!run_col_ena %in% names(meta_df)) {
          cli::cli_abort("Cannot find run column {.val {run_col_ena}} in ENA/SRA metadata.")
        }
        run_list_df <- meta_df %>% distinct(across(all_of(run_col_ena)), .keep_all = TRUE)
      }

      if (nrow(run_list_df) == 0) {
        cli::cli_alert_warning("No runs found in metadata for accession {.val {accession}}.")
        return(NULL)
      }

      cli::cli_alert_info("Found {nrow(run_list_df)} unique run(s) associated with {.val {accession}}.")
      return(list(
        metadata_file = metadata_file,
        metadata_source = metadata_source,
        run_list_df = run_list_df
      ))
    },
    error = function(e) {
      cli::cli_alert_danger("Failed to get or parse metadata for {.val {accession}}: {e$message}")
      return(NULL)
    }
  )
}

#' Process a single run with retries
#' @keywords internal
process_run <- function(run_id, run_info, metadata_source, params, software_paths) {
  current_try <- 0
  download_successful <- FALSE
  validation_successful <- FALSE
  downloaded_info <- NULL

  while (current_try <= params$retries && !validation_successful) {
    if (current_try > 0) {
      cli::cli_alert_info("Retrying run {.val {run_id}} (Attempt {current_try}/{params$retries})...")
      if (!is.null(downloaded_info) && length(downloaded_info$files) > 0) {
        cli::cli_inform("Removing potentially failed files from previous attempt: {.path {downloaded_info$files}}")
        fs::file_delete(downloaded_info$files[fs::file_exists(downloaded_info$files)])
        downloaded_info <- NULL
      }
    }

    # Download
    tryCatch(
      {
        if (metadata_source == "GSA") {
          downloaded_info <- download_gsa_run(run_id, run_info,
            use_aspera = params$aspera,
            parallel_n = params$parallel,
            output_dir = params$output_dir,
            software_paths = software_paths
          )
        } else {
          downloaded_info <- download_ena_sra_run(
            run_id,
            run_info,
            download_fastq_gz = params$gzip,
            use_aspera = params$aspera,
            database_preference = params$database,
            parallel_n = params$parallel,
            output_dir = params$output_dir,
            software_paths = software_paths
          )
        }
        download_successful <- !is.null(downloaded_info) && length(downloaded_info$files) > 0
      },
      error = function(e) {
        cli::cli_alert_danger("Download attempt failed for run {.val {run_id}}: {e$message}")
        download_successful <- FALSE
      }
    )

    if (!download_successful) {
      current_try <- current_try + 1
      next
    }

    # Validate
    cli::cli_alert_info("Validating downloaded file(s) for run {.val {run_id}}...")
    tryCatch(
      {
        if (metadata_source == "GSA") {
          validation_successful <- validate_gsa_run(run_id, downloaded_info, run_info, params$output_dir, software_paths)
        } else {
          validation_successful <- validate_ena_sra_run(run_id, downloaded_info, run_info, software_paths)
        }
      },
      error = function(e) {
        cli::cli_alert_danger("Validation failed for run {.val {run_id}}: {e$message}")
        validation_successful <- FALSE
      }
    )

    if (validation_successful) {
      cli::cli_alert_success("Validation successful for run {.val {run_id}}.")
      break
    } else {
      cli::cli_alert_warning("Validation failed for run {.val {run_id}}.")
      current_try <- current_try + 1
    }
  }

  if (!validation_successful) {
    cli::cli_alert_danger("Run {.val {run_id}} failed download/validation after {params$retries} retries.")
    log_status(run_id, "fail", params$output_dir)
    if (!is.null(downloaded_info) && length(downloaded_info$files) > 0) {
      cli::cli_inform("Removing failed files: {.path {downloaded_info$files}}")
      fs::file_delete(downloaded_info$files[fs::file_exists(downloaded_info$files)])
    }
    return(NULL)
  }

  return(downloaded_info)
}

#' Post-process downloaded files
#' @keywords internal
post_process_files <- function(run_id, downloaded_info, params, software_paths) {
  current_run_files <- downloaded_info$files
  is_sra_file <- !is.null(downloaded_info$type) && downloaded_info$type == "sra"

  # Convert SRA to FASTQ if needed
  if (is_sra_file && (params$fastq || params$gzip || params$merge != "none")) {
    cli::cli_alert_info("Post-processing: Converting SRA to FASTQ for {.val {run_id}}...")
    tryCatch(
      {
        converted_fastq_files <- convert_sra_to_fastq(
          sra_file = current_run_files,
          output_dir = params$output_dir,
          threads = params$threads,
          software_paths = software_paths
        )
        if (length(converted_fastq_files) > 0) {
          if (params$remove_sra) {
            cli::cli_inform("Removing original SRA file: {.path {current_run_files}}")
            fs::file_delete(current_run_files)
          }
          current_run_files <- converted_fastq_files
          is_sra_file <- FALSE
        } else {
          cli::cli_alert_warning("SRA to FASTQ conversion failed or yielded no files for {.val {run_id}}.")
          log_status(run_id, "fail", params$output_dir)
          return(NULL)
        }
      },
      error = function(e) {
        cli::cli_alert_danger("SRA to FASTQ conversion failed for {.val {run_id}}: {e$message}")
        log_status(run_id, "fail", params$output_dir)
        return(NULL)
      }
    )
  }

  # Compress FASTQ if needed
  needs_compression <- params$gzip && !is_sra_file && !all(stringr::str_ends(current_run_files, ".gz"))
  if (needs_compression) {
    cli::cli_alert_info("Post-processing: Compressing FASTQ files for {.val {run_id}}...")
    tryCatch(
      {
        compressed_files <- compress_fastq(
          fastq_files = current_run_files,
          threads = params$threads,
          software_paths = software_paths
        )
        if (length(compressed_files) > 0) {
          current_run_files <- compressed_files
        } else {
          cli::cli_alert_warning("FASTQ compression failed or yielded no files for {.val {run_id}}.")
        }
      },
      error = function(e) {
        cli::cli_alert_danger("FASTQ compression failed for {.val {run_id}}: {e$message}")
      }
    )
  }

  log_status(run_id, "success", params$output_dir)
  return(current_run_files)
}

#' Main function
#' @export
iseq_download <- function(input,
                          output_dir = ".",
                          metadata_only = FALSE,
                          gzip = FALSE,
                          fastq = FALSE,
                          merge = c("none", "ex", "sa", "st"),
                          threads = 8,
                          database = c("ena", "sra"),
                          parallel = 0,
                          aspera = FALSE,
                          retries = 2,
                          remove_sra = TRUE) {
  # Validate merge parameter
  merge <- rlang::arg_match(merge)
  database <- rlang::arg_match(database)

  # 1. Validate and initialize parameters
  params <- validate_and_init_params(input, output_dir, merge, database, threads, parallel, retries)
  params$metadata_only <- metadata_only
  params$gzip <- gzip
  params$fastq <- fastq
  params$aspera <- aspera
  params$remove_sra <- remove_sra

  # 2. Check software dependencies
  software_paths <- init_software_paths(gzip, fastq, merge, parallel, aspera)

  # 3. Process each accession
  processed_metadata_files <- list()

  for (accession in params$accessions) {
    cli::cli_h1("Processing Accession: {.val {accession}}")
    accession_start_time <- Sys.time()

    # Get and process metadata
    metadata_result <- process_metadata(accession, params$output_dir)
    if (is.null(metadata_result)) next

    processed_metadata_files[[accession]] <- metadata_result$metadata_file

    # Skip download if requested
    if (params$metadata_only) {
      cli::cli_alert_info("Metadata only requested, skipping sequence download for {.val {accession}}.")
      next
    }

    # Process runs
    runs_to_process <- if (metadata_result$metadata_source == "GSA") {
      metadata_result$run_list_df$Run
    } else {
      metadata_result$run_list_df$run_accession
    }

    all_runs_successful <- TRUE
    final_run_files <- list()

    for (run_index in seq_along(runs_to_process)) {
      run_id <- runs_to_process[run_index]
      run_info <- metadata_result$run_list_df[run_index, ]

      cli::cli_h2("Processing Run: {.val {run_id}} ({run_index}/{length(runs_to_process)})")

      # Skip if already successful
      if (is_already_successful(run_id, params$output_dir)) {
        cli::cli_alert_success("Run {.val {run_id}} already listed in success.log, skipping.")
        final_run_files[[run_id]] <- fs::dir_ls(params$output_dir, glob = paste0(run_id, "*"), recurse = FALSE)
        next
      }

      # Process run with retries
      downloaded_info <- process_run(
        run_id = run_id,
        run_info = run_info,
        metadata_source = metadata_result$metadata_source,
        params = params,
        software_paths = software_paths
      )
      if (is.null(downloaded_info)) {
        all_runs_successful <- FALSE
        next
      }

      # Post-process files
      final_files <- post_process_files(run_id, downloaded_info, params, software_paths)
      if (!is.null(final_files)) {
        final_run_files[[run_id]] <- final_files
      } else {
        all_runs_successful <- FALSE
      }
    }

    # Merge files if requested
    if (params$merge != "none" && all_runs_successful && length(final_run_files) > 0) {
      cli::cli_h2("Merging files for accession {.val {accession}}")
      first_run_id <- names(final_run_files)[1]
      are_files_gzipped <- any(stringr::str_ends(final_run_files[[first_run_id]], ".gz"))

      tryCatch(
        {
          merge_fastq_files(
            accession = accession,
            metadata_file = metadata_result$metadata_file,
            source = metadata_result$metadata_source,
            merge_level = params$merge,
            output_dir = params$output_dir,
            is_gzipped = are_files_gzipped,
            remove_originals = TRUE
          )
        },
        error = function(e) {
          cli::cli_alert_danger("Merging failed for accession {.val {accession}}: {e$message}")
        }
      )
    } else if (params$merge != "none" && !all_runs_successful) {
      cli::cli_alert_warning("Skipping merge for {.val {accession}} because some runs failed.")
    }

    accession_end_time <- Sys.time()
    accession_duration <- difftime(accession_end_time, accession_start_time)
    cli::cli_alert_success("Finished processing accession {.val {accession}} in {format(accession_duration)}")
    cli::cli_rule()
  }

  cli::cli_h1("All processing complete.")
  invisible(processed_metadata_files)
}
