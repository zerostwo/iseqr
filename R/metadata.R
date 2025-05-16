# --- Metadata Functions (metadata.R) ---

#' Validate Accession and Generate ENA Query
#'
#' @param accession Input accession (Project, Study, Sample, Exp, Run, GSE, GSM).
#' @return ENA query string or throws error.
#' @keywords internal
validate_and_query_ena <- function(accession) {
  accession_upper <- toupper(accession)

  patterns <- list(
    project = "^PRJ[EDN][A-Z][0-9]+$",
    study = "^[EDS]RP[0-9]{6,}$",
    biosample = "^SAM[EDN][A-Z]?[0-9]+$",
    sample = "^[EDS]RS[0-9]{6,}$",
    experiment = "^[EDS]RX[0-9]{6,}$",
    run = "^[EDS]RR[0-9]{6,}$",
    geo_series = "^GSE[0-9]+$",
    geo_sample = "^GSM[0-9]+$"
  )

  type <- case_when(
    str_detect(accession_upper, patterns$project) ~ "project",
    str_detect(accession_upper, patterns$study) ~ "study",
    str_detect(accession_upper, patterns$biosample) ~ "biosample",
    str_detect(accession_upper, patterns$sample) ~ "sample",
    str_detect(accession_upper, patterns$experiment) ~ "experiment",
    str_detect(accession_upper, patterns$run) ~ "run",
    str_detect(accession_upper, patterns$geo_series) ~ "geo_series",
    str_detect(accession_upper, patterns$geo_sample) ~ "geo_sample",
    TRUE ~ "unknown"
  )

  if (type == "unknown") {
    cli::cli_abort("Invalid accession format: {.val {accession}}.")
  }

  # Handle GEO conversion
  if (type == "geo_series") {
    geo_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession)
    cli::cli_alert_info("Input is GEO Series {.val {accession}}, querying NCBI for BioProject...")
    tryCatch(
      {
        resp <- httr::GET(geo_url, httr::timeout(30))
        httr::stop_for_status(resp, task = paste("fetch BioProject for", accession))
        content <- httr::content(resp, "text", encoding = "UTF-8")
        # More robust extraction needed
        bioproject <- stringr::str_extract(content, "PRJ[EDN][A-Z][0-9]+") %>% unique()
        if (length(bioproject) == 0 || is.na(bioproject[1])) {
          cli::cli_abort("Could not find associated BioProject for GEO Series {.val {accession}}.")
        }
        bioproject <- bioproject[1] # Take the first one if multiple? Usually unique.
        cli::cli_alert_success("Found BioProject {.val {bioproject}} for GEO Series {.val {accession}}.")
        accession <- bioproject # Continue using the BioProject ID
        type <- "project"
      },
      error = function(e) {
        cli::cli_abort("Failed to query NCBI for GEO Series {.val {accession}}: {e$message}")
      }
    )
  } else if (type == "geo_sample") {
    geo_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession)
    cli::cli_alert_info("Input is GEO Sample {.val {accession}}, querying NCBI for BioSample...")
    tryCatch(
      {
        resp <- httr::GET(geo_url, httr::timeout(30))
        httr::stop_for_status(resp, task = paste("fetch BioSample for", accession))
        content <- httr::content(resp, "text", encoding = "UTF-8")
        # More robust extraction needed
        biosample <- stringr::str_extract(content, "SAM[EDN][A-Z]?[0-9]+") %>% unique()
        if (length(biosample) == 0 || is.na(biosample[1])) {
          cli::cli_abort("Could not find associated BioSample for GEO Sample {.val {accession}}.")
        }
        biosample <- biosample[1]
        cli::cli_alert_success("Found BioSample {.val {biosample}} for GEO Sample {.val {accession}}.")
        accession <- biosample # Continue using the BioSample ID
        type <- "biosample"
      },
      error = function(e) {
        cli::cli_abort("Failed to query NCBI for GEO Sample {.val {accession}}: {e$message}")
      }
    )
  }

  # Construct ENA query
  query <- switch(type,
    project = , # Fallthrough
    study = paste0("(study_accession=", accession, " OR secondary_study_accession=", accession, ")"),
    biosample = , # Fallthrough
    sample = paste0("(sample_accession=", accession, " OR secondary_sample_accession=", accession, ")"),
    experiment = paste0("experiment_accession=", accession),
    run = paste0("run_accession=", accession)
  )

  return(query)
}


#' Fetch Metadata from ENA/SRA
#'
#' @param accession The accession ID.
#' @param output_dir Directory to save metadata file.
#' @return Path to the metadata TSV file or NULL if failed.
#' @keywords internal
get_ena_sra_metadata <- function(accession, output_dir = ".") {
  metadata_file <- fs::path(output_dir, paste0(accession, ".metadata.tsv"))
  metadata_xml_file <- fs::path(output_dir, paste0(accession, ".xml")) # Temp file

  # 1. Try ENA Portal API
  cli::cli_alert_info("Fetching metadata for {.val {accession}} from ENA Portal API...")
  query <- tryCatch(validate_and_query_ena(accession), error = function(e) {
    cli::cli_warn("Failed validation or GEO lookup for {.val {accession}}: {e$message}")
    return(NULL) # Allow trying SRA if validation fails here? Or abort? Let's abort.
    # cli::cli_abort("Cannot proceed without valid accession/query.")
  })
  if (is.null(query)) cli::cli_abort("Cannot proceed without valid accession/query.")

  ena_url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&query=",
    URLencode(paste0("\"", query, "\"")), # Ensure query is URL encoded
    "&fields=all&limit=0" # Use limit=0 to get all results
  )

  metadata_found <- FALSE
  tryCatch(
    {
      # Use httr::write_disk for reliable download
      timeout_sec <- getOption("iseq.timeout", default = 240)
      resp <- httr::GET(ena_url, httr::write_disk(metadata_file, overwrite = TRUE), httr::progress(), httr::timeout(seconds = timeout_sec))
      httr::stop_for_status(resp, task = "fetch ENA metadata")

      # Check if file has content beyond header
      lines <- readr::read_lines(metadata_file, n_max = 2)
      if (length(lines) > 1 && nchar(lines[2]) > 0) {
        cli::cli_alert_success("Successfully downloaded ENA metadata to {.path {metadata_file}}")
        metadata_found <- TRUE
      } else {
        cli::cli_alert_warning("ENA metadata found for {.val {accession}} but appears empty (only header).")
        fs::file_delete(metadata_file) # Remove empty file
      }
    },
    error = function(e) {
      cli::cli_alert_warning("Failed to fetch metadata from ENA for {.val {accession}}: {e$message}")
      if (fs::file_exists(metadata_file)) fs::file_delete(metadata_file) # Clean up partial download
    }
  )

  # 2. Try SRA Entrez if ENA failed or was empty
  if (!metadata_found) {
    cli::cli_alert_info("ENA metadata not found or empty for {.val {accession}}, trying SRA Entrez...")
    search_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=", accession, "&usehistory=y")
    fetch_url_template <- "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&WebEnv=%s&query_key=%s"

    web_env <- NULL
    query_key <- NULL

    tryCatch(
      {
        resp_search <- httr::GET(search_url, httr::write_disk(metadata_xml_file, overwrite = TRUE), httr::timeout(60))
        httr::stop_for_status(resp_search, task = "SRA Entrez search")
        xml_content <- xml2::read_xml(metadata_xml_file)
        web_env <- xml2::xml_text(xml2::xml_find_first(xml_content, ".//WebEnv"))
        query_key <- xml2::xml_text(xml2::xml_find_first(xml_content, ".//QueryKey"))
        fs::file_delete(metadata_xml_file) # Clean up temp XML

        if (length(web_env) == 0 || length(query_key) == 0 || is.na(web_env) || is.na(query_key)) {
          stop("Could not extract WebEnv/QueryKey from SRA search result.")
        }

        fetch_url <- sprintf(fetch_url_template, web_env, query_key)
        resp_fetch <- httr::GET(fetch_url, httr::write_disk(metadata_file, overwrite = TRUE), httr::progress(), httr::timeout(120))
        httr::stop_for_status(resp_fetch, task = "SRA Entrez fetch runinfo")

        # SRA returns CSV, convert to TSV for consistency? Or keep as CSV?
        # The bash script converts. Let's convert.
        sra_csv <- readr::read_csv(metadata_file, show_col_types = FALSE)
        if (nrow(sra_csv) > 0) {
          # Be careful about quoting issues when converting CSV to TSV
          # Using readr::write_tsv might handle quoting better
          readr::write_tsv(sra_csv, metadata_file)
          cli::cli_alert_success("Successfully downloaded SRA metadata (and converted to TSV) to {.path {metadata_file}}")
          metadata_found <- TRUE
        } else {
          cli::cli_alert_warning("SRA metadata found for {.val {accession}} but appears empty.")
          fs::file_delete(metadata_file) # Remove empty file
        }
      },
      error = function(e) {
        cli::cli_alert_danger("Failed to fetch metadata from SRA for {.val {accession}}: {e$message}")
        if (fs::file_exists(metadata_xml_file)) fs::file_delete(metadata_xml_file)
        if (fs::file_exists(metadata_file)) fs::file_delete(metadata_file)
      }
    )
  }

  if (!metadata_found) {
    cli::cli_abort(c(
      "x" = "Failed to retrieve metadata for {.val {accession}} from both ENA and SRA.",
      "i" = "Check the accession ID and your internet connection."
    ))
  }

  return(metadata_file)
}


#' Fetch Metadata from GSA
#'
#' Handles GSA's specific metadata download process.
#' Note: GSA website structure/API might change. This requires careful implementation and testing.
#'
#' @param accession The GSA accession ID (CRA, CRR, CRX, PRJCA, SAMC).
#' @param output_dir Directory to save metadata file.
#' @return Path to the metadata CSV file or NULL if failed.
#' @keywords internal
get_gsa_metadata <- function(accession, output_dir = ".") {
  metadata_file_csv <- fs::path(output_dir, paste0(accession, ".metadata.csv"))
  metadata_file_xlsx <- fs::path(output_dir, paste0(accession, ".metadata.xlsx")) # For CRA/Project level
  temp_file <- fs::path_temp(paste0(accession, ".tmp"))
  on.exit(fs::file_delete(temp_file), add = TRUE)

  accession_upper <- toupper(accession)
  base_url <- "https://ngdc.cncb.ac.cn/gsa"
  runinfo_url <- ""
  form_data <- list()
  cra_id <- NULL # To store associated CRA for later XLSX download

  # Determine endpoint and parameters based on accession type
  if (stringr::str_detect(accession_upper, "^CRR[0-9]+$|^CRX[0-9]+$")) {
    cli::cli_alert_info("Accession {.val {accession}} is a Run/Experiment, finding associated CRA project...")
    search_url <- paste0(base_url, "/search?searchTerm=", accession)
    tryCatch(
      {
        resp_search <- httr::GET(search_url, httr::timeout(30))
        httr::stop_for_status(resp_search, task = "find CRA for GSA accession")
        content <- httr::content(resp_search, "text", encoding = "UTF-8")

        filtered_lines <- stringr::str_split(content, "\n")[[1]] %>%
          stringr::str_subset("example", negate = TRUE)

        # Extract CRA using regex - might need adjustment if website changes
        # cras <- stringr::str_extract_all(content, "CRA[0-9]+")[[1]] %>% unique()
        cras <- stringr::str_extract_all(filtered_lines, "CRA[0-9]+") %>%
          unlist() %>%
          unique()
        if (length(cras) == 0) stop("Could not find associated CRA.")
        if (length(cras) > 1) cli::cli_alert_warning("Multiple CRAs found ({.val {cras}}), using the first: {.val {cras[1]}}")
        cra_id <- cras[1]
        cli::cli_alert_success("Found associated project {.val {cra_id}} for {.val {accession}}.")

        # Now fetch run info for this specific CRA
        runinfo_url <- paste0(base_url, "/search/getRunInfoByCra")
        form_data <- list(searchTerm = cra_id, totalDatas = "9999", downLoadCount = "9999") # Assuming large enough counts
      },
      error = function(e) {
        cli::cli_abort("Failed to find associated CRA for GSA accession {.val {accession}}: {e$message}")
      }
    )
  } else if (stringr::str_detect(accession_upper, "^PRJC[A-Z][0-9]+$|^SAMC[0-9]+$")) {
    runinfo_url <- paste0(base_url, "/search/getRunInfo")
    # Encoding might be tricky here, check what the bash script implies with %26quot%3B
    # It seems to be URL-encoded quotes around the accession
    form_data <- list(searchTerm = paste0('"', accession, '"'), totalDatas = "9999", downLoadCount = "9999")
  } else if (stringr::str_detect(accession_upper, "^CRA[0-9]+$")) {
    cra_id <- accession # Store CRA for later XLSX download
    runinfo_url <- paste0(base_url, "/search/getRunInfoByCra")
    form_data <- list(searchTerm = accession, totalDatas = "9999", downLoadCount = "9999")
  } else {
    cli::cli_abort("Invalid GSA accession format: {.val {accession}}.")
  }

  # Fetch the main run info (usually CSV-like text)
  metadata_found_csv <- FALSE
  cli::cli_alert_info("Fetching GSA run metadata for {.val {accession}}...")
  tryCatch(
    {
      resp_runinfo <- httr::POST(runinfo_url,
        body = form_data, encode = "form",
        httr::write_disk(temp_file, overwrite = TRUE), httr::progress(), httr::timeout(120)
      )
      httr::stop_for_status(resp_runinfo, task = "fetch GSA runinfo")

      # Check content, GSA might return HTML on error or empty results in various ways
      # A simple check is file size or reading the first few lines
      if (fs::file_size(temp_file) > 100) { # Arbitrary small size check
        # If it was a CRR/CRX, filter only the relevant run/exp
        if (stringr::str_detect(accession_upper, "^CRR[0-9]+$|^CRX[0-9]+$")) {
          # GSA response might be comma-separated text, read and filter
          # Need to know the column structure - assuming run/exp accession is in one of the first columns
          all_data <- readr::read_delim(temp_file, delim = ",", show_col_types = FALSE, guess_max = 5000) # guess_max important
          # Find the column containing the run/exp ID (adjust column name/index as needed)
          # Let's assume it might be 'Run Accession' or 'Experiment Accession'
          # run_col <- names(all_data)[stringr::str_detect(names(all_data), "Run\\s+Accession")]
          # exp_col <- names(all_data)[stringr::str_detect(names(all_data), "Experiment\\s+Accession")]
          run_col <- names(all_data)[stringr::str_detect(names(all_data), "Run")]
          exp_col <- names(all_data)[stringr::str_detect(names(all_data), "Experiment")]
          relevant_data <- all_data %>%
            filter(if (length(run_col) > 0) {
              .data[[run_col]] == accession_upper
            } else {
              FALSE |
                if (length(exp_col) > 0) .data[[exp_col]] == accession_upper else FALSE
            })

          if (nrow(relevant_data) > 0) {
            readr::write_csv(relevant_data, metadata_file_csv)
            metadata_found_csv <- TRUE
          } else {
            cli::cli_warn("Could not find specific run/experiment {.val {accession}} within the fetched CRA data.")
          }
        } else {
          # For Project/Study/Sample, copy the whole file
          fs::file_copy(temp_file, metadata_file_csv, overwrite = TRUE)
          metadata_found_csv <- TRUE
        }

        if (metadata_found_csv) {
          # Verify it's not empty after potential filtering
          if (!check_file_valid(metadata_file_csv, accession)) {
            metadata_found_csv <- FALSE
            fs::file_delete(metadata_file_csv)
          } else {
            cli::cli_alert_success("Successfully downloaded GSA metadata CSV to {.path {metadata_file_csv}}")
          }
        }
      } else {
        cli::cli_alert_warning("Downloaded GSA metadata for {.val {accession}} appears empty or invalid.")
      }
    },
    error = function(e) {
      cli::cli_alert_danger("Failed to fetch GSA metadata for {.val {accession}}: {e$message}")
    }
  )

  if (!metadata_found_csv) {
    cli::cli_abort("Failed to retrieve valid GSA metadata CSV for {.val {accession}}.")
  }

  # Download associated XLSX metadata if applicable (CRA level)
  # Extract unique CRAs from the downloaded CSV
  tryCatch(
    {
      gsa_csv_data <- readr::read_csv(metadata_file_csv, show_col_types = FALSE)
      # Assuming CRA is in a column named 'BioProject Accession' or similar? Adjust based on actual GSA output.
      # Let's guess a column name. This is fragile.
      if ("BioProject" %in% names(gsa_csv_data)) {
        cra_ids_in_csv <- unique(gsa_csv_data$`BioProject`)
        cra_ids_in_csv <- cra_ids_in_csv[stringr::str_detect(cra_ids_in_csv, "^CRA[0-9]+$")]
      } else if (!is.null(cra_id)) {
        cra_ids_in_csv <- cra_id # Use the one identified earlier if column not found
      } else {
        cra_ids_in_csv <- character(0)
      }

      if (length(cra_ids_in_csv) == 0 && !is.null(cra_id)) { # Fallback if extraction failed but we know the CRA
        cra_ids_in_csv <- cra_id
      }

      for (current_cra in cra_ids_in_csv) {
        xlsx_url <- paste0(base_url, "/file/exportExcelFile")
        xlsx_file_path <- fs::path(output_dir, paste0(current_cra, ".metadata.xlsx"))
        if (!fs::file_exists(xlsx_file_path)) {
          cli::cli_alert_info("Downloading associated GSA metadata XLSX for {.val {current_cra}}...")
          # This POST request structure might need adjustment
          resp_xlsx <- httr::POST(xlsx_url,
            body = list(type = "3", dlAcession = current_cra),
            encode = "form",
            httr::write_disk(xlsx_file_path, overwrite = TRUE),
            httr::progress(), httr::timeout(120)
          )
          # GSA might return non-200 even on success with file download, check file existence/size
          if (fs::file_exists(xlsx_file_path) && fs::file_size(xlsx_file_path) > 0) {
            cli::cli_alert_success("Successfully downloaded GSA metadata XLSX to {.path {xlsx_file_path}}")
          } else {
            cli::cli_alert_warning("Failed to download or validate GSA metadata XLSX for {.val {current_cra}} (Status: {httr::status_code(resp_xlsx)}).")
            if (fs::file_exists(xlsx_file_path)) fs::file_delete(xlsx_file_path) # Clean up empty/failed download
          }
        } else {
          cli::cli_alert_info("GSA metadata XLSX already exists: {.path {xlsx_file_path}}")
        }
      }
    },
    error = function(e) {
      cli::cli_warn("Problem fetching/processing associated GSA XLSX metadata: {e$message}")
    }
  )


  return(metadata_file_csv) # Return path to the primary CSV metadata
}


#' Get Metadata Wrapper
#'
#' Determines the source (GSA or ENA/SRA) and calls the appropriate function.
#'
#' @param accession Accession ID.
#' @param output_dir Output directory.
#' @return Path to the metadata file.
#' @keywords internal
get_metadata <- function(accession, output_dir = ".") {
  accession_upper <- toupper(accession)
  metadata_file <- NULL

  # Check if metadata already exists
  gsa_meta_csv <- fs::path(output_dir, paste0(accession, ".metadata.csv"))
  ena_meta_tsv <- fs::path(output_dir, paste0(accession, ".metadata.tsv"))

  if (stringr::str_detect(accession_upper, "^(PRJC|SAMC|CRA|CRX|CRR)")) {
    if (fs::file_exists(gsa_meta_csv)) {
      cli::cli_alert_info("GSA metadata file exists: {.path {gsa_meta_csv}}")
      if (!check_file_valid(gsa_meta_csv, accession)) {
        fs::file_delete(gsa_meta_csv)
        cli::cli_alert_warning("Existing GSA metadata file was empty, attempting re-download.")
        metadata_file <- get_gsa_metadata(accession, output_dir)
      } else {
        metadata_file <- gsa_meta_csv
      }
    } else {
      metadata_file <- get_gsa_metadata(accession, output_dir)
    }
    attr(metadata_file, "source") <- "GSA"
  } else { # Assume ENA/SRA
    if (fs::file_exists(ena_meta_tsv)) {
      cli::cli_alert_info("ENA/SRA metadata file exists: {.path {ena_meta_tsv}}")
      if (!check_file_valid(ena_meta_tsv, accession)) {
        fs::file_delete(ena_meta_tsv)
        cli::cli_alert_warning("Existing ENA/SRA metadata file was empty, attempting re-download.")
        metadata_file <- get_ena_sra_metadata(accession, output_dir)
      } else {
        metadata_file <- ena_meta_tsv
      }
    } else {
      metadata_file <- get_ena_sra_metadata(accession, output_dir)
    }
    attr(metadata_file, "source") <- "ENA/SRA"
  }

  # Final check after potential download
  if (is.null(metadata_file) || !fs::file_exists(metadata_file) || fs::file_size(metadata_file) == 0) {
    cli::cli_abort("Failed to obtain valid metadata for {.val {accession}}.")
  }

  return(metadata_file)
}
