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
#' @param outdir Directory to save metadata file. If NULL, only return data.frame.
#' @return Data frame of metadata, and optionally write TSV file.
#' @keywords internal
get_ena_sra_metadata <- function(accession, outdir = NULL) {
  metadata_file <- if (!is.null(outdir)) fs::path(outdir, paste0(accession, ".metadata.tsv")) else tempfile(fileext = ".tsv")
  metadata_xml_file <- withr::local_tempfile(fileext = ".xml")  # auto-cleaned up

  cli::cli_alert_info("Fetching metadata for {.val {accession}} from ENA Portal API...")
  query <- tryCatch(validate_and_query_ena(accession), error = function(e) {
    cli::cli_warn("Failed validation or GEO lookup for {.val {accession}}: {e$message}")
    return(NULL)
  })
  if (is.null(query)) cli::cli_abort("Cannot proceed without valid accession/query.")

  ena_url <- paste0(
    "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&query=",
    URLencode(paste0("\"", query, "\"")),
    "&fields=all&limit=0"
  )

  metadata_found <- FALSE
  metadata_df <- NULL
  tryCatch({
    timeout_sec <- getOption("iseq.timeout", default = 240)
    resp <- httr::GET(ena_url, httr::write_disk(metadata_file, overwrite = TRUE), httr::progress(), httr::timeout(seconds = timeout_sec))
    httr::stop_for_status(resp, task = "fetch ENA metadata")

    lines <- readr::read_lines(metadata_file, n_max = 2)
    if (length(lines) > 1 && nchar(lines[2]) > 0) {
      metadata_df <- readr::read_tsv(metadata_file, show_col_types = FALSE)
      cli::cli_alert_success("Successfully downloaded ENA metadata for {.val {accession}}")
      metadata_found <- TRUE
    } else {
      cli::cli_alert_warning("ENA metadata for {.val {accession}} appears empty.")
    }
  }, error = function(e) {
    cli::cli_alert_warning("Failed to fetch metadata from ENA: {e$message}")
  })

  if (!metadata_found) {
    cli::cli_alert_info("Trying SRA Entrez for {.val {accession}}...")
    search_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=", accession, "&usehistory=y")
    fetch_url_template <- "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&WebEnv=%s&query_key=%s"

    tryCatch({
      resp_search <- httr::GET(search_url, httr::write_disk(metadata_xml_file, overwrite = TRUE), httr::timeout(60))
      httr::stop_for_status(resp_search, task = "SRA Entrez search")
      xml <- xml2::read_xml(metadata_xml_file)
      web_env <- xml2::xml_text(xml2::xml_find_first(xml, ".//WebEnv"))
      query_key <- xml2::xml_text(xml2::xml_find_first(xml, ".//QueryKey"))
      fetch_url <- sprintf(fetch_url_template, web_env, query_key)

      resp_fetch <- httr::GET(fetch_url, httr::write_disk(metadata_file, overwrite = TRUE), httr::progress(), httr::timeout(120))
      httr::stop_for_status(resp_fetch, task = "SRA Entrez fetch runinfo")

      sra_csv <- readr::read_csv(metadata_file, show_col_types = FALSE)
      if (nrow(sra_csv) > 0) {
        metadata_df <- sra_csv
        cli::cli_alert_success("Successfully downloaded and parsed SRA metadata for {.val {accession}}")
        metadata_found <- TRUE
      } else {
        cli::cli_alert_warning("SRA metadata for {.val {accession}} appears empty.")
      }
    }, error = function(e) {
      cli::cli_alert_danger("Failed to fetch metadata from SRA for {.val {accession}}: {e$message}")
    })
  }

  if (!metadata_found) cli::cli_abort("Failed to retrieve metadata for {.val {accession}} from both ENA and SRA.")

  if (!is.null(outdir)) readr::write_tsv(metadata_df, metadata_file)
  return(metadata_df)
}

#' Fetch Metadata from GSA
#'
#' Handles GSA's specific metadata download process.
#' Note: GSA website structure/API might change. This requires careful implementation and testing.
#'
#' @param accession The GSA accession ID (CRA, CRR, CRX, PRJCA, SAMC).
#' @param outdir Directory to save metadata file. If NULL, only return data.frame.
#' @return Data frame of metadata and optionally write to CSV/XLSX.
#' @keywords internal
get_gsa_metadata <- function(accession, outdir = NULL) {
  metadata_file_csv <- if (!is.null(outdir)) fs::path(outdir, paste0(accession, ".metadata.csv")) else withr::local_tempfile(fileext = ".csv")
  metadata_file_xlsx <- if (!is.null(outdir)) fs::path(outdir, paste0(accession, ".metadata.xlsx")) else withr::local_tempfile(fileext = ".xlsx")
  temp_file <- withr::local_tempfile(fileext = ".tmp")

  accession_upper <- toupper(accession)
  base_url <- "https://ngdc.cncb.ac.cn/gsa"
  runinfo_url <- ""
  form_data <- list()
  cra_id <- NULL

  if (stringr::str_detect(accession_upper, "^CRR[0-9]+$|^CRX[0-9]+$")) {
    cli::cli_alert_info("Accession {.val {accession}} is a Run/Experiment, finding associated CRA project...")
    search_url <- paste0(base_url, "/search?searchTerm=", accession)
    tryCatch({
      resp_search <- httr::GET(search_url, httr::timeout(30))
      httr::stop_for_status(resp_search, task = "find CRA for GSA accession")
      content <- httr::content(resp_search, "text", encoding = "UTF-8")
      filtered_lines <- stringr::str_split(content, "\n")[[1]] %>% stringr::str_subset("example", negate = TRUE)
      cras <- stringr::str_extract_all(filtered_lines, "CRA[0-9]+") %>% unlist() %>% unique()
      if (length(cras) == 0) stop("Could not find associated CRA.")
      if (length(cras) > 1) cli::cli_alert_warning("Multiple CRAs found ({.val {cras}}), using the first: {.val {cras[1]}}")
      cra_id <- cras[1]
      cli::cli_alert_success("Found associated project {.val {cra_id}} for {.val {accession}}.")
      runinfo_url <- paste0(base_url, "/search/getRunInfoByCra")
      form_data <- list(searchTerm = cra_id, totalDatas = "9999", downLoadCount = "9999")
    }, error = function(e) {
      cli::cli_abort("Failed to find associated CRA for GSA accession {.val {accession}}: {e$message}")
    })
  } else if (stringr::str_detect(accession_upper, "^PRJC[A-Z][0-9]+$|^SAMC[0-9]+$")) {
    runinfo_url <- paste0(base_url, "/search/getRunInfo")
    form_data <- list(searchTerm = paste0('"', accession, '"'), totalDatas = "9999", downLoadCount = "9999")
  } else if (stringr::str_detect(accession_upper, "^CRA[0-9]+$")) {
    cra_id <- accession
    runinfo_url <- paste0(base_url, "/search/getRunInfoByCra")
    form_data <- list(searchTerm = accession, totalDatas = "9999", downLoadCount = "9999")
  } else {
    cli::cli_abort("Invalid GSA accession format: {.val {accession}}.")
  }

  metadata_found_csv <- FALSE
  metadata_df <- NULL
  cli::cli_alert_info("Fetching GSA run metadata for {.val {accession}}...")
  tryCatch({
    resp_runinfo <- httr::POST(runinfo_url,
                               body = form_data, encode = "form",
                               httr::write_disk(temp_file, overwrite = TRUE), httr::progress(), httr::timeout(120)
    )
    httr::stop_for_status(resp_runinfo, task = "fetch GSA runinfo")

    if (fs::file_size(temp_file) > 100) {
      if (stringr::str_detect(accession_upper, "^CRR[0-9]+$|^CRX[0-9]+$")) {
        all_data <- readr::read_delim(temp_file, delim = ",", show_col_types = FALSE, guess_max = 5000)
        run_col <- names(all_data)[stringr::str_detect(names(all_data), "Run")]
        exp_col <- names(all_data)[stringr::str_detect(names(all_data), "Experiment")]
        relevant_data <- all_data %>%
          dplyr::filter(if (length(run_col) > 0) .data[[run_col]] == accession_upper else FALSE |
                          if (length(exp_col) > 0) .data[[exp_col]] == accession_upper else FALSE)
        if (nrow(relevant_data) > 0) {
          metadata_df <- relevant_data
          metadata_found_csv <- TRUE
        } else {
          cli::cli_warn("Could not find specific run/experiment {.val {accession}} within the fetched CRA data.")
        }
      } else {
        metadata_df <- readr::read_delim(temp_file, delim = ",", show_col_types = FALSE, guess_max = 5000)
        metadata_found_csv <- TRUE
      }

      if (metadata_found_csv) {
        if (!check_file_valid(temp_file, accession)) {
          metadata_found_csv <- FALSE
        } else {
          cli::cli_alert_success("Successfully downloaded GSA metadata for {.val {accession}}")
          if (!is.null(outdir)) readr::write_csv(metadata_df, metadata_file_csv)
        }
      }
    } else {
      cli::cli_alert_warning("Downloaded GSA metadata for {.val {accession}} appears empty or invalid.")
    }
  }, error = function(e) {
    cli::cli_alert_danger("Failed to fetch GSA metadata for {.val {accession}}: {e$message}")
  })

  if (!metadata_found_csv) cli::cli_abort("Failed to retrieve valid GSA metadata CSV for {.val {accession}}.")

  return(metadata_df)
}


#' Get Metadata Wrapper
#'
#' Determines the source (GSA or ENA/SRA) and calls the appropriate function.
#'
#' @param accession Accession ID.
#' @param outdir Output directory. If NULL, only returns data.frame.
#' @return Data frame of metadata.
#' @keywords internal
get_metadata <- function(accession, outdir = NULL) {
  accession_upper <- toupper(accession)
  metadata <- NULL

  if (stringr::str_detect(accession_upper, "^(PRJC|SAMC|CRA|CRX|CRR)")) {
    if (!is.null(outdir)) {
      gsa_meta_csv <- fs::path(outdir, paste0(accession, ".metadata.csv"))
      if (fs::file_exists(gsa_meta_csv)) {
        cli::cli_alert_info("GSA metadata file exists: {.path {gsa_meta_csv}}")
        if (!check_file_valid(gsa_meta_csv, accession)) {
          fs::file_delete(gsa_meta_csv)
          cli::cli_alert_warning("Existing GSA metadata file was empty, attempting re-download.")
          metadata <- get_gsa_metadata(accession, outdir)
        } else {
          metadata <- readr::read_csv(gsa_meta_csv, show_col_types = FALSE)
        }
      } else {
        metadata <- get_gsa_metadata(accession, outdir)
      }
    } else {
      metadata <- get_gsa_metadata(accession, outdir)
    }
    attr(metadata, "source") <- "GSA"
  } else {
    if (!is.null(outdir)) {
      ena_meta_tsv <- fs::path(outdir, paste0(accession, ".metadata.tsv"))
      if (fs::file_exists(ena_meta_tsv)) {
        cli::cli_alert_info("ENA/SRA metadata file exists: {.path {ena_meta_tsv}}")
        if (!check_file_valid(ena_meta_tsv, accession)) {
          fs::file_delete(ena_meta_tsv)
          cli::cli_alert_warning("Existing ENA/SRA metadata file was empty, attempting re-download.")
          metadata <- get_ena_sra_metadata(accession, outdir)
        } else {
          metadata <- readr::read_tsv(ena_meta_tsv, show_col_types = FALSE)
        }
      } else {
        metadata <- get_ena_sra_metadata(accession, outdir)
      }
    } else {
      metadata <- get_ena_sra_metadata(accession, outdir)
    }
    attr(metadata, "source") <- "ENA/SRA"
  }

  if (is.null(metadata) || nrow(metadata) == 0) {
    cli::cli_abort("Failed to obtain valid metadata for {.val {accession}}.")
  }

  return(metadata)
}
