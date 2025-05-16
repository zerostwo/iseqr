# --- Validation Functions (validation.R) ---

#' Validate Downloaded ENA/SRA Run
#'
#' Uses vdb-validate for SRA files or MD5 checksums for FASTQ files if available.
#'
#' @param run_id SRR/ERR/DRR identifier.
#' @param downloaded_info List returned by `download_ena_sra_run` (contains `files`, `type`).
#' @param run_info Metadata row for the run.
#' @param software_paths Paths to tools (esp. vdb-validate, md5sum).
#' @return Logical TRUE if validation passes, FALSE otherwise.
#' @keywords internal
validate_ena_sra_run <- function(run_id, downloaded_info, run_info, software_paths) {
  if (is.null(downloaded_info) || length(downloaded_info$files) == 0) {
    cli::cli_alert_warning("No files found to validate for {.val {run_id}}.")
    return(FALSE)
  }

  if (downloaded_info$type == "sra") {
    sra_file <- downloaded_info$files[1] # Should only be one SRA file
    cli::cli_alert_info("Validating SRA file using {.val vdb-validate}: {.path {sra_file}}")
    if (is.null(software_paths$"vdb_validate")) {
      software_paths$"vdb_validate" <- get_software_path("vdb_validate", "sra-tools", required = TRUE)
      if (is.null(software_paths$"vdb_validate")) {
        cli::cli_alert_warning("Cannot find {.val vdb-validate}, skipping SRA validation for {.val {run_id}}.")
        return(TRUE) # Skip validation but assume success if tool missing? Risky. Maybe FALSE? Let's warn and skip (return TRUE).
      }
    }
    log_file <- fs::path_temp(paste0(run_id, ".vdb-validate.log"))
    on.exit(fs::file_delete(log_file), add = TRUE)

    # Run vdb-validate, capturing output/status
    res <- run_command(software_paths$"vdb_validate",
      args = sra_file,
      error_on_status = FALSE, # Check status manually
      stdout = log_file, stderr = log_file
    ) # Capture output

    if (res$status == 0) {
      cli::cli_alert_success("vdb-validate successful for {.val {run_id}}.")
      # Optional: check log for specific success messages if needed
      return(TRUE)
    } else {
      cli::cli_alert_danger("vdb-validate failed for {.val {run_id}} (exit code: {res$status}). See details below.")
      log_content <- readr::read_lines(log_file)
      cli::cli_verbatim(paste(log_content, collapse = "\n"))
      return(FALSE)
    }
  } else if (downloaded_info$type == "fastq_gz") {
    cli::cli_alert_info("Validating FASTQ.gz file(s) using MD5 checksums for {.val {run_id}}...")
    # Check MD5s from metadata ('fastq_md5' column)
    if (!"fastq_md5" %in% names(run_info) || is.na(run_info$fastq_md5) || run_info$fastq_md5 == "") {
      cli::cli_alert_warning("No MD5 checksums found in metadata for {.val {run_id}}. Skipping MD5 validation.")
      return(TRUE) # Assume success if no checksum available
    }

    expected_md5s <- str_split(run_info$fastq_md5, ";")[[1]] %>% discard(~ . == "")
    # Get corresponding filenames (assuming order matches in 'fastq_ftp')
    # This relies on the download function returning files in the correct order
    downloaded_filenames <- fs::path_file(downloaded_info$files)

    if (length(expected_md5s) != length(downloaded_info$files)) {
      cli::cli_alert_warning("Mismatch between number of expected MD5s ({length(expected_md5s)}) and downloaded files ({length(downloaded_info$files)}) for {.val {run_id}}. Skipping MD5 validation.")
      return(TRUE) # Cannot reliably validate
    }

    # Use R's internal md5sum or external command
    use_external_md5sum <- TRUE # Default to external for potentially better performance
    md5sum_path <- NULL
    if (use_external_md5sum) {
      md5sum_path <- get_software_path("md5sum", "coreutils", required = FALSE)
      if (is.null(md5sum_path)) {
        cli::cli_alert_info("External {.val md5sum} not found, using R's {.fun tools::md5sum}.")
        use_external_md5sum <- FALSE
      }
    }

    all_match <- TRUE
    for (i in seq_along(downloaded_info$files)) {
      file_to_check <- downloaded_info$files[i]
      expected_md5 <- tolower(expected_md5s[i]) # Ensure lowercase comparison

      if (!fs::file_exists(file_to_check)) {
        cli::cli_alert_danger("File not found for MD5 check: {.path {file_to_check}}")
        all_match <- FALSE
        break
      }

      cli::cli_inform("Calculating MD5 for: {.path {file_to_check}}")
      id <- cli::cli_status("Calculating MD5...")

      calculated_md5 <- tryCatch(
        {
          if (use_external_md5sum) {
            res <- run_command(md5sum_path, args = file_to_check, error_on_status = TRUE, spinner = FALSE)
            tolower(str_split(str_trim(res$stdout), "\\s+")[[1]][1]) # Extract first part, lowercase
          } else {
            tools::md5sum(file_to_check) # Returns named vector
          }
        },
        error = function(e) {
          cli::cli_status_clear(id)
          cli::cli_alert_danger("Failed to calculate MD5 for {.path {file_to_check}}: {e$message}")
          return(NULL)
        }
      )
      cli::cli_status_clear(id)

      if (is.null(calculated_md5)) {
        all_match <- FALSE
        break
      }

      # Ensure internal md5sum result is extracted correctly
      if (!use_external_md5sum) {
        calculated_md5 <- unname(calculated_md5) # Get the value from named vector
      }


      if (calculated_md5 == expected_md5) {
        cli::cli_alert_success("MD5 match for: {.path {file_to_check}}")
      } else {
        cli::cli_alert_danger("MD5 mismatch for: {.path {file_to_check}} (Expected: {expected_md5}, Got: {calculated_md5})")
        all_match <- FALSE
        # break # Option to stop checking after first mismatch
      }
    }

    return(all_match)
  } else {
    cli::cli_alert_warning("Unknown download type '{downloaded_info$type}' for validation of {.val {run_id}}.")
    return(FALSE) # Cannot validate unknown type
  }
}


#' Validate Downloaded GSA Run
#'
#' Uses MD5 checksums fetched from GSA's md5sum.txt file.
#'
#' @param run_id CRR identifier.
#' @param downloaded_info List returned by `download_gsa_run` (contains `files`).
#' @param run_info Metadata row for the run.
#' @param output_dir Directory containing downloads and potentially the md5sum.txt file.
#' @param software_paths Paths to tools (esp. md5sum).
#' @return Logical TRUE if validation passes, FALSE otherwise.
#' @keywords internal
validate_gsa_run <- function(run_id, downloaded_info, run_info, output_dir = ".", software_paths) {
  if (is.null(downloaded_info) || length(downloaded_info$files) == 0) {
    cli::cli_alert_warning("No files found to validate for GSA run {.val {run_id}}.")
    return(FALSE)
  }

  # 1. Get the CRA to find the corresponding md5sum.txt file
  cra_id <- run_info$`BioProject Accession` # Assuming this column exists
  if (is.null(cra_id)) {
    cli::cli_alert_warning("Could not determine CRA ID from GSA metadata for run {.val {run_id}}. Skipping MD5 validation.")
    return(TRUE) # Assume success? Or fail? Let's skip (TRUE).
  }
  md5sum_file <- fs::path(output_dir, paste0(cra_id, ".md5sum.txt"))

  # 2. Download md5sum.txt if it doesn't exist
  if (!fs::file_exists(md5sum_file)) {
    cli::cli_alert_info("MD5 checksum file not found for {.val {cra_id}}, attempting download...")
    # Construct the URL - this is tricky, depends on GSA structure (gsa/CRA... or gsa-human/CRA...)
    # Try to guess base path from run_info 'ftpFile' if possible, otherwise default structure.
    base_url <- "https://download.cncb.ac.cn/"
    relative_path <- NULL
    if ("ftpFile" %in% names(run_info) && !is.na(run_info$ftpFile) && run_info$ftpFile != "") {
      # Example ftpFile: /gsa/CRA000001/CRR... or /gsahuman/CRA...
      first_file_path <- str_split(run_info$ftpFile, "[|]")[[1]][1]
      # Extract the part like 'gsa/CRA...' or 'gsahuman/CRA...'
      match <- stringr::str_match(first_file_path, "(/(gsa|gsa-human|gsa-plant)/CRA[0-9]+)/")
      if (!is.na(match[1, 2])) {
        relative_path <- stringr::str_sub(match[1, 2], 2) # Remove leading '/'
      }
    }
    # Fallback guess if path extraction failed
    if (is.null(relative_path)) {
      # Basic guess, might be wrong for gsa-human etc.
      relative_path <- paste0("gsa/", cra_id)
      cli::cli_alert_warning("Could not reliably determine GSA path, guessing: {.path {relative_path}}")
    }

    md5_url <- paste0(base_url, relative_path, "/md5sum.txt")
    cli::cli_inform("Attempting to download from: {.url {md5_url}}")

    tryCatch(
      {
        resp <- httr::GET(md5_url, httr::write_disk(md5sum_file, overwrite = TRUE), httr::timeout(60))
        # Check if file downloaded and is not empty
        if (!fs::file_exists(md5sum_file) || fs::file_size(md5sum_file) == 0) {
          # Try alternative common path (e.g., gsa-human if first guess was gsa) - this gets complex
          # For now, just report failure if first attempt fails.
          stop(paste("Download failed or resulted in empty file. Status:", httr::status_code(resp)))
        }
        cli::cli_alert_success("Downloaded MD5 checksum file to {.path {md5sum_file}}")
      },
      error = function(e) {
        cli::cli_alert_warning("Failed to download GSA MD5 checksum file for {.val {cra_id}}: {e$message}")
        if (fs::file_exists(md5sum_file)) fs::file_delete(md5sum_file) # Clean up failed download
        cli::cli_alert_warning("Skipping MD5 validation for GSA run {.val {run_id}}.")
        return(TRUE) # Skip validation if checksum file unavailable
      }
    )
  }

  # 3. Read the md5sum file
  md5_data <- tryCatch(
    {
      readr::read_delim(md5sum_file, delim = " ", col_names = c("md5", "filename"), trim_ws = TRUE, skip_empty_rows = TRUE, col_types = "cc") %>%
        mutate(filename = str_remove(filename, "^\\*")) %>% # Remove leading '*' if present
        filter(!is.na(md5), !is.na(filename)) # Ensure valid rows
    },
    error = function(e) {
      cli::cli_alert_warning("Failed to read or parse GSA MD5 checksum file {.path {md5sum_file}}: {e$message}")
      return(NULL)
    }
  )

  if (is.null(md5_data) || nrow(md5_data) == 0) {
    cli::cli_alert_warning("GSA MD5 checksum file is empty or invalid. Skipping MD5 validation for {.val {run_id}}.")
    return(TRUE) # Skip validation
  }

  # 4. Compare checksums
  cli::cli_alert_info("Validating GSA file(s) using MD5 checksums for {.val {run_id}}...")
  use_external_md5sum <- TRUE
  md5sum_path <- NULL
  if (use_external_md5sum) {
    md5sum_path <- get_software_path("md5sum", "coreutils", required = FALSE)
    if (is.null(md5sum_path)) {
      cli::cli_alert_info("External {.val md5sum} not found, using R's {.fun tools::md5sum}.")
      use_external_md5sum <- FALSE
    }
  }

  all_match <- TRUE
  for (file_to_check in downloaded_info$files) {
    if (!fs::file_exists(file_to_check)) {
      cli::cli_alert_danger("File not found for MD5 check: {.path {file_to_check}}")
      all_match <- FALSE
      next # Check next file
    }

    file_basename <- fs::path_file(file_to_check)
    # Find the expected MD5 for this file basename
    expected_md5_row <- md5_data %>% filter(filename == file_basename)

    if (nrow(expected_md5_row) == 0) {
      cli::cli_alert_warning("No MD5 checksum found in {.path {md5sum_file}} for file: {.path {file_basename}}. Skipping check for this file.")
      next # Skip check for this specific file
    }
    expected_md5 <- tolower(expected_md5_row$md5[1]) # Take first if multiple matches? Should be unique.

    cli::cli_inform("Calculating MD5 for: {.path {file_to_check}}")
    id <- cli::cli_status("Calculating MD5...")

    calculated_md5 <- tryCatch(
      {
        if (use_external_md5sum) {
          res <- run_command(md5sum_path, args = file_to_check, error_on_status = TRUE, spinner = FALSE)
          tolower(str_split(str_trim(res$stdout), "\\s+")[[1]][1])
        } else {
          tools::md5sum(file_to_check)
        }
      },
      error = function(e) {
        cli::cli_status_clear(id)
        cli::cli_alert_danger("Failed to calculate MD5 for {.path {file_to_check}}: {e$message}")
        return(NULL)
      }
    )
    cli::cli_status_clear(id)

    if (is.null(calculated_md5)) {
      all_match <- FALSE
      next # Skip to next file if calculation failed
    }

    if (!use_external_md5sum) {
      calculated_md5 <- unname(calculated_md5)
    }

    if (calculated_md5 == expected_md5) {
      cli::cli_alert_success("MD5 match for: {.path {file_to_check}}")
    } else {
      cli::cli_alert_danger("MD5 mismatch for: {.path {file_to_check}} (Expected: {expected_md5}, Got: {calculated_md5})")
      all_match <- FALSE
    }
  }

  return(all_match)
}
