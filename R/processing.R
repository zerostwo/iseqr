# --- Post-processing Functions (processing.R) ---

# Count lines in a file (FASTQ, possibly gzipped)
count_lines <- function(file) {
  if (grepl("\\.gz$", file)) {
    as.integer(system2("zcat", shQuote(file), stdout = TRUE) |> length())
  } else {
    system2("wc", c("-l", shQuote(file)), stdout = TRUE) |>
      strsplit("\\s+") |>
      unlist() |>
      (\(x) x[1])() |>
      as.integer()
  }
}

# Check all files are multiple of 4 and have equal line count
check_fastq_integrity <- function(files) {
  lines <- purrr::map_int(files, count_lines)

  valid <- all(lines %% 4 == 0) && length(unique(lines)) == 1
  list(valid = valid, lines = lines)
}

#' Convert SRA to FASTQ using fasterq-dump
#'
#' @param sra_file Path to the SRA file.
#' @param output_dir Directory to write FASTQ files.
#' @param threads Number of threads to use.
#' @param software_paths Paths to tools (esp. fasterq-dump).
#' @return Character vector of generated FASTQ file paths.
#' @keywords internal
convert_sra_to_fastq <- function(sra_file, output_dir = ".", threads = 8, software_paths) {
  run_id <- fs::path_ext_remove(fs::path_file(sra_file))
  cli::cli_alert_info("Checking existing FASTQ for {.val {run_id}}...")

  fastq_pattern <- paste0(run_id, "*.fastq")
  generated_files <- fs::dir_ls(
    path = output_dir,
    regexp = fastq_pattern
  )

  if (length(generated_files) > 0) {
    cli::cli_alert_info("Found existing FASTQ file(s) for {.val {run_id}}, validating integrity...")

    integrity <- check_fastq_integrity(generated_files)

    if (!integrity$valid) {
      cli::cli_alert_warning("FASTQ files for {.val {run_id}} failed integrity check (line count issue). Re-generating...")
      fs::file_delete(generated_files)
    } else {
      cli::cli_alert_success("FASTQ files passed integrity check (all files: {integrity$lines[1]} lines).")
      cli::cli_ul(generated_files)
      return(generated_files)
    }
  }

  cli::cli_alert_info("Converting SRA {.val {run_id}} to FASTQ using {.val fasterq-dump}...")

  if (is.null(software_paths$"fasterq_dump")) {
    software_paths$"fasterq_dump" <- get_software_path("fasterq_dump", "sra-tools", required = TRUE)
    if (is.null(software_paths$"fasterq_dump")) {
      cli::cli_abort("Cannot find {.val fasterq-dump} for SRA conversion.")
    }
  }

  fs::dir_create(output_dir)

  fastq_dump_args <- c(
    "-p", "-S", "--include-technical",
    "-e", as.character(threads),
    "-O", output_dir,
    sra_file
  )

  on.exit(
    {
      temp_dump_files <- fs::dir_ls(path = output_dir, glob = "fasterq.tmp*", type = "any")
      if (length(temp_dump_files) > 0) {
        cli::cli_inform("Cleaning up temporary fasterq-dump files...")
        fs::file_delete(temp_dump_files)
      }
    },
    add = TRUE
  )

  tryCatch(
    {
      run_command(software_paths$"fasterq_dump", args = fastq_dump_args, spinner = TRUE)
    },
    error = function(e) {
      cli::cli_abort("fasterq-dump failed for {.val {run_id}}.")
    }
  )

  generated_files <- fs::dir_ls(
    path = output_dir,
    regexp = fastq_pattern
  )

  if (length(generated_files) == 0) {
    cli::cli_alert_warning("fasterq-dump ran but no FASTQ files found for {.val {run_id}} in {.path {output_dir}}.")
  } else {
    cli::cli_alert_success("Successfully converted {.val {run_id}} to {length(generated_files)} FASTQ file(s):")
    cli::cli_ul(generated_files)
  }

  return(generated_files)
}

#' Compress FASTQ files using pigz
#'
#' @param fastq_files Character vector of FASTQ file paths.
#' @param threads Number of threads.
#' @param remove_original Remove original FASTQ after compression?
#' @param software_paths Paths to tools (esp. pigz).
#' @return Character vector of compressed file paths (.fastq.gz).
#' @keywords internal
#' Compress FASTQ files using pigz
#'
#' @param fastq_files Character vector of FASTQ file paths (.fastq).
#' @param threads Number of threads.
#' @param software_paths Named list of software paths, especially `pigz`.
#' @return Character vector of successfully compressed file paths (.fastq.gz).
#' @keywords internal
compress_fastq <- function(fastq_files, threads = 8, software_paths) {
  if (length(fastq_files) == 0) {
    return(character())
  }

  if (is.null(software_paths$pigz)) {
    software_paths$pigz <- get_software_path("pigz", "pigz", required = TRUE)
    if (is.null(software_paths$pigz)) {
      cli::cli_abort("Cannot find {.val pigz} for compression.")
    }
  }

  # Determine which files still need compression
  compressed_files <- paste0(fastq_files, ".gz")
  needs_compression <- !fs::file_exists(compressed_files)
  files_to_compress <- fastq_files[needs_compression]
  compressed_to_check <- compressed_files[needs_compression]

  if (length(files_to_compress) == 0) {
    cli::cli_alert_info("All FASTQ files are already compressed. Skipping.")
    return(compressed_files)
  }

  cli::cli_alert_info("Compressing {length(files_to_compress)} FASTQ file(s) using {.val pigz}...")

  success_flags <- rep(TRUE, length(files_to_compress))

  tryCatch(
    {
      run_command(software_paths$pigz, args = c("-p", as.character(threads), files_to_compress), spinner = TRUE)

      # Confirm results
      if (!all(fs::file_exists(compressed_to_check))) {
        cli::cli_warn("Some compressed files were not created. Possible pigz failure.")
        success_flags <- fs::file_exists(compressed_to_check)
      } else {
        cli::cli_alert_success("Successfully compressed FASTQ files.")
      }
    },
    error = function(e) {
      cli::cli_alert_danger("pigz compression failed.")
      success_flags <<- rep(FALSE, length(files_to_compress))
    }
  )

  # Return only successful compressed files (both newly and previously existing)
  final_success <- fs::file_exists(compressed_files)
  return(compressed_files[final_success])
}


#' Merge FASTQ Files
#'
#' Merges FASTQ files based on experiment, sample, or study level.
#' Handles both SRA/ENA (.fastq/.fastq.gz) and GSA (.fastq.gz assumed) files.
#'
#' @param accession The main accession ID used for the download batch.
#' @param metadata_file Path to the metadata file (TSV or CSV).
#' @param source "ENA/SRA" or "GSA".
#' @param merge_level "ex" (experiment), "sa" (sample), or "st" (study).
#' @param output_dir Directory containing the FASTQ files.
#' @param is_gzipped Are the input FASTQ files gzipped?
#' @param remove_originals Remove original run-level FASTQ files after merging?
#' @keywords internal
merge_fastq_files <- function(accession, metadata_file, source, merge_level, output_dir = ".", is_gzipped = TRUE, remove_originals = TRUE) {
  cli::cli_alert_info("Merging FASTQ files for {.val {accession}} at the {.val {merge_level}} level...")

  # Read metadata
  meta <- tryCatch(
    {
      if (source == "GSA") {
        readr::read_csv(metadata_file, show_col_types = FALSE, guess_max = 10000)
      } else { # ENA/SRA TSV
        readr::read_tsv(metadata_file, show_col_types = FALSE, guess_max = 10000)
      }
    },
    error = function(e) {
      cli::cli_abort("Failed to read metadata file {.path {metadata_file}} for merging: {e$message}")
    }
  )

  # Identify grouping column and ID prefix based on source and merge level
  grouping_col <- NULL
  id_prefix <- NULL # Expected prefix for the group ID (e.g., SRX, SAMEA)
  run_col <- NULL # Column containing the run ID (e.g., SRR, CRR)

  if (source == "ENA/SRA") {
    run_col <- "run_accession"
    grouping_col <- switch(merge_level,
      "ex" = "experiment_accession", # SRX, ERX, DRX
      "sa" = "sample_accession", # SAMEA, SAMN, SAMD, ERS, DRS
      "st" = "study_accession" # PRJNA, PRJEB, PRJD, SRP, ERP, DRP
      # Add secondary accessions if needed? e.g. secondary_sample_accession
    )
    # Prefixes are more complex for ENA/SRA due to multiple types per level
    # We might just group by the value in the column directly.
  } else { # GSA
    run_col <- "Run" # Check actual GSA CSV header
    grouping_col <- switch(merge_level,
      "ex" = "Experiment", # CRX
      "sa" = "Sample", # SAMC
      "st" = "BioProject" # PRJCA, CRA
    )
    id_prefix <- switch(merge_level,
      "ex" = "CRX",
      "sa" = "SAMC",
      "st" = "PRJC|CRA"
    )
  }

  if (is.null(grouping_col) || !grouping_col %in% names(meta)) {
    cli::cli_abort("Cannot find grouping column {.val {grouping_col}} in metadata for merge level {.val {merge_level}}.")
  }
  if (is.null(run_col) || !run_col %in% names(meta)) {
    cli::cli_abort("Cannot find run column {.val {run_col}} in metadata.")
  }


  # Group runs by the chosen level
  grouped_runs <- meta %>%
    select(group_id = all_of(grouping_col), run_id = all_of(run_col)) %>%
    distinct() %>%
    group_by(group_id) %>%
    summarise(runs = list(run_id), .groups = "drop")

  file_suffix <- if (is_gzipped) ".fastq.gz" else ".fastq"

  # Process each group
  for (i in 1:nrow(grouped_runs)) {
    group_id <- grouped_runs$group_id[i]
    runs_in_group <- grouped_runs$runs[[i]]

    cli::cli_alert("Processing group: {.val {group_id}} with runs: {.val {runs_in_group}}")

    # Find FASTQ files for these runs
    # Expected patterns: {run_id}{file_suffix}, {run_id}_1{file_suffix}, {run_id}_2{file_suffix}, etc.
    # GSA files might have different naming conventions based on 'ftpFile' column.
    # Let's assume standard SRA-Toolkit naming convention first.

    files_to_merge_r1 <- character()
    files_to_merge_r2 <- character()
    files_to_merge_se <- character() # Single-end or potentially technical reads
    all_original_files <- character()

    is_paired <- FALSE # Track if we find paired files

    for (run_id in runs_in_group) {
      # Find files for this run in the output directory
      # run_files <- fs::dir_ls(output_dir, glob = paste0(run_id, "*", file_suffix), recurse = FALSE)
      run_files <- fs::dir_ls(output_dir, regexp = paste0(run_id, file_suffix), recurse = FALSE)
      all_original_files <- c(all_original_files, run_files)

      # r1_file <- run_files[stringr::str_detect(fs::path_file(run_files), "_1" %R% file_suffix %R% "$")]
      # r2_file <- run_files[stringr::str_detect(fs::path_file(run_files), "_2" %R% file_suffix %R% "$")]
      escaped_suffix <- gsub("\\.", "\\\\.", file_suffix) # Replaces . with \\.
      pattern_r1 <- paste0("_1", escaped_suffix, "$")
      pattern_r2 <- paste0("_2", escaped_suffix, "$")
      r1_file <- run_files[stringr::str_detect(fs::path_file(run_files), pattern_r1)]
      r2_file <- run_files[stringr::str_detect(fs::path_file(run_files), pattern_r2)]

      # SE file: ends with suffix but not _1 or _2 suffix
      pattern_se_end <- paste0(escaped_suffix, "$")
      pattern_se_not_12 <- paste0("_[12]", escaped_suffix, "$") # Pattern to exclude
      # se_file <- run_files[stringr::str_detect(fs::path_file(run_files), file_suffix %R% "$") &
      #                        !stringr::str_detect(fs::path_file(run_files), "_[12]" %R% file_suffix %R% "$")]
      # se_file <- run_files[stringr::str_detect(fs::path_file(run_files), file_suffix %R% "$") &
      #                        !stringr::str_detect(fs::path_file(run_files), "_[12]" %R% file_suffix %R% "$")]
      se_file <- run_files[stringr::str_detect(fs::path_file(run_files), pattern_se_end) &
        !stringr::str_detect(fs::path_file(run_files), pattern_se_not_12)]
      # Other files (e.g., _3, _I1, etc.) - treat as single-end for now? Or handle separately?
      other_files <- run_files[!run_files %in% c(r1_file, r2_file, se_file)]


      if (length(r1_file) == 1 && length(r2_file) == 1) {
        files_to_merge_r1 <- c(files_to_merge_r1, r1_file)
        files_to_merge_r2 <- c(files_to_merge_r2, r2_file)
        is_paired <- TRUE
        # Add any 'other' files to the single-end list? Or ignore? Let's add them.
        files_to_merge_se <- c(files_to_merge_se, other_files)
      } else if (length(se_file) >= 1) { # Can be multiple SE files (e.g. technical + biological)
        files_to_merge_se <- c(files_to_merge_se, se_file, r1_file, r2_file, other_files) # Add everything found
      } else if (length(r1_file) == 1 && length(r2_file) == 0) { # R1 only
        files_to_merge_se <- c(files_to_merge_se, r1_file, other_files)
      } else {
        cli::cli_alert_warning("Could not determine file layout (paired/single) for run {.val {run_id}} in group {.val {group_id}}. Found: {.path {run_files}}. Adding all to single-end merge.")
        files_to_merge_se <- c(files_to_merge_se, run_files)
      }
    }

    # Perform the merge using `cat` via system2 (safer for large files)
    perform_merge <- function(input_files, output_file) {
      if (length(input_files) == 0) {
        return(FALSE)
      }
      if (length(input_files) == 1) {
        cli::cli_inform("Only one file for group {.val {group_id}} ({output_file}), renaming {.path {input_files[1]}} -> {.path {output_file}}")
        fs::file_move(input_files[1], output_file)
        return(TRUE)
      }

      cli::cli_inform("Merging {length(input_files)} files into {.path {output_file}}...")
      cli::cli_ul(input_files)
      id <- cli::cli_status("Merging files...")

      # Use processx::run with output redirection
      # Safer than system2("cat ... > ...") as it handles arguments better
      # Create a temporary file list if too many files for command line limits? Usually not an issue.
      # processx can handle redirection.

      # Using file connections in R (less safe for huge files):
      # tryCatch({
      #     out_con <- file(output_file, "wb") # Binary mode if gzipped
      #     on.exit(close(out_con))
      #     for (infile in input_files) {
      #         in_con <- file(infile, "rb")
      #         # Read and write in chunks - requires more code
      #         # Simple version (reads all into memory per file!)
      #         # content <- readBin(in_con, "raw", n = fs::file_size(infile))
      #         # writeBin(content, out_con)
      #         # close(in_con)
      #         # SAFER: Use cat command
      #     }
      # }, error = ...)

      # Use system `cat` - most robust for large files
      # Need to handle potential spaces/special chars in filenames
      # Quote filenames for safety
      # quoted_files <- shQuote(input_files)

      # Construct command. Redirection needs careful handling.
      # Let processx handle redirection.
      # This writes stdout to the output file.
      res <- tryCatch(
        {
          processx::run("cat", args = input_files, stdout = output_file, error_on_status = FALSE, spinner = FALSE, echo_cmd = TRUE)
        },
        error = function(e) {
          cli::cli_status_clear(id)
          cli::cli_alert_danger("Failed to run cat command: {e$message}")
          return(NULL)
        }
      )

      cli::cli_status_clear(id)

      if (is.null(res) || res$status != 0) {
        cli::cli_alert_danger("Merging failed for {.path {output_file}} (exit code: {res$status}).")
        if (!is.null(res$stderr) && nzchar(res$stderr)) cli::cli_verbatim(res$stderr)
        if (fs::file_exists(output_file)) fs::file_delete(output_file) # Clean up partial merge
        return(FALSE)
      } else {
        cli::cli_alert_success("Successfully merged files into {.path {output_file}}")
        return(TRUE)
      }
    }

    merge_successful <- TRUE
    merged_files_list <- character()

    # Merge R1 files if any
    if (length(files_to_merge_r1) > 0) {
      merged_r1_file <- fs::path(output_dir, paste0(group_id, "_1", file_suffix))
      # Check if target exists and skip?
      if (fs::file_exists(merged_r1_file)) {
        cli::cli_alert_warning("Merged file already exists, skipping: {.path {merged_r1_file}}")
      } else {
        if (perform_merge(files_to_merge_r1, merged_r1_file)) {
          merged_files_list <- c(merged_files_list, merged_r1_file)
        } else {
          merge_successful <- FALSE
        }
      }
    }

    # Merge R2 files if any (only if R1 merge was attempted/successful)
    if (length(files_to_merge_r2) > 0 && merge_successful) {
      merged_r2_file <- fs::path(output_dir, paste0(group_id, "_2", file_suffix))
      if (fs::file_exists(merged_r2_file)) {
        cli::cli_alert_warning("Merged file already exists, skipping: {.path {merged_r2_file}}")
      } else {
        if (perform_merge(files_to_merge_r2, merged_r2_file)) {
          merged_files_list <- c(merged_files_list, merged_r2_file)
        } else {
          merge_successful <- FALSE
        }
      }
    }

    # Merge SE files if any (only if paired merges weren't done or failed, or if only SE files exist)
    if (length(files_to_merge_se) > 0 && merge_successful) {
      # If paired data was successfully merged, we usually don't merge SE files unless they are technical reads
      # The current logic adds _all_ non-R1/R2 files here. Decide if that's correct.
      # Let's merge them into a single SE file for the group if they exist.
      merged_se_file <- fs::path(output_dir, paste0(group_id, file_suffix)) # No _1/_2 suffix

      # Avoid overwriting if only R1 existed and was renamed to this name
      is_renamed_r1 <- length(files_to_merge_r1) == 1 && length(files_to_merge_r2) == 0 && length(files_to_merge_se) == 0

      if (fs::file_exists(merged_se_file) && !is_renamed_r1) {
        cli::cli_alert_warning("Merged file already exists, skipping: {.path {merged_se_file}}")
      } else if (!is_renamed_r1) {
        if (perform_merge(files_to_merge_se, merged_se_file)) {
          merged_files_list <- c(merged_files_list, merged_se_file)
        } else {
          merge_successful <- FALSE
        }
      }
    }


    # Clean up original run files if merge was successful
    if (merge_successful && remove_originals && length(merged_files_list) > 0) {
      # Only remove files that were actually part of a successful merge group
      files_to_remove <- unique(all_original_files)
      # Exclude files that are identical to the merged output (in case of single-file rename)
      files_to_remove <- files_to_remove[!files_to_remove %in% merged_files_list]

      if (length(files_to_remove) > 0) {
        cli::cli_inform("Removing {length(files_to_remove)} original run file(s) for group {.val {group_id}}...")
        fs::file_delete(files_to_remove)
      }
    } else if (!merge_successful) {
      cli::cli_alert_danger("Merging failed for group {.val {group_id}}, original files kept.")
      # Clean up any partially created merged files
      fs::file_delete(merged_files_list[fs::file_exists(merged_files_list)])
    }
  } # End loop over groups
}
