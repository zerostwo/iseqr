#' Get Path to External Software
#'
#' Checks PATH environment variable and R options for software executable.
#'
#' @param software_name The name of the executable (e.g., "wget").
#' @param package_name Optional: Name of the package/tool for error messages (e.g., "wget", "sra-tools").
#' @param required Logical: If TRUE, throws an error if not found. If FALSE, returns NULL.
#' @return The path to the executable or NULL/Error.
#' @keywords internal
get_software_path <- function(software_name, package_name = software_name, required = TRUE) {
  path <- Sys.which(software_name)

  if (path == "") {
    option_name <- paste0("iseq.", software_name)
    path_option <- getOption(option_name, "")
    if (path_option != "" && fs::file_exists(path_option)) {
      path <- path_option
      cli::cli_inform("Using {.val {software_name}} from R option {.val {option_name}}: {.path {path}}")
    }
  }

  if (path == "") {
    if (required) {
      conda_install_cmd <- {
        if (grepl("[><=]", package_name)) {
          paste("conda install -c bioconda", shQuote(package_name), "-y")
        } else {
          paste("conda install -c conda-forge -c bioconda", package_name, "-y")
        }
      }

      cli::cli_abort(c(
        "x" = "{.val {software_name}} not found in PATH or {.val {option_name}} option.",
        "i" = "Please install {.val {package_name}} and add it to your PATH.",
        "i" = "Example installation using conda: {.code {conda_install_cmd}}",
        "i" = "Alternatively, set the R option: {.code options({option_name} = \"/path/to/{software_name}\")}"
      ))
    } else {
      return(NULL)
    }
  }
  # Return normalized path, especially important on Windows
  # return(fs::path_real(path))
  return(path)
}

#' Run External Command
#'
#' Wrapper around processx::run for consistency and error handling.
#'
#' @param command Path to the command.
#' @param args Character vector of arguments.
#' @param error_on_status Should an error be thrown on non-zero exit status?
#' @param spinner Show a spinner?
#' @param ... Additional arguments passed to processx::run.
#' @return Result from processx::run.
#' @keywords internal
run_command <- function(command, args = character(), error_on_status = TRUE, spinner = FALSE, ...) {
  command_name <- fs::path_file(command)
  cli::cli_alert_info("Running: {.file {command_name}} {paste(args, collapse = ' ')}")

  # processx needs command and args separately
  # Ensure command path is quoted if it contains spaces (less common on Linux/Mac)
  # processx handles this reasonably well, but good practice

  # Don't show spinner if stdout/stderr are piped or not interactive
  show_spinner <- spinner && interactive() && is.null(list(...)$stdout) && is.null(list(...)$stderr)

  id <- NULL
  if(show_spinner) id <- cli::cli_status("Running {.file {command_name}}...")

  res <- tryCatch({
    processx::run(
      command = command,
      args = args,
      error_on_status = FALSE, # Check status manually for better error message
      echo_cmd = FALSE,       # We print it ourselves
      echo = FALSE,           # Don't echo output by default
      spinner = FALSE,        # Use cli spinner instead
      ...
    )
  }, error = function(e) {
    if(!is.null(id)) cli::cli_status_clear(id)
    cli::cli_abort("Failed to execute {.file {command_name}}: {e$message}")
  })

  if(!is.null(id)) cli::cli_status_clear(id)

  if (error_on_status && res$status != 0) {
    cli::cli_abort(c(
      "x" = "{.file {command_name}} failed (exit code: {res$status}).",
      "i" = "Command: {.file {command}} {paste(args, collapse = ' ')}",
      if (nzchar(res$stdout)) c(">" = "Standard Output:", cli::cli_verbatim(res$stdout)),
      if (nzchar(res$stderr)) c("!" = "Standard Error:", cli::cli_verbatim(res$stderr))
    ))
  }
  return(res)
}


#' Check File Existence and Size
#'
#' @param file_path Path to the file.
#' @param accession Accession ID for context in error messages.
#' @return TRUE if file exists and is not empty, FALSE otherwise (with warning).
#' @keywords internal
check_file_valid <- function(file_path, accession) {
  if (!fs::file_exists(file_path)) {
    cli::cli_warn("File not found: {.path {file_path}} for accession {.val {accession}}")
    return(FALSE)
  }
  if (fs::file_size(file_path) == 0) {
    cli::cli_warn("File is empty: {.path {file_path}} for accession {.val {accession}}")
    # Consider removing the empty file? fs::file_delete(file_path)
    return(FALSE)
  }
  return(TRUE)
}


#' Log Download Status
#'
#' Appends run ID to success or fail log file.
#'
#' @param run_id The Run accession (e.g., SRR12345).
#' @param status "success" or "fail".
#' @param log_dir Directory where logs are stored.
#' @keywords internal
log_status <- function(run_id, status = "success", log_dir = ".") {
  log_file <- fs::path(log_dir, paste0(status, ".log"))
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste(timestamp, run_id, sep = "\t")
  # Use cat for appending
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

#' Check if Run is Already Successfully Downloaded
#'
#' @param run_id The Run accession.
#' @param log_dir Directory where logs are stored.
#' @return Logical TRUE if found in success.log.
#' @keywords internal
is_already_successful <- function(run_id, log_dir = ".") {
  log_file <- fs::path(log_dir, "success.log")
  if (!fs::file_exists(log_file)) {
    return(FALSE)
  }
  # Read lines and check efficiently
  # Avoid reading huge logs into memory if possible, but for typical use, readLines is fine.
  log_content <- readr::read_lines(log_file, lazy = TRUE) # lazy = TRUE might help slightly
  # Search for the run_id, assuming tab separation and run_id is the second field
  # Using fixed = TRUE for speed if run_id patterns are simple
  any(stringr::str_detect(log_content, paste0("\t", run_id), negate = FALSE))
}
