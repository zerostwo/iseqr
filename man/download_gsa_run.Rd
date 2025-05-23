% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{download_gsa_run}
\alias{download_gsa_run}
\title{Download a Single GSA Run}
\usage{
download_gsa_run(
  run_id,
  run_info,
  use_aspera = FALSE,
  parallel_n = 0,
  output_dir = ".",
  software_paths
)
}
\arguments{
\item{run_id}{CRR identifier.}

\item{run_info}{A single row tibble/data.frame from the metadata file for this run.}

\item{use_aspera}{Should Aspera be used if available?}

\item{parallel_n}{Parallel connections for axel.}

\item{output_dir}{Destination directory.}

\item{software_paths}{Paths to tools.}
}
\value{
A list containing paths to downloaded files. Key: \code{files}. Type is implicit (various).
}
\description{
Handles logic for downloading GSA files (FASTQ, BAM, etc.) using FTP/HTTPS/Aspera/Huawei.
}
\keyword{internal}
