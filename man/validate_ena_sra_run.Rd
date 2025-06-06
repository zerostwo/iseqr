% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{validate_ena_sra_run}
\alias{validate_ena_sra_run}
\title{Validate Downloaded ENA/SRA Run}
\usage{
validate_ena_sra_run(run_id, downloaded_info, run_info, software_paths)
}
\arguments{
\item{run_id}{SRR/ERR/DRR identifier.}

\item{downloaded_info}{List returned by \code{download_ena_sra_run} (contains \code{files}, \code{type}).}

\item{run_info}{Metadata row for the run.}

\item{software_paths}{Paths to tools (esp. vdb-validate, md5sum).}
}
\value{
Logical TRUE if validation passes, FALSE otherwise.
}
\description{
Uses vdb-validate for SRA files or MD5 checksums for FASTQ files if available.
}
\keyword{internal}
