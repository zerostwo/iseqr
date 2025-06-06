% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{get_ena_sra_metadata}
\alias{get_ena_sra_metadata}
\title{Fetch Metadata from ENA/SRA}
\usage{
get_ena_sra_metadata(accession, output_dir = ".")
}
\arguments{
\item{accession}{The accession ID.}

\item{output_dir}{Directory to save metadata file.}
}
\value{
Path to the metadata TSV file or NULL if failed.
}
\description{
Fetch Metadata from ENA/SRA
}
\keyword{internal}
