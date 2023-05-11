\name{filter_variants}
\alias{filter_variants}
\title{
Filter variants from CSV files based on variant type and various
filtering parameters
}
\usage{
filter_variants(subdir_list, variant_type, maxpopfreq = NULL,
gatkcounts = NULL, MPCscore = NULL, filename_filter)
}
\arguments{
\item{subdir_list}{a character vector specifying the subdirectories
containing the CSV files to be filtered}
\item{variant_type}{a character string specifying the type of variant to
filter for. Allowed values are "PTV" (for protein-truncating
variants) and "missense" (for missense variants)}
\item{maxpopfreq}{an optional numeric value specifying the maximum allowed
population frequency of a variant. If specified, variants with
a maximum population frequency greater than this value will be
filtered out.}
\item{gatkcounts}{an optional integer value specifying the maximum allowed
GATK count of a variant. If specified, variants with a GATK count
greater than this value will be filtered out.}
\item{MPCscore}{an optional numeric value specifying the maximum allowed MPC 
score of a variant. If specified, variants with an MPC score greater than
this value will be filtered.}
\item{filename_filter}{a character string to be included in the filtered file
names to distinguish them from the original files}
}
\description{
This function reads in CSV files containing variant data, and filters the
variants based on the specified type (either "PTV" or "\emph{missense}") and 
various optional filtering parameters. The filtered data is saved as a new CSV 
file.
}
\details{
For each specified subdirectory, the function looks for CSV files with names
starting with "HI_" and filters the variants based on the input parameters.
The filtered data is saved as new CSV files with the same names as the original
files, with the addition of the string "filtered_" and the input parameters to
the file name.
}
\examples{

Filter for missense variants with a maximum population frequency of 0.01
and a GATK count less than 5 in the specified subdirectories
subdir_list <- c("path/to/subdir1", "path/to/subdir2")
filter_variants(subdir_list, "missense", maxpopfreq = 0.01, gatkcounts = 5,
filename_filter = "example")
}

\author{Nermina Logo Lendo}
\keyword{filter, variant}