\name{intersect_func}
\alias{intersect_func}
\title{
  Perform intersect for a list of subdirectories
}
\description{
  This function performs an intersect for each file with the pattern
  "^filtered_de_novo_PTV_" in each subdirectory specified in the input
  list. It reads in the file, writes it as a bed file, performs the
  intersect using bedtools, and writes the output to file.
}
\usage{
  intersect_func(subdirlist)
}
\arguments{
  \item{subdirlist}{
    A character vector of subdirectory paths.
  }
}
\details{
  This function performs an intersect for each file with the pattern
  "^filtered_de_novo_PTV_" in each subdirectory specified in the input
  list. It reads in the file, writes it as a bed file, performs the
  intersect using bedtools, and writes the output to file. The bedtools
  intersect command is used for the intersection, with the input bed file
  as the first argument (-a) and the pLI_intersect.bed file as the second
  argument (-b). The output is written to a file with the same name as the
  input file, with "_intersect.csv" appended to the end.
}
\examples{
  # Define a list of subdirectories
  subdir_list <- c("path/to/subdir1", "path/to/subdir2", "path/to/subdir3")
  
  # Run the intersect function
  intersect_func(subdir_list)
}

\author{Nermina Logo Lendo}