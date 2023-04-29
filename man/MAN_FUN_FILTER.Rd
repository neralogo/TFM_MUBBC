\name{filter_variants}
\alias{filter_variants}
\title{Filter \emph{De Novo} Variants in Annotated Dataframes

}
\description{This function takes several dataframes within a subdirectories list with annotated data and filters \emph{de novo} variants based on the variant type, maximum population frequency, inheritance pattern, GATK count, and output filename.

}
\usage{
filter_variants(subdir_list, variant_type, maxpopfreq = NULL, inheritance = NULL, gatkcounts = NULL, filename_hom_het)
}
}

\arguments{
  \item{subdir_list}{ A list of subdirectory paths where the input CSV files are located.
}
  \item{variant_type}{ A character string indicating the variant type to filter for. Possible values are "PTV" (for protein-truncating variants) or "missense" (for missense variants).

}
  \item{maxpopfreq}{A numeric value indicating the maximum population frequency allowed for a variant to pass the filter. Variants with a population frequency greater than this value will be filtered out. Default is \code{NULL}.
}
  \item{inheritance}{A character string indicating a regular expression pattern used to filter variants based on inheritance pattern. Variants with an inheritance pattern matching this pattern will be filtered out. Default is \code{NULL}.
}
  \item{gatkcounts}{An integer indicating the maximum GATK count allowed for a variant to pass the filter. Variants with a GATK count greater than this value will be filtered out. Default is \code{NULL}.
}
  \item{filename_hom_het}{A character string indicating the filename prefix for the output files. This is used to differentiate files based on homozygous or heterozygous calls.}
}


\author{Nermina Logo Lendo
}

\examples{
Define input arguments
subdir_list <- list("./data/subdir1", "./data/subdir2")
variant_type <- "PTV"
maxpopfreq <- 0.01
inheritance <- "^de novo$"
gatkcounts <- 10
filename_hom_het <- "het"

{filter_variants(subdir_list, variant_type, maxpopfreq, inheritance, gatkcounts, filename_hom_het)}
}

\keyword{filtering}

