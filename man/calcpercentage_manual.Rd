\name{calc_percentage}
\alias{calc_percentage}
\title{Calculate the percentage of positive cases for groups based on patient data}
\usage{
  calc_percentage(data)
}
\arguments{
  \item{data}{a data frame containing patient data}
}
\description{
  This function splits the patients into 7 groups based on their age at walking,
  age at first birth, and history of epilepsy, and calculates the percentage of
  positive cases for each group based on their exome data.
}
\details{
  The function first splits the patients into 7 groups based on different split
  criteria.
  
  The function then creates a new data frame with the patient groups and their
  exome data, and calculates the percentage of positive cases for each group.
  
  The output is a vector of percentages, one for each group.
}

\seealso{
  \code{\link{aggregate}}, \code{\link{tidyr::complete}}
}
\author{
  Nermina Logo Lendo
}

