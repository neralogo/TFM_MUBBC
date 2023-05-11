\name{calc_percentage}
\alias{calc_percentage}
\title{Calculate the percentage of positive cases for groups based on patient 
data}
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
  The function first splits the patients into 7 groups based on the following 
  criteria:
  
  \begin{enumerate}
    \item Age at walking <= 12 and age at first birth <= 35 (Group1)
    \item Age at walking <= 10.7 and age at first birth > 35 (Group2)
    \item Age at walking > 12 and age at first birth <= 29.5 (Group3)
    \item Age at walking > 10.7 and age at first birth > 29.5 (Group4)
    \item Age at first birth <= 39.7 and history of epilepsy (Group5)
    \item Age at first birth <= 39.7 and no history of epilepsy (Group6)
    \item Age at first birth > 39.7 (Group7)
  \end{enumerate}
  
  The function then creates a new data frame with the patient groups and their
  exome data, and calculates the percentage of positive cases for each group.
  
  The output is a vector of percentages, one for each group.
}
\examples{
  # Create a sample data frame
  data <- data.frame(
    Age_walking = c(10.5, 12, 14, 11, 9),
    Age_F_Birth = c(32, 38, 25, 30, 40),
    Age_M_Birth = c(29, 25, 30, 34, 38),
    Epilepsy = c("No", "No", "No", "Yes", "Yes"),
    Exome = c("Positive", "Negative", "Positive", "Negative", "Positive")
  )
  
  # Calculate the percentage of positive cases for each group
  calc_percentage(data)
}
\seealso{
  \code{\link{aggregate}}, \code{\link{tidyr::complete}}
}
\author{
  Nermina Logo Lendo
}

