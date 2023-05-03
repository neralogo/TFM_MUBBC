\name{keep_de_novo}
\alias{keep_de_novo}
\title{
Keep de novo variants
}
\usage{
keep_de_novo(df, id)
}
\arguments{
\item{df}{A data frame containing the de novo mutations and parent IDs}
\item{id}{A character string specifying the ID to be searched for}
}
\description{
This function removes rows from a data frame containing parent IDs and
returns a subset that only includes de novo mutations. It searches for
rows where the parent IDs match those of the child and removes them.
}
\details{
The function takes in a data frame that contains columns for the ID of
the child, the ID of the parents, and information about mutations. The
function looks for parent IDs that match the ID of the child and removes
the rows that contain these parent IDs. This is because the purpose of
this function is to filter out parent mutations and return only de novo
mutations.
}
\value{
A data frame containing only de novo mutations, i.e. mutations where
neither parent has the mutation.
}
\examples{

df <- data.frame(IDs = c("G01-GEA-1234-MA,G01-GEA-1234-HI", "G01-GEA-5678-MA,G01-GEA-5678-PA,G01-GEA-5678-HI", "G01-GEA-9101-HI"))
# Keep only rows where ID "G01-GEA-1234" has "-HI" and not "-MA" or "-PA"
keep_de_novo(df, "G01-GEA-1234")
}
\author{Nermina Logo Lendo}
\keyword{data}
