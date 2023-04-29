#' .
#'
#' This function takes a data frame and a specified ID as inputs. It removes any rows where the specified ID is present in the Padded_IDs column along with "-MA" or "-PA".
#'
#' @param df A data frame.
#' @param id A character string specifying the ID to be searched for.
#' 
#' @return A data frame with the specified ID present in the Padded_IDs column only if it is followed by "-HI".
#' 
#' @examples
#' # Given the following data frame:
#' df <- data.frame(Padded_IDs = c("G01-GEA-1234-MA,G01-GEA-1234-HI", "G01-GEA-5678-MA,G01-GEA-5678-PA,G01-GEA-5678-HI", "G01-GEA-9101-HI"))
#' 
#' # Keep only rows where ID "G01-GEA-1234" has "-HI" and not "-MA" or "-PA"
#' keep_de_novo(df, "G01-GEA-1234")
#'
#' @export
keep_de_novo <- function(df, id) {
  parent_ids <- paste0(id, "-MA|", id, "-PA")
  keep_rows <- sapply(strsplit(df$Padded_IDs, ","), function(x) {
    !any(grepl(parent_ids, x))
  })
  df[keep_rows, ]
