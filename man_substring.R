"Function name: extract_substrings

Description: This function takes a string and a substring as input and returns 
a new string containing only the substrings that match the given input substring. 
It first splits the input string into a vector of individual samples using the 
comma as a separator. It then selects only the samples that start with the given 
substring by comparing the first 11 characters of each sample to the input substring. 
Finally, it combines the selected samples into a single string separated by commas.

Arguments:
  
  x: a string containing multiple substrings separated by commas
  substring: a string representing the substring to be extracted from the input string
Returns: a string containing only the substrings that match the given input substring, 
separated by commas. If no substrings match the given input substring, the function 
returns an empty string."


extract_substrings <- function(x, substring) {
  # Split the string into a vector of samples
  samples <- unlist(strsplit(x, ","))
  # Select only samples that start with the given substring
  samples_to_keep <- samples[substr(samples, 1, 11) == substring]
  # Combine selected samples into a single string
  return(paste(samples_to_keep, collapse = ","))
}