#Hyphen function

"""This is a function named add_hyphen that takes a character vector as input, 
and uses the stringr package function str_replace_all to replace all occurrences 
of a pattern in the input strings. The pattern to be replaced is a string that 
consists of one or more digits, followed by one or more uppercase letters. 
The replacement pattern is the original pattern with a hyphen - inserted between 
the digits and letters.

For example, if the input vector is c("A1B2C3", "D123EFG", "HIJKLM5N"), 
the output would be c("A1-B2-C3", "D123-EFG", "HIJKLM5-N").

This function is useful when you need to add hyphens to separate numeric and 
alphabetic characters in a string, which can be helpful for formatting or parsing data."""

add_hyphen <- function(x) {
  str_replace_all(x, "(\\d+)([A-Z]+)", "\\1-\\2")
}