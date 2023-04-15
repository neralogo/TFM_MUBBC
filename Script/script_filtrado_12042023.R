# Set the working directory to the specified path
setwd("~/NERA/TFM/SCRIPT")

#Libraries:
library(dplyr)
library(stringr)


#####################GENE SELECTION########################################


# Reads the content of the file "ASC_Spain.PASS.annotated.txt" into a character
# vector "text", where each line of the file becomes an element of the vector
text <- readLines("EXOMAS_HGUGM_AUTISMO1Y2/Exomas_HGUGM_autismo.tab.gz")

# Calculates the maximum number of columns in the tab-separated values file by 
# splitting each line of the file by tabs and finding the maximum number of 
# resulting elements
max_cols <- max(sapply(strsplit(text, "\t"), length))

# Creates an empty matrix with the same number of rows as "text" and the same 
# number of columns as "max_cols"
data_matrix <- matrix(nrow = length(text), ncol = max_cols)

# Loops through each element of "text" and splits it by tabs to create a vector 
# "line". The values of "line" are then assigned to the corresponding row of 
# "data_matrix"
for (i in 1:length(text)) {
  line <- strsplit(text[i], "\t")[[1]]
  data_matrix[i, 1:length(line)] <- line}

# Assigns the values of the first row of "data_matrix" as column names of the 
# matrix. The first row is then removed from "data_matrix", and the matrix is 
# converted to a data frame. Finally, the matrix is removed from the environment.
colnames(data_matrix) <- data_matrix[1, ]
data_matrix <- data_matrix[-1,]
data_df <- as.data.frame(data_matrix)
rm(data_matrix)


data_df_NA <- data.frame(lapply(data_df, function(x) ifelse(is.na(x) | x == "", "NA", x)))
data_df_NA$Chr <- paste0("chr", data_df_NA$Chr)
data_df_NA$Start <- as.numeric(data_df_NA$Start)
data_df_NA$Start <- data_df_NA$Start - 1


write.table(data_df_NA, "Exoma_HGUGM_NA.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


system("bedtools intersect -a Exoma_HGUGM_NA.bed -b GENE_BBDD.bed > intersect.bed")

text <- readLines("intersect.bed")
max_cols <- max(sapply(strsplit(text, "\t"), length))
data_matrix <- matrix(nrow = length(text), ncol = max_cols)
for (i in 1:length(text)) {
  line <- strsplit(text[i], "\t")[[1]]
  data_matrix[i, 1:length(line)] <- line}
colnames(data_matrix) <- colnames(data_df)
data_matrix <- data_matrix[-1,]
data_df_ASD_genes <- as.data.frame(data_matrix)
rm(data_matrix)
rm(data_df)


#####################DATA PREPROCESSING########################################


# Replaces the "|" character with "," in the "GATK.samples" column of "data_df_test"
# and splits each element of the column by "," to create a list.
data_df_ASD_genes$GATK.samples <- gsub("\\|", ",", data_df_ASD_genes$GATK.samples)
data_df_ASD_genes$GATK.samples <- gsub("_", "-", data_df_ASD_genes$GATK.samples)

#Loops through each line and each ID of the row, checks whether the ID already
#contains the substring "G01-GEA-"; if not it replaces "GEA" with "G01-GEA-".
#Finally, it will join the modified IDs back together and collapse them with a comma.
data_df_ASD_genes <- data_df_ASD_genes %>%
  mutate(GATK.samples = strsplit(as.character(GATK.samples), ",")) %>%
  mutate(GATK.samples = lapply(GATK.samples, function(ids) {
    sapply(ids, function(id) {
      if (!grepl("G01-GEA-", id)) {
        id <- gsub("^GEA", "G01-GEA-", id)
      }
      id
    })
  })) %>%
  mutate(GATK.samples = sapply(GATK.samples, function(ids) {
    paste(ids, collapse = ",")
  }))

# Define function to add hyphen after the number in the ID
add_hyphen <- function(x) {
  str_replace_all(x, "(\\d+)([A-Z]+)", "\\1-\\2")
}

# Homogenize the sample names
data_df_ASD_genes <- data_df_ASD_genes %>%
  mutate(GATK.samples = str_pad(GATK.samples, width = 15, side = "right", pad = " ")) %>%
  mutate(GATK.samples = str_split(GATK.samples, ",")) %>%
  mutate(GATK.samples = lapply(GATK.samples, add_hyphen)) %>%
  mutate(GATK.samples = sapply(GATK.samples, function(x) paste(x, collapse = ","))) 

# Split the sample IDs into separate columns
sample_cols <- strsplit(data_df_ASD_genes$GATK.samples, ",")
max_cols <- max(lengths(data_df_ASD_genes))
sample_df <- matrix(NA, nrow = max_cols, ncol = length(sample_cols))
for (i in seq_along(sample_cols)) {
  sample_df[1:length(sample_cols[[i]]), i] <- sample_cols[[i]]
}

sample_df <- gsub("(?<=-)(\\d)(?=-|$)", "00\\1", sample_df, perl = TRUE)
sample_df <- gsub("(?<=-)(\\d{2})(?=-|$)", "0\\1", sample_df, perl = TRUE)

dataframe_concat <- apply(sample_df, 2, paste0, collapse = ",")
data_df_ASD_genes$Padded_IDs <- dataframe_concat


#Remove all "NA" instances from rows of Padded IDs
data_df_ASD_genes$Padded_IDs <- gsub("NA,", "", data_df_ASD_genes$Padded_IDs)
data_df_ASD_genes$Padded_IDs <- gsub(",NA", "", data_df_ASD_genes$Padded_IDs)

# Loop through the IDs column of data_df
for (i in 1:nrow(data_df_ASD_genes)) {
  
  # Split IDs into a vector by ";"
  ids <- unlist(strsplit(as.character(data_df_ASD_genes$Padded_IDs[i]), ","))
  
  # Loop through the vector of IDs
  for (j in 1:length(ids)) {
    
    # Check if ID is missing "-HI" substring
    if (!grepl("-HI", ids[j]) && !grepl("-MA|-PA", ids[j]) && grepl("^G01-GEA-\\d+", ids[j])) {
      
      # Add "-HI" substring to ID
      ids[j] <- sub("(^G01-GEA-\\d{3})([^;]+)(;.*)", "\\1-HI\\2\\3", ids[j])
    }
  }
  
  # Combine modified IDs into a single string separated by ";"
  data_df_ASD_genes$Padded_IDs[i] <- paste(ids, collapse = ",")
}






###############################################################################


# Define function to extract substrings from a string
extract_substrings <- function(x, substring) {
  # Split the string into a vector of samples
  samples <- unlist(strsplit(x, ","))
  # Select only samples that start with the given substring
  samples_to_keep <- samples[substr(samples, 1, 11) == substring]
  # Combine selected samples into a single string
  return(paste(samples_to_keep, collapse = ","))
}

# Get unique substrings from the Padded_IDs column of the data_df_test data frame
unique_substrings <- unique(substr(unlist(strsplit(data_df_test$Padded_IDs, ",")), 1, 11))

# Create a list of data frames, with each data frame containing rows from data_df_test where the Padded_IDs column contains a particular substring
substring_df_list <- lapply(unique_substrings, function(substring) {
  rows <- which(sapply(data_df_test$Padded_IDs, function(x) any(grepl(substring, x))))
  substring_df <- data_df_test[rows, ]
  substring_df$Padded_IDs <- sapply(substring_df$Padded_IDs, extract_substrings, substring = substring)
  return(substring_df)
})

# Create a new directory called "Split"
dir.create("Split", showWarnings = FALSE)

# Write each data frame in substring_df_list to a CSV file in the "Split" directory
for (i in seq_along(substring_df_list)) {
  substring <- unique_substrings[i]
  filename <- paste0("Split/substring_", substring, ".csv")
  write.csv(substring_df_list[[i]], file = filename, row.names = FALSE)
}



# Set the directory path where your Excel files are located
dir_path <- "Split/"

# Get a list of all Excel files in the directory
file_list <- list.files(dir_path, pattern = ".csv")

# Create a function to remove rows that don't contain "PA" in the 8th column
remove_non_HI_rows <- function(df) {
  df[grepl("HI", df[,8]), ]
}

# Use lapply and a for loop to iterate through all Excel files
for (file in file_list) {
  # Read in the Excel file
  df <- read.csv(paste0(dir_path, "/", file))
  
  # Remove rows that don't contain "HI" in the 8th column
  df_HI <- remove_non_HI_rows(df)
  
  # Get the number from the 8th column
  sample_id <- unique(substring(df_HI[, 132], 9, 11))
  
  # Create the new filename
  new_file <- paste0("Split_", sample_id, ".csv")
  
  # Save the modified Excel file with the new filename
  write.csv(df_HI, paste0(dir_path, "/", new_file), row.names = FALSE)
}


#####################VARIANT FILTERING########################################

#Filtering homozygous 

#filtered_df_test <- subset(data_df_test, Zygosity = "HOMOCIGOTO") ADAPTARLO A LA PALABRA QUE SALGA QUE NO SÉ CUAL ES EXACTAMENTE

#Filtering maximum frequency on any population HAY QUE MODIFICARLO PARA QUE LO HAGA DE CADA ARCHIVO
filtered_df_test <- subset(filtered_df_test, MaxPopFreq <= 0.05)
"filtered_df_test <- subset(filtered_df_test, MaxPopSemiFreq <= 0.05)" #No ejecutamos porque no está en este archivo

#Filtering out the variants with no OMIM fenotype
filtered_df_test <- subset(filtered_df_test, OMIM_Disorder != "")

#Filtering out all results non related to neurodeveolpmental disorders
filtered_df_test <- subset(filtered_df_test, CGD_Category != "Neuro") #Hay que ajustarlo para que sea cualquier palabra que contenga "Neuro"

#Filtering out all autosomal dominant variants 
filtered_df_test <- subset(filtered_df_test, CGD_Inheritance != "AD") #Comprobar como aparece autosómico dominante

#Remova all variants that appear in more than 50 subjects (internal DB)
filtered_df_test <- subset(filtered_df_test, vardb_Illumina_GATK <= 50)

#Remove all variants that we are not interested in (synonymous)
filtered_df_test <- subset(filtered_df_test, Annotation.RefSeq != "synonymous") #Determinar qué variantes son de interés para nosotros 


#FALTA CADD Y DANN. CADD corte en 15; DANN