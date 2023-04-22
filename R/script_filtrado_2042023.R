#Libraries:
library(dplyr)
library(stringr)


#####################GENE SELECTION########################################


# Reads the content of the file "ASC_Spain.PASS.annotated.txt" into a character
# vector "text", where each line of the file becomes an element of the vector
text <- readLines("home/nodotea/NERA/Exomas_HGUGM_autismo/")

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
max_cols <- max(lengths(sample_cols))
sample_df <- matrix("", nrow = max_cols, ncol = length(sample_cols))

for (i in seq_along(sample_cols)) {
  # pad the numeric part of each ID
  padded_ids <- gsub("(?<=-)(\\d)(?=-|$)", "00\\1", sample_cols[[i]], perl = TRUE)
  padded_ids <- gsub("(?<=-)(\\d{2})(?=-|$)", "0\\1", padded_ids, perl = TRUE)
  
  # assign padded IDs to the matrix
  sample_df[1:length(padded_ids), i] <- padded_ids
}

# combine the padded IDs back into a single column
data_df_ASD_genes$Padded_IDs <- apply(sample_df, 2, paste0, collapse = ",")


data_df_ASD_genes$Padded_IDs <- gsub(",+", ",", data_df_ASD_genes$Padded_IDs)
data_df_ASD_genes$Padded_IDs <- gsub("[:,]+$", "", data_df_ASD_genes$Padded_IDs)



# Loop through the IDs column of data_df
for (i in 1:nrow(data_df_ASD_genes)) {
  
  # Split IDs into a vector by ";"
  ids <- unlist(strsplit(as.character(data_df_ASD_genes$Padded_IDs[i]), ","))
  
  # Loop through the vector of IDs
  for (j in 1:length(ids)) {
    
    # Check if ID is missing "-HI", "-MA" or "-PA" substring
    if (grepl("^G01-GEA-\\d+", ids[j]) && !grepl("-HI|-MA|-PA", ids[j])) {
      
      # Keep only the first 14 characters of the ID
      ids[j] <- substr(ids[j], start = 1, stop = 14)
      
    } else if (grepl("^G01-GEA-\\d+", ids[j]) && grepl("-HI|-MA|-PA", ids[j])) {
      
      # Remove everything after "-HI", "-MA" or "-PA"
      ids[j] <- substr(ids[j], start = 1, stop = 14)
    }
  }
  
  # Combine modified IDs into a single string separated by ";"
  data_df_ASD_genes$Padded_IDs[i] <- paste(ids, collapse = ",")
}



###############################################################################

# Get unique substrings from the Padded_IDs column of the data_df_test data frame
unique_ids <- unique(str_extract(data_df_ASD_genes$Padded_IDs, "^G01-GEA-\\d+"))

# Create a new directory called "Split"
dir.create("Split", showWarnings = FALSE)

# Loop through unique IDs
for (id in unique_ids) {
  
  # Create dataframe containing only rows with current ID
  id_df <- data_df_ASD_genes[grepl(id, data_df_ASD_genes$Padded_IDs),]
  
  # Create folder for the ID if it doesn't already exist
  folder_name <- paste0("Split/", id)
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Write dataframe to file
  filename <- paste0(folder_name, "/", id, ".csv")
  write.csv(id_df, file = filename, row.names = FALSE)
}


# Define a function to remove rows that don't contain the ID substring + "-HI"
remove_non_HI_rows <- function(df, id) {
  df[grepl(paste0(id, "-HI"), df[, "Padded_IDs"]), ]
}

# Loop through each subdirectory
for (subdir in subdir_list) {
  # Get the unique ID from the subdirectory name
  id <- basename(subdir)
  
  # Get a list of all CSV files in the subdirectory
  file_list <- list.files(subdir, pattern = ".csv")
  
  # Loop through the file list and apply the function to each file
  for (file in file_list) {
    # Read in the CSV file as a data frame
    df <- read.csv(paste0(subdir, "/", file), header = TRUE, stringsAsFactors = FALSE)
    
    # Remove rows that don't contain the ID substring + "-HI" in the Padded_IDs column
    df <- remove_non_HI_rows(df, id)
    
    # Create a new file with the modified data frame
    write.csv(df, paste0(subdir, "/", "HI_only_", file), row.names = FALSE)
  }
}



#####################VARIANT FILTERING########################################
        # FILTERING CLINVAR PATHOGENIC VARIANTS

# Get a list of all subdirectories (i.e., unique IDs) in the directory
subdir_list <- list.dirs(dir_path, recursive = FALSE)

# Loop through each subdirectory
for (subdir in subdir_list) {
  
  # Get a list of all "HI_only_" files in the subdirectory
  file_list <- list.files(paste0(subdir, "/"), pattern = "^HI_only_", full.names = TRUE)
  
  # Loop through each "HI_only_" file in the subdirectory and apply the function
  for (file in file_list) {
    # Read in the CSV file as a data frame
    df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    
    # Filter rows that contain "Pathogenic" in the "ClinVar_Significance" column
    df_pathogenic <- subset(df, ClinVar_Significance == "Pathogenic")
    
    # Create a new file name without the "HI_only_" prefix
    new_filename <- gsub("^HI_only_", "", basename(file))
    new_filename <- paste0(subdir, "/", new_filename, "_path_var_clinvar.csv")
    
    # Write a new file with the filtered data frame and the new file name
    write.csv(df_pathogenic, file = new_filename, row.names = FALSE)
  }
}

      # FILTERING HOMOZYGOUS

# Loop through each subdirectory
for (subdir in subdir_list) {
  
  # Get a list of all "HI_only_" files in the subdirectory
  file_list <- list.files(paste0(subdir, "/"), pattern = "^HI_only_", full.names = TRUE)
  
  # Loop through each "HI_only_" file in the subdirectory and filter homozygous variants
  for (file in file_list) {
    # Read in the CSV file as a data frame
    df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    
    # Filter rows that have MaxPopFreq smaller than 0.05
    df <- subset(df, MaxPopFreq < 0.05)
    
    # Filter rows that do not have "AD" in the CGD_Inheritance column
    df <- df[!grepl("\\bAD\\b", df$CGD_Inheritance), ]
    
    # Filter rows that have a GATK.counts value smaller than 50
    df <- df[gsub("\\..*$", "", df$GATK.counts) < 50, ]
    
    # Filter rows that do not have "synonymous" in the Annotation.RefSeq column
    df <- df[!grepl("\\bsynonymous\\b", df$Annotation.RefSeq), ]
    
    # Create a new file with the filtered data frame
    filename <- paste0(subdir, "/", "HI_only_homozygous_", basename(file))
    write.csv(df, file = filename, row.names = FALSE)
  }
}




#FALTA CADD Y DANN. CADD corte en 15; DANN
#Mirar si falta añadir OMIM Disorder not empty y CGD category



        # FILTERING HETEROZYGOUS



