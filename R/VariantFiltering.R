#Libraries:
library(dplyr)
library(stringr)


#####################GENE SELECTION########################################

# Read in file as character vector
text <- readLines("NERA/Exomas_HGUGM_autismo/Exomas_HGUGM_autismo.tab.gz")

# Determine the maximum number of columns in the text data
max_cols <- max(sapply(strsplit(text, "\t"), length))

# Create an empty matrix with the correct dimensions
data_matrix <- matrix(nrow = length(text), ncol = max_cols)

# Fill the matrix with data from the text file
for (i in 1:length(text)) {
  line <- strsplit(text[i], "\t")[[1]]
  data_matrix[i, 1:length(line)] <- line
}

# Assign column names to the matrix using the first row
colnames(data_matrix) <- data_matrix[1, ]

# Remove the first row of the matrix (since it contains column names)
data_matrix <- data_matrix[-1,]

# Convert the matrix to a data frame
data_df <- as.data.frame(data_matrix)

# Remove the matrix to save memory
rm(data_matrix)

# Replace empty values and NAs with "NA" in the data frame
data_df_NA <- data.frame(lapply(data_df, function(x) ifelse(is.na(x) | x == "",
                                                            "NA", x)))

# Prepend "chr" to the chromosome column
data_df_NA$Chr <- paste0("chr", data_df_NA$Chr)

# Convert start positions to numeric values and subtract 1 (bedtools requires 
# 0-based start positions)
data_df_NA$Start <- as.numeric(data_df_NA$Start)
data_df_NA$Start <- data_df_NA$Start - 1

# Write the modified data frame to a file in BED format
write.table(data_df_NA, "Exoma_HGUGM_NA.bed", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Use bedtools to find the intersection of the input BED file with a gene 
# annotation file
system("bedtools intersect -a Exoma_HGUGM_NA.bed -b GENE_BBDD.bed > intersect.bed")

# Read in the intersection file as a character vector
text <- readLines("intersect.bed")

# Determine the maximum number of columns in the text data
max_cols <- max(sapply(strsplit(text, "\t"), length))

# Create an empty matrix with the correct dimensions
data_matrix <- matrix(nrow = length(text), ncol = max_cols)

# Fill the matrix with data from the text file
for (i in 1:length(text)) {
  line <- strsplit(text[i], "\t")[[1]]
  data_matrix[i, 1:length(line)] <- line
}

# Assign column names to the matrix using the original data frame
colnames(data_matrix) <- colnames(data_df)

# Remove the first row of the matrix (since it contains column names)
data_matrix <- data_matrix[-1,]

# Convert the matrix to a data frame
data_df_ASD_genes <- as.data.frame(data_matrix)

# Remove the matrix to save memory
rm(data_matrix)

# Remove the original data frame to save memory
rm(data_df)

# Remove the NA filled dataframe to save memory
rm(data_df_NA)

#####################DATA PREPROCESSING########################################

# Replace "|" with "," in GATK.samples column
# Replace "_" with "-" in GATK.samples column
data_df_ASD_genes$GATK.samples <- gsub("\\|", ",", 
                                       data_df_ASD_genes$GATK.samples)
data_df_ASD_genes$GATK.samples <- gsub("_", "-", 
                                       data_df_ASD_genes$GATK.samples)

# Split comma-separated values in GATK.samples into separate rows
# If an ID doesn't start with "G01-GEA-", add it to the beginning
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

# Add hyphen between numbers and letters in GATK.samples column
add_hyphen <- function(x) {
  str_replace_all(x, "(\\d+)([A-Z]+)", "\\1-\\2")
}

# Pad the GATK samples column with spaces to ensure consistent length.
# Split the GATK samples column by comma into a list.
# Apply the add_hyphen function to add hyphens between the numbers and letters 
# in each sample ID.
# Join the list back into a comma-separated string for each row.
data_df_ASD_genes <- data_df_ASD_genes %>%
  mutate(GATK.samples = str_pad(GATK.samples, width = 15, side = "right", 
                                pad = " ")) %>%
  mutate(GATK.samples = str_split(GATK.samples, ",")) %>%
  mutate(GATK.samples = lapply(GATK.samples, add_hyphen)) %>%
  mutate(GATK.samples = sapply(GATK.samples, function(x) paste(x, 
                                                               collapse = ","))) 

# Create a matrix of padded IDs from GATK.samples column
# Pad IDs with leading zeros if needed
sample_cols <- strsplit(data_df_ASD_genes$GATK.samples, ",")
max_cols <- max(lengths(sample_cols))
sample_df <- matrix("", nrow = max_cols, ncol = length(sample_cols))

for (i in seq_along(sample_cols)) {
  padded_ids <- gsub("(?<=-)(\\d)(?=-|$)", "00\\1", sample_cols[[i]], 
                     perl = TRUE)
  padded_ids <- gsub("(?<=-)(\\d{2})(?=-|$)", "0\\1", padded_ids, perl = TRUE)
  
  sample_df[1:length(padded_ids), i] <- padded_ids
}

# Combine padded IDs into a comma-separated string
data_df_ASD_genes$Padded_IDs <- apply(sample_df, 2, paste0, collapse = ",")

# Remove any consecutive commas or trailing commas
data_df_ASD_genes$Padded_IDs <- gsub(",+", ",", data_df_ASD_genes$Padded_IDs)
data_df_ASD_genes$Padded_IDs <- gsub("[:,]+$", "", data_df_ASD_genes$Padded_IDs)

data_df_ASD_genes$Padded_IDs_short <- ""

# Loop through each row in the dataframe 
for (i in 1:nrow(data_df_ASD_genes)) {
  
  # Split the Padded_IDs string into individual IDs and loop through each one
  ids <- unlist(strsplit(as.character(data_df_ASD_genes$Padded_IDs[i]), ","))
  
  for (j in 1:length(ids)) {
    
    # If the ID starts with "G01-GEA-" and does not contain "-HI", "-MA", or 
    # "-PA" then add "-HI" and keep only the first 14 characters
    if (grepl("^G01-GEA-\\d{3}", ids[j]) && !grepl("-MA|-PA|-HI", ids[j])) {
      ids[j] <- paste0(substr(ids[j], start = 1, stop = 11), "-HI")
      
      # If the ID starts with "G01-GEA-" and contains "-HI", "-MA", or "-PA"
      # then keep only the first 14 characters
    } else if (grepl("^G01-GEA-\\d+", ids[j]) && grepl("-HI|-MA|-PA", ids[j])) {
      
      ids[j] <- substr(ids[j], start = 1, stop = 14)
    }
  }
  # Update the Padded_IDs column with the modified IDs
  data_df_ASD_genes$Padded_IDs_short[i] <- paste(ids, collapse = ", ")
}


############################CREATE TRIO FILES##################################

# Get unique substrings from the Padded_IDs column of the data_df_test dataframe
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
#remove_non_HI_rows <- function(df, id) {
#  df[grepl(paste0(id, "-HI"), df[, "Padded_IDs"]), ]
#}

##################seguir aqui##############!!!!!!!!!!!!!!!!!!!!!!!

# Define a function to remove rows that contain the parents ID substring so we 
# only keep de novo mutations
keep_de_novo <- function(df, id) {
  parent_ids <- paste0(id, "-MA|", id, "-PA")
  keep_rows <- sapply(strsplit(df$Padded_IDs, ","), function(x) {
    !any(grepl(parent_ids, x))
  })
  df[keep_rows, ]
}

dir_path <- "Split/"

# Get a list of all subdirectories (i.e., unique IDs) in the directory
subdir_list <- list.dirs(dir_path, recursive = FALSE)

# Loop through each subdirectory
for (subdir in subdir_list) {
  # Get the unique ID from the subdirectory name
  id <- basename(subdir)
  
  # Get a list of all CSV files in the subdirectory
  file_list <- list.files(subdir, pattern = ".csv")
  
  # Loop through the file list and apply the function to each file
  for (file in file_list) {
    # Read in the CSV file as a data frame
    df <- read.csv(paste0(subdir, "/", file), header = TRUE, 
                   stringsAsFactors = FALSE)
    
    # Remove rows that don't contain the ID substring + "-HI" in the Padded_IDs 
    # column
    df <- keep_de_novo(df, id)
    
    # Create a new file with the modified data frame
    write.csv(df, paste0(subdir, "/", "HI_", file), row.names = FALSE)
  }
}



#####################VARIANT FILTERING########################################

#··················Filtering CLINVAR pathogenic variants······················

# Loop through each subdirectory
for (subdir in subdir_list) {
  
  # Get a list of all "HI_only_" files in the subdirectory
  file_list <- list.files(paste0(subdir, "/"), pattern = "^HI_only_", 
                          full.names = TRUE)
  
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

#············Filtering Homozygous and Heterozygous pathogenic variants··········

filter_variants <- function(subdir_list, variant_type, maxpopfreq = NULL, 
                            inheritance = NULL, gatkcounts = NULL, 
                            filename_hom_het) {
  
  # Loop through each subdirectory
  for (subdir in subdir_list) {
    
    # Get a list of all "HI_only_" files in the subdirectory
    file_list <- list.files(paste0(subdir, "/"), pattern = "^HI_only_", 
                            full.names = TRUE)
    
    # Loop through each "HI_only_" file in the subdirectory and filter variants
    for (file in file_list) {
      # Read in the CSV file as a data frame
      df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
      
      # Filter for PTV variants
      if (variant_type == "PTV") {
        # Remove 3' and 5' variants, downstream, intronic, missense, and synonymous
        df <- subset(df, !(grepl("3'", df$Annotation.RefSeq) | 
                             grepl("5'", df$Annotation.RefSeq) | 
                             grepl("downstream", df$Annotation.RefSeq) |
                             grepl("intronic", df$Annotation.RefSeq) |
                             grepl("missense", df$Annotation.RefSeq) |
                             grepl("synonymous", df$Annotation.RefSeq)))
      } 
      
      # Filter for missense variants
      if (variant_type == "missense") {
        # Keep only missense variants
        df <- subset(df, grepl("missense", df$Annotation.RefSeq))
      }
      
      # Filters the data frame to remove any variants with a maximum population 
      # frequency greater than the input value (if provided).
      if (!is.null(maxpopfreq)) {
        df <- subset(df, MaxPopFreq < hom_maxpopfreq)
      }
      
      # Filters the data frame to remove any variants with an inheritance 
      # pattern matching the input regular expression (if provided).
      if (!is.null(inheritance)) {
        df <- df[!grepl(hom_inheritance, df$CGD_Inheritance), ]
      }
      
      # Filters the data frame to remove any variants with a GATK count greater 
      # than the input value (if provided).
      if (!is.null(gatkcounts)) {
        df <- df[gsub("\\..*$", "", df$GATK.counts) < hom_gatkcounts, ]
      }
      
      # Creates a new file with the filtered data and writes it to disk. The 
      # new file is saved with a name based on the original file name and the 
      # input parameters.
      filename <- paste0(subdir, "/", "filtered_", filename_hom_het, "_", 
                         variant_type, "_", basename(file))
      write.csv(df, file = filename, row.names = FALSE)
    }
  }
}


# Homozygote PTV variants
filter_variants(variant_type = "PTV", maxpopfreq = 0.05, inheritance = "AD",
                gatkcounts = 50, filename_hom_het = "HOM")

# Homozygote missense variants
filter_variants(variant_type = "missense", maxpopfreq = 0.05, inheritance = "AD",
                gatkcounts = 50, filename_hom_het = "HOM")

# Compound heterozygote missense variants
#filter_variants(variant_type = "missense", maxpopfreq = 0.05, inheritance = "AD",
#               gatkcounts = 50, filename_hom_het = "HOM")

# Heterozygote with AD inheritance missense variants
filter_variants(variant_type = "PTV", maxpopfreq = 0.01, inheritance = "AR",
                gatkcounts = 25, filename_hom_het = "HET")

# Heterozygote with AD inheritance missense variants
filter_variants(variant_type = "missense", maxpopfreq = 0.01, inheritance = "AR",
                gatkcounts = 25, filename_hom_het = "HOM")

# Filter out the PTV whose gene pLI score is lower than 0.9

