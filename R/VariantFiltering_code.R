#Libraries:
library(dplyr)
library(stringr)


#####################GENE SELECTION########################################

# Read file and save as data frame
Exomas_HGUGM_autismo <- read.delim("Exomas_HGUGM_autismo.tab")

Exomas_HGUGM_autismo <-  data.frame(lapply(Exomas_HGUGM_autismo, 
                   function(x) ifelse(is.na(x) | x == "","NA", x)))

# Prepend "chr" to the chromosome column
Exomas_HGUGM_autismo$Chr <- paste0("chr", Exomas_HGUGM_autismo$Chr)

# Convert start positions to numeric values and subtract 1 (bedtools requires 
# 0-based start positions)
Exomas_HGUGM_autismo$Start <- as.numeric(Exomas_HGUGM_autismo$Start)
Exomas_HGUGM_autismo$Start <- Exomas_HGUGM_autismo$Start - 1

# Write the modified data frame to a file in BED format
write.table(Exomas_HGUGM_autismo, "Exoma_HGUGM.bed", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Use bedtools to find the intersection of the input BED file with a gene 
# annotation file
system("bedtools intersect -a Exoma_HGUGM.bed -b GENE_BBDD.bed > intersect.bed")

# Read in the intersection file as a character vector
Exomas_ASD_autismo <- read.delim("intersect.bed", header=FALSE)

# Assign column names to the matrix using the original data frame
colnames(Exomas_ASD_autismo) <- colnames(Exomas_HGUGM_autismo)

# Remove the original data frame to save memory
rm(Exomas_HGUGM_autismo)

#####################DATA PREPROCESSING########################################

# Replace "|" with "," in GATK.samples column
# Replace "_" with "-" in GATK.samples column
Exomas_ASD_autismo$GATK.samples <- gsub("\\|", ",", 
                                       Exomas_ASD_autismo$GATK.samples)
Exomas_ASD_autismo$GATK.samples <- gsub("_", "-", 
                                       Exomas_ASD_autismo$GATK.samples)

# Split comma-separated values in GATK.samples into separate rows
# If an ID doesn't start with "G01-GEA-", add it to the beginning
Exomas_ASD_autismo <- Exomas_ASD_autismo %>%
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
Exomas_ASD_autismo <- Exomas_ASD_autismo %>%
  mutate(GATK.samples = str_pad(GATK.samples, width = 15, side = "right", 
                                pad = " ")) %>%
  mutate(GATK.samples = str_split(GATK.samples, ",")) %>%
  mutate(GATK.samples = lapply(GATK.samples, add_hyphen)) %>%
  mutate(GATK.samples = sapply(GATK.samples, function(x) paste(x, 
                                                               collapse = ","))) 

# Create a matrix of padded IDs from GATK.samples column
# Pad IDs with leading zeros if needed
sample_cols <- strsplit(Exomas_ASD_autismo$GATK.samples, ",")
max_cols <- max(lengths(sample_cols))
sample_df <- matrix("", nrow = max_cols, ncol = length(sample_cols))

for (i in seq_along(sample_cols)) {
  padded_ids <- gsub("(?<=-)(\\d)(?=-|$)", "00\\1", sample_cols[[i]], 
                     perl = TRUE)
  padded_ids <- gsub("(?<=-)(\\d{2})(?=-|$)", "0\\1", padded_ids, perl = TRUE)
  
  sample_df[1:length(padded_ids), i] <- padded_ids
}

# Combine padded IDs into a comma-separated string
Exomas_ASD_autismo$Padded_IDs <- apply(sample_df, 2, paste0, collapse = ",")

# Remove any consecutive commas or trailing commas
Exomas_ASD_autismo$Padded_IDs <- gsub(",+", ",", 
                                      Exomas_ASD_autismo$Padded_IDs)
Exomas_ASD_autismo$Padded_IDs <- gsub("[:,]+$", "", 
                                      Exomas_ASD_autismo$Padded_IDs)

Exomas_ASD_autismo$Padded_IDs_short <- ""

# Loop through each row in the dataframe 
for (i in 1:nrow(Exomas_ASD_autismo)) {
  
  # Split the Padded_IDs string into individual IDs and loop through each one
  ids <- unlist(strsplit(as.character(Exomas_ASD_autismo$Padded_IDs[i]), ","))
  
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
  Exomas_ASD_autismo$Padded_IDs_short[i] <- paste(ids, collapse = ", ")
}


############################CREATE TRIO FILES##################################

# Get unique substrings from the Padded_IDs column of the data_df_test dataframe
unique_ids <- unique(str_extract(Exomas_ASD_autismo$Padded_IDs, 
                                 "^G01-GEA-\\d+"))

# Create a new directory called "Split"
dir.create("Split", showWarnings = FALSE)

# Loop through unique IDs
for (id in unique_ids) {
  
  # Create dataframe containing only rows with current ID
  id_df <- Exomas_ASD_autismo[grepl(id, Exomas_ASD_autismo$Padded_IDs),]
  
  # Create folder for the ID if it doesn't already exist
  folder_name <- paste0("Split/", id)
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  
  # Write dataframe to file
  filename <- paste0(folder_name, "/", id, ".csv")
  write.csv(id_df, file = filename, row.names = FALSE)
}


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
  file_list <- list.files(paste0(subdir, "/"), pattern = "^HI_", 
                          full.names = TRUE)
  
  # Loop through each "HI_only_" file in the subdirectory and apply the function
  for (file in file_list) {
    # Read in the CSV file as a data frame
    df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
    
    # Filter rows that contain "Pathogenic" in the "ClinVar_Significance" column
    df_pathogenic <- subset(df, ClinVar_Significance == "Pathogenic")
    
    # Create a new file name without the "HI_only_" prefix
    new_filename <- gsub("^HI_", "", basename(file))
    new_filename <- paste0(subdir, "/", new_filename, "_path_var_clinvar.csv")
    
    # Write a new file with the filtered data frame and the new file name
    write.csv(df_pathogenic, file = new_filename, row.names = FALSE)
  }
}

#···················Filtering de novo pathogenic variants······················

filter_variants <- function(subdir_list, variant_type, maxpopfreq = NULL, 
                            MPCscore = NULL, filename_filter) {
  
  # Loop through each subdirectory
  for (subdir in subdir_list) {
    
    # Get a list of all "HI_" files in the subdirectory
    file_list <- list.files(paste0(subdir, "/"), pattern = "^HI_", 
                            full.names = TRUE)
    
    # Loop through each "HI_only_" file in the subdirectory and filter variants
    for (file in file_list) {
      # Read in the CSV file as a data frame
      df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
      
      # Filter for PTV variants
      if (variant_type == "PTV") {
        # Remove 3' and 5' variants, downstream, intronic, missense, and 
        # synonymous
        df <- subset(df, (grepl("frameshift", df$Annotation.RefSeq) | 
                            grepl("stop_gained", df$Annotation.RefSeq) | 
                            grepl("stop_lost", df$Annotation.RefSeq) |
                            grepl("start_lost", df$Annotation.RefSeq)))
      } 
      
      # Filter for missense variants
      if (variant_type == "missense") {
        # Keep only missense variants
        df <- subset(df, grepl("missense", df$Annotation.RefSeq))
      }
      
      # Filters the data frame to remove any variants with a maximum population 
      # frequency greater than the input value (if provided).
      if (!is.null(maxpopfreq)) {
        df <- subset(df, MaxPopFreq < maxpopfreq)
      }
      
      # Filters the data frame to remove any variants with a MPC score greater
      # than the input value (if provided).
      if(!is.null(MPCscore)) {
        df$MPC_score <- as.numeric(df$MPC_score)
        df <- subset(df, MPC_score >= MPCscore)
      }
      
      # Keeps only the rows where OMIM_ID is not NA
      df <- df[!is.na(df$OMIM_ID),]
      
      # Filters the data frame to remove any rows containing "MA" or "PA" in the 
      # "Padded_IDs_short" column.
      df <- subset(df, !grepl("MA|PA", Padded_IDs_short))
      
      # Creates a new file with the filtered data and writes it to disk. The 
      # new file is saved with a name based on the original file name and the 
      # input parameters.
      filename <- paste0(subdir, "/", "filtered_", filename_filter, "_", 
                         variant_type, "_", basename(file))
      write.csv(df, file = filename, row.names = FALSE)
    }
  }
}


# De novo PTV variants
filter_variants(subdir_list = subdir_list, variant_type = "PTV", 
                maxpopfreq = 0.01, gatkcounts = 25, filename_filter = "de_novo")

# De novo missense variants
filter_variants(subdir_list = subdir_list, variant_type = "missense", 
maxpopfreq = 0.01, gatkcounts = 25, MPCscore = 2, filename_filter =  "de_novo")

# Read in the data file
data_pLI <- read.delim("~/NERA/gnomad.v2.1.1.lof_metrics.by_gene.txt")

# Convert pLI column to numeric and filter based on pLI value
data_pLI$pLI <- as.numeric(data_pLI$pLI)
data_pLI_filtered <- subset(data_pLI, pLI >= 0.9)

# Subset data to only relevant columns and format chromosome column (chr, start,
# end)
data_pLI_filtered_intersectq <- data_pLI_filtered[,75:77]
# Format chromosome column (add "chr" prefix)
data_pLI_filtered_intersect$chromosome <- paste0("chr", 
                                        data_pLI_filtered_intersect$chromosome)


# Convert start_position column to numeric and subtract 1 from each value
data_pLI_filtered_intersect$start_position <- as.numeric(
  data_pLI_filtered_intersect$start_position)
data_pLI_filtered_intersect$start_position <- 
  data_pLI_filtered_intersect$start_position - 1

# Write intersect file to disk
write.table(data_pLI_filtered_intersect, "pLI_intersect.bed", sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Loop through each subdirectory
for (subdir in subdir_list) {
  # Copy the pLI_intersect.bed file to the current subdirectory
  file.copy("pLI_intersect.bed", paste0(subdir, "/"))
}

# Define a function to perform intersect for a list of subdirectories
intersect_func <- function(subdirlist) {
  
  # Loop through each subdirectory
  for (subdir in subdir_list) {
    
    # Find all files in subdirectory with a specific pattern
    file_list <- list.files(paste0(subdir, "/"),
                          pattern = "^filtered_de_novo_PTV_", full.names = TRUE)
    
    # Loop through each file and perform intersect
    for (file in file_list) {
      
      # Check if the file is empty
      if (file.info(file)$size > 0) {
        
        # Read in the file and write it as a bed file
        df <- read.csv(file, header = T, stringsAsFactors = FALSE)
        bed_file <- paste0(subdir, "/", 
                           tools::file_path_sans_ext(basename(file)), ".bed")
        write.table(df, file = bed_file, sep = "\t", quote = FALSE, 
                    col.names = FALSE, row.names = FALSE)
        
        # Perform intersect using bedtools and write output to file
        output_file <- paste0(bed_file, "_intersect.csv")
        system(paste("bedtools intersect -a", bed_file, 
                     "-b pLI_intersect.bed >", output_file))
      }
    }
  }
}


# Call intersect function with a list of subdirectories
intersect_func(subdirlist = subdir_list)
