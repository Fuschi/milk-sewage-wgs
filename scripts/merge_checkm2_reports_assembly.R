# Load necessary libraries
library(tidyverse)

# Define the input file paths and output path
input_files <- snakemake@params$reports
output_file <- snakemake@output[[1]]

# Function to process each CheckM2 file
process_checkm2_file <- function(file_path) {
  df <- read_tsv(file_path, col_types = cols())

  df_clean <- df %>%
    mutate(
      completeness = if_else(Completeness_Model_Used == "Neural Network (Specific Model)",
                             Completeness_Specific, Completeness_General)
    ) %>%
    select(genome = Name, completeness, contamination = Contamination)%>%
    mutate(genome = paste0(genome, ".fasta"))
  return(df_clean)
}

# Apply the function to all input files and combine the results
checkm2_data <- map_dfr(input_files, process_checkm2_file)

# Write the combined table to the output file as TSV (tab-separated)
write_csv(checkm2_data, output_file)

# Log success and list processed files
cat("CheckM2 table for dRep successfully created.\n")
cat("Processed files:\n")
print(input_files)

