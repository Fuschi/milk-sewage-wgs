# Load necessary libraries
library(tidyverse)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

df <- read_tsv(input_file, col_types = cols())

df_clean <- df %>%
  mutate(
    completeness = if_else(Completeness_Model_Used == "Neural Network (Specific Model)",
                           Completeness_Specific, Completeness_General)
  ) %>%
  select(genome = Name, completeness, contamination = Contamination) %>%
  mutate(genome = paste0(genome, ".fasta"))

write_csv(df_clean, output_file)

cat("Single CheckM2 report processed for dRep.\n")
cat("Processed file:\n")
print(input_file)


