
library(tidyverse)

coverm_files <- snakemake@input[["coverm"]]
output_file <- snakemake@output[["summary"]]

read_coverm <- function(file) {
  sample <- str_remove(basename(file), "\\.tsv$")

  read_tsv(file, show_col_types = FALSE) %>%
    rename(mag_id = Genome) %>%
    rename_with(~ str_remove(.x, "^.*?\\.fastq\\.gz "), -mag_id) %>%
    rename_with(~ str_to_lower(str_replace_all(.x, "[^A-Za-z0-9]+", "_"))) %>%
    rename_with(~ str_remove(.x, "_$")) %>%
    mutate(sample = sample, .before = 1)
}

summary <- map_dfr(coverm_files, read_coverm)

write_tsv(summary, output_file, na = "NA")

