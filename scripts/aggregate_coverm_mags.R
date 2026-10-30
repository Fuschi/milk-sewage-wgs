library(dplyr)
library(purrr)
library(readr)
library(stringr)

coverm_files <- snakemake@input[["coverm"]]
output_file <- snakemake@output[["summary"]]

read_coverm <- function(file) {
  sample <- str_remove(basename(file), "\\.tsv$")

  read_tsv(file, show_col_types = FALSE) |>
    rename(mag_id = Genome) |>
    rename_with(
      ~ str_remove(.x, "^.*?\\.fastq\\.gz "),
      -mag_id
    ) |>
    rename(
      relative_abundance = `Relative Abundance (%)`,
      mean = Mean,
      trimmed_mean = `Trimmed Mean`,
      covered_fraction = `Covered Fraction`,
      covered_bases = `Covered Bases`,
      variance = Variance,
      rpkm = RPKM,
      tpm = TPM
    ) |>
    mutate(sample = sample, .before = 1)
}

summary <- map_dfr(coverm_files, read_coverm)

write_tsv(summary, output_file, na = "NA")
