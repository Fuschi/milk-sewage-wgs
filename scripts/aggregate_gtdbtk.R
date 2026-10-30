library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)

taxonomy_files <- unlist(snakemake@input, use.names = FALSE)
output_file <- snakemake@output[["summary"]]

taxonomy <- map_dfr(
  taxonomy_files,
  ~ read_tsv(
    .x,
    na = c("", "N/A"),
    show_col_types = FALSE
  )
) |>
  rename(mag_id = user_genome) |>
  separate_wider_delim(
    classification,
    delim = ";",
    names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    cols_remove = FALSE
  ) |>
  mutate(
    across(
      domain:species,
      ~ str_remove(.x, "^[dpcofgs]__") |>
        na_if("")
    )
  ) |>
  relocate(domain:species, .after = mag_id) |>
  arrange(mag_id)

write_tsv(taxonomy, output_file, na = "NA")
