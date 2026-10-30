library(dplyr)
library(readr)
library(stringr)

# Read the majority-vote classifications.
taxonomy <- read_tsv(
  snakemake@input[[1]],
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) |>
  rename(
    mag_id = `Genome ID`,
    gtdb_classification = `GTDB classification`,
    ncbi_classification = `Majority vote NCBI classification`
  )

stopifnot(!anyNA(taxonomy$mag_id), !anyDuplicated(taxonomy$mag_id))

# Extract the NCBI taxonomic ranks.
ranks <- c(
  domain = "d",
  phylum = "p",
  class = "c",
  order = "o",
  family = "f",
  genus = "g",
  species = "s"
)

for (rank in names(ranks)) {
  pattern <- paste0("(?:^|;)\\s*", ranks[[rank]], "__([^;]*)")

  taxonomy[[paste0("ncbi_", rank)]] <- taxonomy$ncbi_classification |>
    str_match(pattern) |>
    (\(x) x[, 2])() |>
    str_trim() |>
    na_if("")
}

# Select the most specific assigned taxon.
rank_columns <- paste0("ncbi_", names(ranks))

taxonomy$ncbi_most_specific <- do.call(
  coalesce,
  rev(taxonomy[rank_columns])
)

taxonomy |>
  arrange(mag_id) |>
  write_tsv(snakemake@output[[1]], na = "NA")
