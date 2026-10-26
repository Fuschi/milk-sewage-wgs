library(tidyverse)

collapse_values <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) str_c(x, collapse = " | ") else NA_character_
}

aggregate_resfinder <- function(files, phenotypes, abundance_output, metadata_output) {
  abundance <- map(files, function(f) {
    df <- read_tsv(f, col_types = cols(reference = col_character(), .default = col_double()), show_col_types = FALSE)
    stop_for_problems(df)
    df |> mutate(sample_id = basename(dirname(f)), .before = 1)
  }) |>
    list_rbind() |>
    rename(gene_id = reference) |>
    arrange(sample_id, gene_id)

  metadata <- read_tsv(phenotypes, col_types = cols(.default = col_character())) |>
    rename(gene_id = `Gene_accession no.`, antibiotic_class = Class, phenotype = Phenotype, pmid = PMID,
           resistance_mechanism = `Mechanism of resistance`, notes = Notes, required_gene = Required_gene) |>
    group_by(gene_id) |>
    summarise(across(everything(), collapse_values), .groups = "drop")

  metadata <- tibble(gene_id = union(metadata$gene_id, abundance$gene_id)) |>
    left_join(metadata, by = "gene_id", relationship = "one-to-one") |>
    arrange(gene_id)

  write_tsv(abundance, abundance_output, na = "NA")
  write_tsv(metadata, metadata_output, na = "NA")
}

aggregate_resfinder(snakemake@input[["quantifications"]], snakemake@input[["phenotypes"]],
                   snakemake@output[["abundance"]], snakemake@output[["metadata"]])

