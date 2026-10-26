library(tidyverse)
library(readxl)

read_headers <- function(path) {
  con <- file(path, "rt")
  on.exit(close(con))
  headers <- character()
  repeat {
    lines <- readLines(con, n = 100000, warn = FALSE)
    if (!length(lines)) break
    headers <- c(headers, lines[startsWith(lines, ">")])
  }
  str_remove(headers, "^>")
}

collapse_values <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) str_c(x, collapse = " | ") else NA_character_
}

aggregate_vfdb <- function(files, fasta, annotations, abundance_output, metadata_output) {
  abundance <- map(files, function(f) {
    df <- read_tsv(f, col_types = cols(reference = col_character(), .default = col_double()), show_col_types = FALSE)
    stop_for_problems(df)
    df |>
      mutate(sample_id = basename(dirname(f)), gene_id = str_extract(reference, "^VFG\\d+")) |>
      select(sample_id, gene_id, -reference, everything())
  }) |>
    list_rbind() |>
    arrange(sample_id, gene_id)

  headers <- read_headers(fasta)
  metadata <- tibble(
    gene_id = str_extract(headers, "^VFG\\d+"),
    accession = str_match(headers, "\\(gb\\|([^)]*)\\)")[, 2],
    gene_name = str_match(headers, "\\)\\s+\\(([^)]*)\\)")[, 2],
    vf_id = str_extract(headers, "\\bVF\\d+\\b"),
    vfc_id = str_extract(headers, "\\bVFC\\d+\\b"),
    organism = str_match(headers, "\\[([^\\[\\]]+)\\]$")[, 2],
    reference = headers
  )

  factors <- read_excel(annotations, skip = 1, col_types = "text") |>
    rename(vf_id = VFID, vf_name = VF_Name, vf_full_name = VF_FullName, bacteria = Bacteria,
           annotation_vfc_id = VFCID, vf_category = VFcategory, characteristics = Characteristics,
           structure = Structure, function_description = Function, mechanism = Mechanism, literature = Reference) |>
    filter(!is.na(vf_id)) |>
    group_by(vf_id) |>
    summarise(across(everything(), collapse_values), .groups = "drop")

  stopifnot(!anyNA(abundance$gene_id), !anyNA(metadata$gene_id), !anyDuplicated(metadata$gene_id))
  stopifnot(all(abundance$gene_id %in% metadata$gene_id))

  metadata <- metadata |>
    left_join(factors, by = "vf_id", relationship = "many-to-one") |>
    arrange(gene_id)

  write_tsv(abundance, abundance_output, na = "NA")
  write_tsv(metadata, metadata_output, na = "NA")
}

aggregate_vfdb(snakemake@input[["quantifications"]], snakemake@input[["fasta"]], snakemake@input[["annotations"]],
               snakemake@output[["abundance"]], snakemake@output[["metadata"]])
