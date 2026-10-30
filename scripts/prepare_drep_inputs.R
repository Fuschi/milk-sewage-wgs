#!/usr/bin/env Rscript

library(tidyverse)

# Selection thresholds
min_completeness <- 50
max_contamination <- 10

# Find DAS Tool bin FASTA files
bin_files <- list.files("data/refinement/dastool", pattern = "\\.fa$",
                        recursive = TRUE, full.names = TRUE)

bins <- tibble(
  fasta = normalizePath(bin_files, mustWork = TRUE),
  fasta_name = basename(fasta),
  bin_id = str_remove(fasta_name, "\\.fa$"),
  unit = basename(dirname(dirname(fasta))),
  mode = basename(dirname(dirname(dirname(fasta))))) |>
  filter(mode %in% c("sample", "coassembly"),
         basename(dirname(fasta)) == paste0(unit, "_DASTool_bins"))

# Read CheckM2 and GUNC reports for each unit
units <- bins |> distinct(mode, unit)

checkm2 <- pmap_dfr(units, function(mode, unit) {
  read_tsv(file.path("data/bin_quality/checkm2", mode, unit, "quality_report.tsv"),
           show_col_types = FALSE) |>
    transmute(mode = mode, unit = unit, bin_id = Name,
              completeness = Completeness, contamination = Contamination)
})

gunc <- pmap_dfr(units, function(mode, unit) {
  read_tsv(file.path("data/bin_quality/gunc", mode, unit,
                     "GUNC.progenomes_3.maxCSS_level.tsv"),
           show_col_types = FALSE) |>
    transmute(mode = mode, unit = unit, bin_id = genome,
              css = clade_separation_score, rrs = reference_representation_score,
              pass_gunc = as.logical(pass.GUNC))
})

# Join quality results and select bins for dRep
summary <- bins |>
  left_join(checkm2, by = c("mode", "unit", "bin_id")) |>
  left_join(gunc, by = c("mode", "unit", "bin_id")) |>
  mutate(
    selection_status = case_when(
      is.na(completeness) | is.na(contamination) ~ "missing_checkm2",
      completeness < min_completeness ~ "low_completeness",
      contamination >= max_contamination ~ "high_contamination",
      is.na(pass_gunc) ~ "missing_gunc",
      !pass_gunc ~ "failed_gunc",
      TRUE ~ "selected"),
    selected_for_drep = selection_status == "selected") |>
  arrange(mode, unit, bin_id)

selected <- summary |> filter(selected_for_drep)

# Write the summary and dRep input files
write_tsv(summary, "tables/bin_quality_summary.tsv", na = "NA")
write_lines(selected$fasta, "data/dereplication/inputs/genomes.txt")

selected |>
  transmute(genome = fasta_name, completeness, contamination) |>
  write_csv("data/dereplication/inputs/genome_info.csv")

# Print the number of bins in each selection category
summary |> count(selection_status, name = "n_bins") |> print()
