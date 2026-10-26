library(tidyverse)

metadata <- read_tsv(
  snakemake@input[["metadata"]],
  col_types = cols(.default = col_character())
)

reports <- bind_rows(
  tibble(
    report_file = snakemake@input[["sample_reports"]],
    assembly_type = "sample"
  ),
  tibble(
    report_file = snakemake@input[["coassembly_reports"]],
    assembly_type = "coassembly"
  )
) |>
  mutate(
    id = basename(dirname(report_file)),
    sample_id = if_else(assembly_type == "sample", id, NA_character_),
    biome = if_else(
      assembly_type == "sample",
      metadata$biome[match(id, metadata$sample_id)],
      id
    )
  ) |>
  select(assembly_type, sample_id, biome, report_file)

summary <- reports |>
  mutate(
    metrics = map(report_file, \(path) {
      read_tsv(
        path,
        col_names = c("metric", "value"),
        col_types = cols(.default = col_character())
      ) |>
        pivot_wider(names_from = metric, values_from = value)
    })
  ) |>
  unnest(metrics) |>
  type_convert(col_types = cols(
    assembly_type = col_character(),
    sample_id = col_character(),
    biome = col_character(),
    report_file = col_character(),
    Assembly = col_character()
  ))

write_tsv(summary, snakemake@output[["summary"]], na = "")
