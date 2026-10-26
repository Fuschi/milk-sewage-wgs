library(tidyverse)

samples_file <- snakemake@input[["samples"]]
samples <- read_tsv(
  samples_file,
  col_types = cols(.default = col_character())
)
root <- dirname(dirname(normalizePath(samples_file)))

# Find FASTA files recursively.
bin_files <- function(path) {
  list.files(
    path,
    pattern = "\\.(fa|fna|fasta)(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

# Read FASTA files in chunks and calculate assembly statistics.
fasta_stats <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(con))

  n_contigs <- 0L
  total_bp <- sequence_bp <- 0
  min_bp <- Inf

  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (!length(lines)) break

    for (line in lines) {
      if (startsWith(line, ">")) {
        if (sequence_bp > 0) min_bp <- min(min_bp, sequence_bp)
        n_contigs <- n_contigs + 1L
        sequence_bp <- 0
      } else {
        line_bp <- nchar(gsub("[[:space:]]", "", line))
        sequence_bp <- sequence_bp + line_bp
        total_bp <- total_bp + line_bp
      }
    }
  }

  if (sequence_bp > 0) min_bp <- min(min_bp, sequence_bp)

  tibble(
    contigs = n_contigs,
    total_bp = total_bp,
    min_contig_bp = min_bp
  )
}

# Define bin directories for sample and multisample runs.
sample_runs <- samples |>
  transmute(sample = sample_id, biome) |>
  crossing(tool = c("comebin", "semibin2", "vamb")) |>
  mutate(
    mode = if_else(tool == "comebin", "sample", "multisample"),
    directory = case_when(
      tool == "comebin" ~ file.path(
        root, "data/binning/comebin/sample",
        sample, "comebin_res/comebin_res_bins"
      ),
      tool == "semibin2" ~ file.path(
        root, "data/binning/semibin2/multisample",
        biome, "samples", sample, "output_bins"
      ),
      tool == "vamb" ~ file.path(
        root, "data/binning/vamb/multisample",
        biome, "run/bins"
      )
    )
  )

# Define bin directories for coassembly runs.
coassembly_runs <- samples |>
  distinct(biome) |>
  crossing(tool = c("comebin", "semibin2")) |>
  mutate(
    sample = NA_character_,
    mode = "coassembly",
    directory = case_when(
      tool == "comebin" ~ file.path(
        root, "data/binning/comebin/coassembly",
        biome, "comebin_res/comebin_res_bins"
      ),
      tool == "semibin2" ~ file.path(
        root, "data/binning/semibin2/coassembly",
        biome, "output_bins"
      )
    )
  )

# List FASTA files and assign VAMB bins to their source sample.
inventory <- bind_rows(sample_runs, coassembly_runs) |>
  mutate(fasta = map(directory, bin_files)) |>
  unnest_longer(fasta) |>
  filter(
    tool != "vamb" |
      startsWith(basename(fasta), paste0(sample, "_k"))
  ) |>
  mutate(
    bin_id = if_else(
      mode == "coassembly",
      str_c(tool, mode, biome, basename(fasta), sep = "__"),
      str_c(tool, mode, biome, sample, basename(fasta), sep = "__")
    ),
    fasta = map_chr(fasta, normalizePath)
  ) |>
  select(bin_id, tool, mode, biome, sample, fasta)

# Add FASTA statistics and write the inventory.
inventory <- inventory |>
  mutate(stats = map(fasta, fasta_stats)) |>
  unnest(stats)

write_tsv(inventory, snakemake@output[["inventory"]], na = "")
