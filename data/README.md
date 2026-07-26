# Data

## Project status

The bioinformatics pipeline has already been run on the complete dataset. Its large inputs and outputs—including FASTQ, FASTA, BAM files, reference databases, intermediate files, and complete tool directories—are stored in the local pipeline workspace and are intentionally excluded from Git.

The pipeline code is currently being reorganized and should therefore be considered a **work in progress**. However, the principal results required for the current analyses have already been generated and are available in this directory as lightweight tables.

The repository currently includes:

- MAG abundance and coverage estimates generated with **CoverM**;
- MAG taxonomic classifications generated with **GTDB-Tk**;
- MAG quality, completeness, and contamination statistics generated with **CheckM2**;
- reference-based taxonomic read counts generated with **Kaiju** and provided by **Fulvia Troja**;
- sample metadata and cleaned-read summary statistics;
- mappings between GTDB and NCBI taxonomies used to compare assembly-based and reference-based results.

The files under `data/` are sufficient for the downstream R and Quarto analyses currently included in the repository. The complete pipeline workspace is not required for these analyses and is not committed to Git.

## Compression

Tabular files are stored individually as gzip-compressed files, using extensions such as `.tsv.gz` and `.tab.gz`. Packages including `readr` can read these files directly; manual decompression is not required.

For example:

```r
library(readr)

abundances <- read_tsv(
  "data/organized/abundances_coverm.tsv.gz",
  show_col_types = FALSE
)
```

The Quarto source file `organize_data.qmd` and this README remain uncompressed so that they can be inspected directly on GitHub. TIFF figures are also retained in their original format because they are already internally compressed.

## Directory structure

The `raw/` directory contains lightweight outputs exported from the completed pipeline or supplied by collaborators. It does **not** contain raw sequencing files or the complete pipeline workspace.

The `organized/` directory contains harmonized, analysis-ready tables produced from the files under `raw/`. These are the preferred inputs for downstream analyses.

## Analysis-ready tables

### `organized/meta_samples.tsv.gz`

Sample metadata combined with sequencing statistics calculated from cleaned reads.

Main columns include:

- `sample_id`;
- `biome`;
- `num_paired_seqs`;
- `num_seqs`;
- `sum_len`;
- `avg_len`;
- `min_len` and `max_len`;
- additional sample-level descriptive variables when available.

Sources:

- `raw/samples/sample_metadata.tsv.gz`;
- `raw/samples/seqkit_cleaned_reads.tsv.gz`.

### `organized/gtbtk.tsv.gz`

Taxonomic classification of dereplicated MAGs generated with GTDB-Tk.

Main columns include:

- `taxa_id`;
- `classification`;
- `domain`;
- `phylum`;
- `class`;
- `order`;
- `family`;
- `genus`;
- `species`.

Sources:

- `raw/gtdbtk/gtdbtk.bac120.summary.tsv.gz`;
- `raw/gtdbtk/gtdbtk.ar53.summary.tsv.gz`.

### `organized/checkm2.tsv.gz`

Quality assessment of dereplicated MAGs generated with CheckM2.

The final `Completeness` field is selected according to the CheckM2 model used. When `Completeness_Model_Used` is `Gradient Boost (General Model)`, `Completeness_General` is used; otherwise, `Completeness_Specific` is used.

Main columns include:

- `sample_id`;
- `taxa_id`;
- `Completeness`;
- `Contamination`;
- `Strain_heterogeneity`;
- `GC_content`;
- `Genome_size`.

Sources:

- individual CheckM2 reports under `raw/checkm2/`.

### `organized/abundances_coverm.tsv.gz`

MAG abundance and coverage estimates generated with CoverM for each sample.

Main columns include:

- `sample_id`;
- `taxa_id`;
- `read_count`;
- `reads_per_base`;
- `relative_abundance`;
- `rpkm`;
- `tpm`;
- `mean`;
- `trimmed_mean`;
- `variance`;
- `covered_fraction`;
- `covered_bases`;
- `length`.

Sources:

- per-sample CoverM outputs under `raw/coverm/`.

### `organized/taxa_kaiju.tsv.gz`

Taxonomic metadata associated with the reference-based Kaiju classifications.

The Kaiju read-count tables under `raw/kaiju_sewage-milk_mags/` were provided by Fulvia Troja and include counts and relative abundances summarized at species, genus, and phylum levels.

### Taxonomy mapping tables

The following tables harmonize GTDB and NCBI taxonomic identifiers and names:

- `organized/map_gtdb_to_ncbi.tsv.gz`;
- `organized/map_ncbi_to_gtdb.tsv.gz`;
- `organized/map_gtdbtk_ncbi.tsv.gz`.

They are used when comparing MAG-based GTDB-Tk results with reference-based Kaiju classifications.

## GTDB metadata

The files under `raw/gtdb-metadata/` are the complete GTDB release 226 metadata tables used only to construct GTDB–NCBI mappings. They are **not analysis-ready project results** and do not need to be preserved in Git because they can be downloaded again from GTDB.

In particular, `bac120_metadata_r226.tsv.gz` remains larger than GitHub's normal per-file limit even after gzip compression. The entire directory should therefore be ignored by Git:

```gitignore
/data/raw/gtdb-metadata/
```

The metadata can be restored locally with the download commands contained in `organize_data.qmd`.

## Notes

- Only dereplicated MAGs are included in the assembly-based tables.
- Sample sequencing statistics are based on cleaned reads.
- Large sequencing files and complete pipeline outputs are not tracked by Git.
- `raw/gtdb-metadata/` is retained only in the complete local archive and omitted from the Git-ready archive.
- The tables under `organized/` are the preferred inputs for downstream analyses.
