# Results

Tables are tab-separated unless stated otherwise. `NA` and empty fields denote
missing values.
Sample identifiers correspond to `sample_id` in `config/samples.tsv`. MAG
identifiers refer to dereplicated genome filenames without the `.fa` extension,
except in `bin_quality_summary.tsv`, which describes bins before dereplication.

## Read processing and assembly

### `kneaddata_read_counts.tsv`

One row per sample. Counts of paired and orphan reads at the raw, trimmed,
host-decontamination and final stages. `Sample` identifies the sample.

### `metaquast_summary.tsv`

One row per individual assembly or biome coassembly. `assembly_type`,
`sample_id`, `biome` and `report_file` identify the assembly and source
report; a coassembly has no individual `sample_id`. Remaining columns contain
MetaQUAST statistics, including contig counts, total length, GC, N50, N90, auN,
L50 and L90.

## Bin quality and MAG taxonomy

### `bin_quality_summary.tsv`

One row per DAS Tool bin, before dereplication. `bin_id`, `unit` and `mode`
identify the bin and its sample or coassembly; `fasta` and `fasta_name`
identify its sequence file. `completeness` and `contamination` are CheckM2
percentages; `css`, `rrs` and `pass_gunc` come from GUNC.
`selection_status` gives the selection outcome and `selected_for_drep` marks
bins with completeness ≥50%, contamination <10% and a GUNC pass.

### `gtdbtk_summary.tsv`

One row per classified dereplicated MAG. `mag_id` identifies the MAG;
`domain` through `species` split the GTDB classification into ranks, while
`classification` retains the original lineage. The remaining GTDB-Tk fields
include reference matches, ANI, alignment fractions, classification method,
marker information and warnings. Unassigned ranks are `NA`.

### `ncbi_taxonomy.tsv`

One row per MAG with GTDB and inferred NCBI classifications. `ncbi_domain`
through `ncbi_species` contain the separate NCBI ranks; `ncbi_most_specific`
is the most specific assigned rank. The NCBI lineage is inferred by majority
vote and may be incomplete. Join to other MAG tables on `mag_id`.

## MAG abundance and gene screening

### `mags_abundance.tsv`

CoverM estimates in long format, one row per `sample` and `mag_id`:

- `relative_abundance`: percentage relative abundance.
- `mean`, `trimmed_mean`: mean and trimmed mean coverage depth.
- `covered_fraction`, `covered_bases`: fraction and number of covered bases.
- `variance`, `rpkm`, `tpm`: coverage variance and normalized estimates.

The `unmapped` row reports the percentage of reads not mapped to the MAG
collection; its other metrics are `NA`. CoverM uses its default minimum
covered fraction of 10%; values below this threshold can be reported as zero,
so a zero does not establish complete absence. Only paired cleaned reads are
used for this profiling.

### `abricate_resfinder_2.6.0.tsv`

ABRicate hits against ResFinder in dereplicated MAGs. Each row is a reported
gene hit: `#FILE` is the MAG FASTA filename, `SEQUENCE` is the contig,
`START`, `END` and `STRAND` locate the hit, and `GENE` identifies the
reference. `%IDENTITY` and `%COVERAGE` summarize the alignment; other
columns retain the ABRicate output. Hits are filtered at ≥90% identity and
≥60% reference coverage.

### `abricate_vfdb_2026-10-26.tsv`

ABRicate hits against the VFDB nucleotide Set A in dereplicated MAGs, with the
same columns and thresholds as the ResFinder table. `GENE` contains the VFG
identifier and `PRODUCT` retains the reference description.

The custom ABRicate database headers use `.` as a placeholder for fields
such as `ACCESSION` and `RESISTANCE`; consult the metadata tables below for
annotations. An individual hit can cover only part of its reference gene.

## Read-based resistance and virulence profiling

### `resfinder_reads_abundance.tsv`

KMA-based ResFinder quantification in long format, one row per `sample_id` and
reported `gene_id`. Results retain references passing ≥90% `query_identity_pct`
and ≥60% `template_coverage_pct`. The columns include:

- `length_bp`, `template_identity_pct`, `query_identity_pct`,
  `template_coverage_pct`, `depth`: reference and alignment statistics.
- `aligned_fragments`, `aligned_bases`, `total_fragments`: mapping counts.
- `fpm`: aligned fragments per million input fragments; `fpkm` additionally
  accounts for reference length in kilobases.
- `tpm_filtered`: length-normalized fragments per million, calculated over
  references passing the filters within each sample.
- `fragment_proportion_pct`: aligned fragments as a percentage of input
  fragments.
- `pass_filters`: filter status; rows in this table have passed the filters.

### `vfdb_reads_abundance.tsv`

KMA-based VFDB Set A quantification with the same sample-level statistics and
filters as the ResFinder table. `gene_id` is the VFG identifier and
`reference` retains the VFDB reference header.

### `resfinder_metadata.tsv`

ResFinder annotations indexed by `gene_id`: antibiotic class, phenotype,
publication identifier, resistance mechanism, notes and required gene. Join
with `resfinder_reads_abundance.tsv` on `gene_id`, or with the ABRicate
ResFinder table using `GENE`.

### `vfdb_metadata.tsv`

VFDB annotations indexed by VFG `gene_id`: accession, gene name, organism,
virulence factor and category identifiers, names, characteristics, function,
mechanism and literature. Join with `vfdb_reads_abundance.tsv` on `gene_id`,
or with the ABRicate VFDB table using `GENE`.

## Sourmash

### `sourmash/{all,milk,sewage}_{angular,jaccard}_distance.csv`

Six sample-by-sample distance matrices: one for each group and metric. The
header lists sample IDs; rows follow the same order and have no separate row
labels. Angular distance uses abundance-aware sketches; Jaccard ignores
abundance. Corresponding matrices, labels and plots are in
`sourmash/{group}/{metric}/`.

Within each group and metric directory, `labels.csv` records sketch labels,
order and source paths; `distance.mat.labels.txt` lists sample labels.
Source paths in existing label files reflect the run that produced them and
may predate the workflow's path migration.

## Linking tables

- Join `mags_abundance.tsv` to `gtdbtk_summary.tsv` on `mag_id`; exclude the
  special `unmapped` row when a MAG identifier is required.
- To link ABRicate hits to a MAG table, remove the `.fa` suffix from `#FILE`
  and match it to `mag_id`. A MAG can have multiple hits.
- Read-based abundance rows and ABRicate hits use different units: a
  sample–gene quantification and a gene hit on a MAG contig, respectively.
