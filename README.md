# Milk–sewage WGS

Snakemake workflow for paired-end metagenomic read cleaning on the OPH cluster.
The workflow currently stops after KneadData; assembly and binning will be added
in subsequent steps.

## Setup

- Place reads in `data/reads_raw/{sample_id}_R1.fastq.gz` and
  `data/reads_raw/{sample_id}_R2.fastq.gz`.
- List samples in `config/samples.tsv`; `sample_id` must be non-empty and unique.
  IDs may contain letters, digits, underscores and hyphens. `biome` is retained
  as metadata but does not affect cleaning.
- `--sequencer-source NexteraPE` is set in `rules/cleaning.smk`, based on the Nextera
  adapter signal in raw FastQC reports for both mates of samples 3550 and 3554.
  This identifies a compatible adapter family, not the exact commercial kit.
  The setting applies to all samples: check the remaining batches before a full
  run. Mixed or custom kits require adapting the rule.
- Verify the Bos taurus index prefix. The rule expects the six `.bt2` index files.

Create the controller environment following [these instructions](docs/conda.md).

## Cleaning

Each sample is processed with KneadData 0.12.4:

1. FastQC on raw reads.
2. Standard KneadData/Trimmomatic adapter and quality trimming for the selected
   kit: initial minimum length 60, adapter clipping, `SLIDINGWINDOW:4:20`,
   and a final minimum length of 50% of the first input read length.
3. Bowtie2 host removal against the configured Bos taurus reference, using
   `--very-sensitive-local`. Strict pair filtering removes both mates when
   either aligns to the host. Surviving trimming singletons are also filtered.
4. FastQC on final reads and a per-sample read-count table.

Tandem-repeat filtering is explicitly disabled (`--bypass-trf`). Intermediate
files are removed by KneadData. Trimming options are not overridden, because a
custom `--trimmomatic-options` string replaces the adapter-trimming defaults too.
Resource limits and cleaning options are defined directly in `rules/cleaning.smk`.
`config/config.yaml` contains only the sample sheet, reference and container paths.

## Outputs

Native outputs are retained under `data/kneaddata/{sample_id}/`:

- `{sample_id}_kneaddata_paired_1.fastq` and `_paired_2.fastq`: cleaned pairs;
- `{sample_id}_kneaddata_unmatched_1.fastq` and `_unmatched_2.fastq`: cleaned
  singletons, when produced; they remain separate and are not discarded;
- `{sample_id}_kneaddata.log`: native log with parameters and read counts;
- `fastqc/`: raw and final quality reports.

FASTQ files are uncompressed. Empty/missing read categories can occur when no
reads survive the corresponding stage; the native run directory is tracked as
a whole. Do not manually modify its contents, as Snakemake tracks directory
completion rather than every contained file. A rerun replaces that directory.

Read-count tables are written to `tables/read_stats/kneaddata/{sample_id}.tsv`.
Counts refer to individual reads; the sample label is `{sample_id}_kneaddata`.
Console logs and resource measurements are in `logs/kneaddata/` and
`benchmarks/kneaddata/`.

Old `data/reads_trim/`, `data/reads_clean/` and assembly outputs are not consumed
or deleted by this workflow. The new output namespace prevents their reuse as
KneadData results. The [singleton correction note](docs/singleton-host-filtering.md)
refers to the retired BBDuk/Bowtie2 workflow.

## Run

Activate the environment and inspect the workflow:

```bash
conda activate snakemake-slurm
snakemake --cores 1 --dry-run
```

On the cluster, first run the two samples whose raw FastQC reports were checked:

```bash
snakemake --profile profiles/slurm data/kneaddata/3550 data/kneaddata/3554
```

Inspect their final FastQC reports, read-count tables and logs, and check adapter
compatibility in the other batches before processing all samples:

```bash
./submit.sh
```

`submit.sh` uses `~/miniforge3` and the `snakemake-slurm` environment.
SLURM execution settings are in `profiles/slurm/config.yaml`.
