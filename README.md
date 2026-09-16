# Milk–sewage WGS

Snakemake workflow for metagenomic data processing on the OPH cluster.

## Structure

```text
config/       Sample sheet and workflow configuration.
rules/        Snakemake rules.
scripts/      Helper scripts.
docs/         Setup instructions.
profiles/     Cluster execution settings.
data/         Workflow data (not tracked by Git).
logs/         Execution logs.
benchmarks/   Resource usage records.
```

## Setup

- Place reads in `data/reads_raw/{sample_id}_R1.fastq.gz` and
  `data/reads_raw/{sample_id}_R2.fastq.gz`.
- List samples in the `sample_id` column of `config/samples.tsv`.
- Check container and reference paths in `config/config.yaml`.

Create the Conda environment following [these instructions](docs/conda.md).

## Run

From the repository root, activate the environment and inspect the workflow
without processing reads:

```bash
conda activate snakemake-slurm
snakemake --cores 1 --dry-run
```

To launch the pipeline on the OPH cluster:

```bash
./submit.sh
```

`submit.sh` uses `~/miniforge3` and the `snakemake-slurm` environment.
SLURM and Apptainer settings are in `profiles/slurm/config.yaml`.
