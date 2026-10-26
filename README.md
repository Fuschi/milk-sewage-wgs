# milk-sewage WGS

Snakemake workflow for shotgun metagenomic analysis of milk and sewage samples,
configured for the OPH SLURM cluster.

The workflow covers read quality control and host removal, metagenome assembly,
genome binning, refinement and quality assessment, and read-based profiling of
antimicrobial resistance and virulence genes.

## Configuration

- Input reads: `data/reads_raw/{sample_id}_R1.fastq.gz` and
  `data/reads_raw/{sample_id}_R2.fastq.gz`.
- Sample metadata: `config/samples.tsv`.
- References and containers: `config/config.yaml`.
- Cluster settings: `profiles/slurm/config.yaml`.

See [environment setup](docs/conda.md) for installation instructions.

## Run

Activate the environment and preview the workflow:

```bash
conda activate snakemake-slurm
snakemake --profile profiles/slurm --dry-run
```

Submit the workflow:

```bash
./submit.sh
```

## Outputs

- `data/`: processed reads and analysis outputs (ignored by Git).
- `tables/`: summary tables.
- `results/`: plots and other analysis results.
- `logs/`: execution logs (ignored by Git).
- `benchmarks/`: runtime and resource usage.
