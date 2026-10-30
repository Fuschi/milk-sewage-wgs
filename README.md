# milk-sewage WGS

Snakemake workflow for shotgun metagenomic analysis of milk and sewage samples
on the OPH Slurm cluster. It covers read processing, assembly, genome binning,
quality assessment, dereplication, taxonomy, and resistance and virulence
profiling from reads and MAGs.

## Configuration

- Reads: `data/reads_raw/{sample_id}_R1.fastq.gz` and `{sample_id}_R2.fastq.gz`.
- Samples: `config/samples.tsv`.
- References and containers: `config/config.yaml`.
- Cluster profile: `profiles/slurm/config.yaml`.

See [environment setup](docs/conda.md) for installation instructions.

## Run

```bash
conda activate snakemake-slurm
snakemake --profile profiles/slurm --dry-run
./submit.sh
```

## Outputs

- `data/`: processed reads and intermediate outputs; ignored by Git.
- `results/`: analysis tables, metadata, and Sourmash results. See
  [results documentation](results/README.md).
- `logs/`: execution logs; ignored by Git.
- `benchmarks/`: runtime and resource usage.

Existing OPH runs require path migration before using this layout.
