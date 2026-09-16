# Snakemake environment

Snakemake runs from a Conda environment named `snakemake-slurm`. This environment
contains the workflow controller; bioinformatics tools run in Apptainer containers.

Create it once:

```bash
conda create -n snakemake-slurm --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda \
  "python>=3.11" "snakemake>=9" "snakemake-executor-plugin-slurm>=1" pandas
```

Activate it before using Snakemake:

```bash
conda activate snakemake-slurm
snakemake --cores 1 --dry-run
```

Run the dry-run from the repository root. `submit.sh` activates the same
environment automatically, using `~/miniforge3/etc/profile.d/conda.sh`.
Adjust that path if Conda is installed elsewhere. Apptainer must be available
on the cluster execution nodes.
