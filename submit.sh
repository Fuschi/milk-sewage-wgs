#!/usr/bin/env bash
set -euo pipefail

# Run this launcher directly: ./submit.sh (not sbatch submit.sh).
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

# SLURM opens controller logs before the job starts.
mkdir -p logs/controller

# Configure controller resources.
sbatch \
    --job-name=snakemake \
    --account=applicata \
    --qos=normal \
    --nodes=1 \
    --mem=12G \
    --time=23:59:00 \
    --cpus-per-task=2 \
    --chdir="$PROJECT_DIR" \
    --output="$PROJECT_DIR/logs/controller/%x_%j.out" \
    --error="$PROJECT_DIR/logs/controller/%x_%j.err" <<'JOB'
#!/usr/bin/env bash
set -euo pipefail

# Adapt this initialization path to your Conda installation.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-slurm

snakemake --profile profiles/slurm
JOB
