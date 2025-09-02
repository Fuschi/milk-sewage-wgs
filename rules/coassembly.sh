#!/bin/bash
#SBATCH --account=applicata
#SBATCH --reservation=prj-trigger
#SBATCH --qos=long
#SBATCH --partition=m8
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --requeue
#SBATCH --job-name=coassembly
#SBATCH --output=logs/slurm/coassembly_megahit_genome/%x_%j.out
#SBATCH --error=logs/slurm/coassembly_megahit_genome/%x_%j.err
 
set -euo pipefail
 
# Verifica se la variabile BIOME è stata passata
if [ -z "${BIOME:-}" ]; then
    echo "ERRORE: devi passare BIOME=Milk o BIOME=Sewage"
    exit 1
fi
 
echo "=== Co-assembly per $BIOME ==="
 
# Attiva Conda
source /home/PERSONALE/niccolo.barbieri3/miniforge3/etc/profile.d/conda.sh
conda activate megahit
 
# Input files
R1="snakestream/coassembly_megahit_genome/conc_reads/${BIOME}/all_R1.fastq.gz"
R2="snakestream/coassembly_megahit_genome/conc_reads/${BIOME}/all_R2.fastq.gz"
SINGLE="snakestream/coassembly_megahit_genome/conc_reads/${BIOME}/all_sing.fastq.gz"
 
# Output directory
OUTDIR="snakestream/coassembly_megahit_genome/${BIOME}"
#mkdir -p "$OUTDIR"
 
# Log file
LOG="logs/coassembly_megahit_genome/${BIOME}.log"
 
# Lancio Megahit
srun megahit \
    -t 64 \
    --memory 100000 \
    --verbose \
    -1 "$R1" \
    -2 "$R2" \
    -r "$SINGLE" \
    --continue \
    -o "$OUTDIR" \
> "$LOG" 2>&1
mv "${OUTDIR}"/final.contigs.fa "${OUTDIR}"/coassembly.final.contigs.fa
echo "OK -> ${OUTDIR}/coassembly/final.contigs.fa"
