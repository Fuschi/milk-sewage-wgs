#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_DIR="${PIPELINE_DIR}/workspace"

mkdir -p \
  "${WORKSPACE_DIR}/reads_raw" \
  "${WORKSPACE_DIR}/reads_trim" \
  "${WORKSPACE_DIR}/reads_clean" \
  "${WORKSPACE_DIR}/host" \
  "${WORKSPACE_DIR}/qc" \
  "${WORKSPACE_DIR}/logs"

snakemake \
  --snakefile "${PIPELINE_DIR}/Snakefile" \
  --directory "${WORKSPACE_DIR}" \
  --executor cluster-generic \
  --jobs 100 \
  --default-resources qos=normal mem_mb=5000 time="00:10:00" requeue=0 trigger=0 \
  --latency-wait 60 \
  --keep-going \
  --printshellcmds \
  --scheduler greedy \
  --use-conda \
  --local-cores 1 \
  --max-jobs-per-second 10 \
  --max-status-checks-per-second 1 \
  --cluster-generic-submit-cmd 'mkdir -p logs/slurm/{rule} && sbatch \
    --account=applicata \
    --nodes=1 \
    --qos={resources.qos} \
    --cpus-per-task={threads} \
    --mem={resources.mem_mb} \
    --job-name=smk-{rule}-{jobid} \
    --output=logs/slurm/{rule}/{rule}-{jobid}-%j.out \
    --error=logs/slurm/{rule}/{rule}-{jobid}-%j.err \
    $( [ "{resources.requeue}" -eq "1" ] && echo "--requeue" ) \
    $( [ "{resources.trigger}" -eq "1" ] && echo "--reservation=prj-trigger --nodelist=mtx30" ) \
    --time={resources.time}' \
  "$@"
