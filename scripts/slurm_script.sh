#!/usr/bin/env bash
set -eu

JOB_NAME="snakemake_master_process"
LOG_DIR="logs"

if [[ ! -d "$LOG_DIR" ]]; then
    echo "Error: Log directory $LOG_DIR does not exist"
    exit 1
fi

MEMORY="32G"
TIME="${TIME:-3d}"
THREADS=2
BINDS="/scratch/user/s4897040"
SINGULARITY_ARGS="-B $BINDS"
DEFAULT_TMP="slurm_account=a_uqccr"
CMD="snakemake --logger snakesee --sdm conda apptainer --executor slurm --jobs 500 --default-resources $DEFAULT_TMP --slurm-init-seconds-before-status-checks=20 --rerun-incomplete --local-cores $THREADS $* --singularity-args '$SINGULARITY_ARGS'"

ssubmit -t "$TIME" -m "$MEMORY" -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e "$JOB_NAME" "$CMD" -- -c "$THREADS" 
