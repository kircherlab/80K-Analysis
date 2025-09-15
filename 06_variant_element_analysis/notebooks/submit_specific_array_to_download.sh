#!/bin/bash
# submit_gwas_array.sh - This script SUBMITS the array job

# Configuration
MANIFEST_FILE="UKBB_GWAS_Imputed_v3_manifest_release_20180731_filtered.tsv"
ARRAY_PROCESSOR_SCRIPT="./submit_gwas_one_task.sh" # Path to the worker script

echo "Preparing to submit GWAS download and filter array job..."
echo "Manifest file: ${MANIFEST_FILE}"
echo "Array processor script: ${ARRAY_PROCESSOR_SCRIPT}"

# Determine the total number of lines in the manifest (minus header)
NUM_TASKS=$(wc -l < "$MANIFEST_FILE") # Use < for efficiency
NUM_TASKS=$((NUM_TASKS - 1))          # Subtract 1 for the header line

if [[ "$NUM_TASKS" -le 0 ]]; then
    echo "Error: Manifest file is empty or contains only a header. No tasks to submit." >&2
    exit 1
fi

echo "Total number of tasks to submit: ${NUM_TASKS}"

# Submit the array job using sbatch.
# The --array directive is now correctly determined BEFORE sbatch is called.
# You can add the %LIMIT here if you want to throttle, e.g., %50
sbatch --array=1-${NUM_TASKS}%50 "$ARRAY_PROCESSOR_SCRIPT"

if [[ $? -eq 0 ]]; then
    echo "Slurm job array submitted successfully."
    echo "You can monitor jobs with 'squeue -u \$USER' or 'scontrol show job <jobid>'."
    echo "Check 'slurm_logs' directory for job outputs and errors."
else
    echo "Error: Slurm job array submission failed." >&2
    exit 1
fi