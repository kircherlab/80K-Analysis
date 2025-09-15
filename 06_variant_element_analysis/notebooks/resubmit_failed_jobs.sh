#!/bin/bash
# resubmit_failed_jobs.sh

RESTART_JOB_FILE="restart_jobs.txt"        # Output from extract_failed_tasks.sh
GWAS_PROCESSOR_SCRIPT="./submit_gwas_one_task.sh" # Your main worker script

if [ ! -f "$RESTART_JOB_FILE" ]; then
    echo "Error: Restart file '$RESTART_JOB_FILE' not found. Please run extract_failed_tasks.sh first." >&2
    exit 1
fi

echo "Resubmitting failed jobs from $RESTART_JOB_FILE..."

while IFS='=' read -r JOB_ID TASK_LIST; do
    echo "Resubmitting tasks for array job ID $JOB_ID: $TASK_LIST"
    # Submit with the --array option using the comma-separated list of task IDs
    # You can still add the %LIMIT here if you want to throttle the re-runs too
    sbatch --array="$TASK_LIST"%50 "$GWAS_PROCESSOR_SCRIPT"
    if [[ $? -ne 0 ]]; then
        echo "Warning: Failed to submit tasks for job ID $JOB_ID." >&2
    fi
done

echo "Finished attempting to resubmit all failed jobs."
echo "Use 'squeue -u \$USER' to monitor the newly submitted jobs."