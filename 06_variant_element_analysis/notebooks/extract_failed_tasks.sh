#!/bin/bash
# extract_failed_tasks.sh

FAILED_LOG_FILE="failed_16799067_jobs.txt"
RESTART_JOB_FILE="restart_jobs.txt" # Output file for sbatch --array ranges

# Clear the output file from previous runs
> "$RESTART_JOB_FILE"

# Use awk to parse the log file paths and group task IDs by array Job ID
# This command will output lines like:
# <JOB_ID>=<TASK_ID>,<TASK_ID>,<TASK_ID>...
awk -F'[_.]' '{
    job_id = $4; # e.g., 16799067
    task_id = $5; # e.g., 2736
    if (job_id && task_id) {
        if (!(job_id in tasks)) {
            tasks[job_id] = "";
        }
        if (tasks[job_id] != "") {
            tasks[job_id] = tasks[job_id] "," task_id;
        } else {
            tasks[job_id] = task_id;
        }
    }
} END {
    for (job_id in tasks) {
        print job_id "=" tasks[job_id];
    }
}' "$FAILED_LOG_FILE" | while IFS='=' read -r JOB_ID TASK_LIST; do
    # For each Job ID, construct the --array string and append to the restart file
    echo "$JOB_ID=$TASK_LIST" >> "$RESTART_JOB_FILE"
done

echo "Extracted failed tasks and prepared restart file: $RESTART_JOB_FILE"
echo "Example content of $RESTART_JOB_FILE:"
head "$RESTART_JOB_FILE" # Show first few lines