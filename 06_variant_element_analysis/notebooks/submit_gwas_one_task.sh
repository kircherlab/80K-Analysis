#!/bin/bash
# gwas_array_processor.sh - This script is EXECUTED by each array task

# Slurm directives (these are fixed for all tasks in the array, set here)
#SBATCH --job-name=gwas_array_task   # A placeholder name for individual tasks
#SBATCH --output=slurm_logs/gwas_array_%A_%a.out  # Standard output file
#SBATCH --error=slurm_logs/gwas_array_%A_%a.err   # Standard error file
#SBATCH --partition=short            # REPLACE with your cluster's partition
#SBATCH --time=01:00:00              # Max Wall-clock time (HH:MM:SS) per task
#SBATCH --mem=1G                     # Memory per task
#SBATCH --cpus-per-task=1            # Number of CPUs/cores per task

# Configuration
MANIFEST_FILE="UKBB_GWAS_Imputed_v3_manifest_release_20180731_filtered.tsv"
OUTPUT_DIR="filtered_gwas_results_slurm"
PVAL_COLUMN_IN_GWAS_FILE=11

# --- Main script logic for EACH array task ---
# Ensure directories exist (only create if they don't, idempotent)
# These mkdirs will run for every task, which is harmless due to -p
mkdir -p "$OUTPUT_DIR"
mkdir -p "slurm_logs"

# Get the current array task ID (1-indexed)
TASK_ID=$SLURM_ARRAY_TASK_ID

# Echo a starting message to the .out file for tracking
# These variables will now be correctly populated by Slurm for each task
echo "$(date): Slurm Array Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $TASK_ID - Starting processing."

# Extract the specific line from the manifest file for this task
# This will correctly fetch the Nth line of the manifest based on TASK_ID
MANIFEST_LINE=$(sed -n "$((TASK_ID + 1))p" "$MANIFEST_FILE")

# Parse the columns from the manifest line
IFS=$'\t' read -r PHENOTYPE_CODE PHENOTYPE_DESC UKB_LINK SEX FILE WGET_COMMAND AWS_FILE DROPBOX_FILE MD5S WGET_CMD_VALID <<< "$MANIFEST_LINE"

# Basic validation of the wget command extracted
if [[ -z "$WGET_COMMAND" ]]; then
    echo "$(date): Error for Task $TASK_ID: WGET_COMMAND is empty. Skipping this task." >&2
    exit 1 # Indicate failure for this specific task
fi

# Sanitize PHENOTYPE_CODE and PHENOTYPE_DESC for use in filenames
SANITIZED_PHENOTYPE_CODE=$(echo "$PHENOTYPE_CODE" | tr ' ' '_' | tr -cd '[:alnum:]_.-')
SANITIZED_PHENETYPE_DESC=$(echo "$PHENOTYPE_DESC" | tr ' ' '_' | tr -cd '[:alnum:]_.-')

OUTPUT_FILENAME="${SANITIZED_PHENOTYPE_CODE}_${SANITIZED_PHENOTYPE_DESC}_significant_gwas.tsv.gz"
OUTPUT_FILEPATH="${OUTPUT_DIR}/${OUTPUT_FILENAME}"

echo "$(date): Task $TASK_ID - Processing Phenotype Code: ${PHENOTYPE_CODE}, Description: ${PHENOTYPE_DESC} - Expected output: ${OUTPUT_FILEPATH}"

# Execute the download and filter
if ! eval "${WGET_COMMAND} -qO- 2>/dev/null | zcat 2>/dev/null | awk -F'\t' 'NR==1 || (\$${PVAL_COLUMN_IN_GWAS_FILE} != \"NaN\" && \$${PVAL_COLUMN_IN_GWAS_FILE} <= 5e-8) {print}' | gzip > \"${OUTPUT_FILEPATH}\""; then
    echo "$(date): Error for Task $TASK_ID: Failed to download/process ${PHENOTYPE_CODE} (${PHENOTYPE_DESC}). The wget command was: ${WGET_COMMAND}" >&2
    rm -f "$OUTPUT_FILEPATH" # Clean up partially downloaded/filtered file
    exit 1 # Indicate failure for this specific task
else
    # Check if the output file is empty (other than header)
    # Note: wc -l on a gzipped file needs zcat
    NUM_LINES_OUT=$(zcat "$OUTPUT_FILEPATH" | wc -l) # Use zcat for accurate line count on gzipped file
    if [[ "$NUM_LINES_OUT" -eq 1 ]]; then
        echo "$(date): Task $TASK_ID - Phenotype Code: ${PHENOTYPE_CODE} (${PHENOTYPE_DESC}) - No significant hits found below threshold. File contains only header."
    else
        echo "$(date): Task $TASK_ID - Phenotype Code: ${PHENOTYPE_CODE} (${PHENOTYPE_DESC}) - Done. Filtered results saved to: ${OUTPUT_FILEPATH}"
    fi
fi

echo "$(date): Task $TASK_ID finished."
exit 0