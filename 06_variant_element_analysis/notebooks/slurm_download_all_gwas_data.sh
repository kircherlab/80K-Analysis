#!/bin/bash
# Slurm directives
#SBATCH --job-name=gwas_table_download  # Job array name
#SBATCH --output=slurm_logs/gwas_table_download_%A_%a.out
#SBATCH --error=slurm_logs/gwas_table_download_%A_%a.err
#SBATCH --array=0-4587        # 1-4 Adjust to number of files - 1
#SBATCH --cpus-per-task=1   # Number of CPUs/cores per task
#SBATCH --mem=1G            # Memory per task
#SBATCH --time=01:00:00     # Max Wall-clock time (HH:MM:SS) per task
#SBATCH --partition=short   # REPLACE with your cluster's partition (e.g., "shared", "compute")


# Configuration
MANIFEST_FILE="UKBB_GWAS_Imputed_v3_manifest_release_20180731_filtered.tsv"           # Your input manifest table file
OUTPUT_DIR="filtered_gwas_results_slurm_with_filtering"    # Directory to store results

# --- Main script logic for EACH array task ---
# Ensure directories exist (only create if they don't, idempotent)
mkdir -p "$OUTPUT_DIR"
mkdir -p "slurm_logs"

# Get the current array task ID (1-indexed)
TASK_ID=$SLURM_ARRAY_TASK_ID

# Echo a starting message to the .out file for tracking
echo "$(date): Slurm Array Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $TASK_ID - Starting processing."

# Extract the specific line from the manifest file for this task
MANIFEST_LINE=$(sed -n "$((TASK_ID + 1))p" "$MANIFEST_FILE")

# Parse the columns from the manifest line
IFS=$'\t' read -r PHENOTYPE_CODE PHENOTYPE_DESC UKB_LINK SEX FILE WGET_COMMAND AWS_FILE DROPBOX_FILE MD5S WGET_CMD_VALID <<< "$MANIFEST_LINE"

# Basic validation of the wget command extracted
if [[ -z "$WGET_COMMAND" ]]; then
    echo "$(date): Error for Task $TASK_ID: WGET_COMMAND is empty. Skipping this task." >&2
    exit 1 # Indicate failure for this specific task
fi

# --- NEW FILENAME GENERATION LOGIC ---
# Sanitize PHENOTYPE_CODE and PHENOTYPE_DESC for use in filenames
# Replace spaces with underscores, and remove characters not safe for filenames.
# This assumes common filename restrictions. Adapt if your cluster has unusual rules.
SANITIZED_PHENOTYPE_CODE=$(echo "$PHENOTYPE_CODE" | tr ' ' '_' | tr -cd '[:alnum:]_.-')
SANITIZED_PHENOTYPE_DESC=$(echo "$PHENOTYPE_DESC" | tr ' ' '_' | tr -cd '[:alnum:]_.-')

# Construct the output filename using both sanitized fields
OUTPUT_FILENAME_BASE="${SANITIZED_PHENOTYPE_CODE}_${SANITIZED_PHENOTYPE_DESC}_significant_gwas.tsv.gz"
OUTPUT_FILEPATH_TEMP="${OUTPUT_DIR}/${OUTPUT_FILENAME_BASE}" # Temporary name before potential rename

echo "$(date): Task $TASK_ID - Processing Phenotype Code: ${PHENOTYPE_CODE}, Description: ${PHENOTYPE_DESC} - Expected output: ${OUTPUT_FILEPATH_TEMP}"

# Execute the download and filter
if $WGET_COMMAND -qO- 2>/dev/null | zcat 2>/dev/null | \
    awk -F'\t' -v threshold=5e-8 '
    BEGIN {
        PVAL_COL = -1;
    }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            colname = tolower($i);
            if (colname == "pval" || colname ~ /pval/) {
                PVAL_COL = i;
                break;
            }
        }
        if (PVAL_COL == -1) {
            print "Error: P-value column not found in header." > "/dev/stderr";
            exit 1;
        }
        print;
    }
    NR > 1 {
        if (PVAL_COL != -1 && $PVAL_COL != "NaN" && $PVAL_COL <= threshold) {
            print;
        }
    }
    ' | gzip > "$OUTPUT_FILEPATH_TEMP";then

    # Check if the output file is empty (other than header)
    # dd reads the first 10k bytes, gzip -dc decompresses, head -n 2 gets the first two lines, wc -l counts them.
    # || echo 0 ensures LINE_COUNT is 0 if the file is truly empty or error occurs.
    LINE_COUNT=$(dd if="$OUTPUT_FILEPATH_TEMP" bs=1k count=10 2>/dev/null | gzip -dc 2>/dev/null | head -n 2 | wc -l || echo 0)

    if [[ "$LINE_COUNT" =~ ^[0-9]+$ && "$LINE_COUNT" -gt 1 ]]; then
        echo "$(date): Task $TASK_ID - Phenotype Code: ${PHENOTYPE_CODE} (${PHENOTYPE_DESC}) - Done. Filtered results saved to: ${OUTPUT_FILEPATH_TEMP}"
    else # If only 1 line (header) is found in the first 2 lines
        NEW_OUTPUT_FILENAME="empty_${OUTPUT_FILENAME_BASE}"
        NEW_OUTPUT_FILEPATH="${OUTPUT_DIR}/${NEW_OUTPUT_FILENAME}"
        mv "$OUTPUT_FILEPATH_TEMP" "$NEW_OUTPUT_FILEPATH"
        echo "$(date): Task $TASK_ID - Phenotype Code: ${PHENOTYPE_CODE} (${PHENOTYPE_DESC}) - No significant hits found below threshold. File contains only header. Renamed to: ${NEW_OUTPUT_FILEPATH}"
    fi
else
    echo "$(date): Error for Task $TASK_ID: Failed to download/process ${PHENOTYPE_CODE} (${PHENOTYPE_DESC}). The wget command was: ${WGET_COMMAND}" >&2
    rm -f "$OUTPUT_FILEPATH_TEMP" # Clean up partially downloaded/filtered file
    exit 1
fi

echo "$(date): Task $TASK_ID finished."
exit 0