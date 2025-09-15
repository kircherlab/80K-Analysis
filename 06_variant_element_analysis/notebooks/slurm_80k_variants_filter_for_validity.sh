#!/bin/bash
#SBATCH --job-name=check_gz
#SBATCH --output=logs/check_gz_%A_%a.out
#SBATCH --error=logs/check_gz_%A_%a.err
#SBATCH --array=0-4323        # Adjust to number of files - 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00       # Adjust as needed
#SBATCH --partition=short     # Use a suitable partition
#SBATCH --mail-type=BEGIN,END,FAIL         # optional: email on end/failure
#SBATCH --mail-user=kilian.salomon@bih-charite.de  # optional: your email

# find filtered_gwas_results_slurm/ -name "*.tsv.gz" > all_files_to_check_250625.txt

# find filtered_gwas_results_slurm -name "*.tsv.gz" > all_files_to_check.txt
source /data/cephfs-1/home/users/kisa11_c/work/miniforge3/etc/profile.d/conda.sh
conda activate mobil
# Directories
VALID_DIR="gwas_significant_data"
INVALID_DIR="gwas_header_only"
VALID_LIST="valid_gwas_files250625.txt"

mkdir -p "$VALID_DIR" "$INVALID_DIR" logs

# Get the current file from the list
FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" all_files_to_check_250625.txt)
BASENAME=$(basename "$FILE")
echo "Processing file: $BASENAME"
# Skip if already processed
if [[ -e "$VALID_DIR/$BASENAME" || -e "$INVALID_DIR/$BASENAME" ]]; then
    echo "[SKIP] Already processed: $BASENAME"
    exit 0
fi

# Check gzip integrity
if gzip -t "$FILE" 2>/dev/null; then
    # LINE_COUNT=$(dd if="$FILE" bs=1k count=10 2>/dev/null | gzip -dc | head -n 2 | wc -l)
    LINE_COUNT=$(dd if="$FILE" bs=1k count=10 2>/dev/null | gzip -dc 2>/dev/null | head -n 2 | wc -l || echo 0)

    if [[ "$LINE_COUNT" =~ ^[0-9]+$ && "$LINE_COUNT" -gt 1 ]]; then
        ln -s "$(realpath "$FILE")" "$VALID_DIR/$BASENAME"
        echo "$VALID_DIR/$BASENAME" >> "$VALID_LIST"
        echo "[OK] Valid and complete: $BASENAME"
    else
        mv "$FILE" "$INVALID_DIR/"
        echo "[SKIP] Header-only or empty: $BASENAME -> moved to $INVALID_DIR"
    fi
else
    mv "$FILE" "$INVALID_DIR/"
    echo "[FAIL] Corrupted gzip: $BASENAME -> moved to $INVALID_DIR"
fi
