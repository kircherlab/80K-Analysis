#!/bin/bash
# run_gwas_processing.sh

# SLURM Directives:
#SBATCH --job-name=gwas_process    # Job name
#SBATCH --output=logs/gwas_process_%A_%a.out # Standard output and error log
#SBATCH --error=logs/gwas_process_%A_%a.err  # %A is job ID, %a is array task ID
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --cpus-per-task=1          # Request 1 CPU core
#SBATCH --mem=4G                   # Request 4 GB of memory per task
#SBATCH --time=00:30:00            # Max runtime for each task (HH:MM:SS)
#SBATCH --array=0-4218%100         # Array of jobs. Adjust to your total number of files - 1.
                                   # %100 limits to 100 concurrent jobs. Adjust as needed.

# --- Configuration Variables ---
# Directory containing your gzipped GWAS files
GWAS_FILES_DIR="filtered_gwas_results_slurm"
GWAS_FILES_DIR="filtered_gwas_results_slurm_with_filtering_one_at_a_time"
GWAS_FILES_DIR="gwas_significant_data"

# Path to the variant information file
VARIANT_INFO_FILE="80k_neuro_variant_info_filtered_spdi_tested.tsv"

# Directory where individual processed output files will be saved
OUTPUT_DIR="processed_gwas_overlaps"

# Python script to run for each GWAS file
PYTHON_SCRIPT="process_gwas.py"

# --- Create necessary directories ---
mkdir -p "$OUTPUT_DIR"
mkdir -p "logs" # For SLURM logs

# Conda activate
source /data/cephfs-1/home/users/kisa11_c/work/miniforge3/etc/profile.d/conda.sh
conda activate mobil

# # --- Get the list of GWAS files ---
# # Using `find` to ensure robust handling of filenames and to filter by extension
# # Sorted to ensure consistent array indexing if files are added/removed.
# # `readarray` reads lines from stdin into an array.
readarray -t GWAS_FILES < <(find "$GWAS_FILES_DIR" -maxdepth 1 -name "*_significant_gwas.tsv.gz" | sort)

# GWAS_LIST_FILE="/data/cephfs-1/scratch/groups/kircher/MPRA/IGVF_Y1_design/projects/80K_MPRA/ludwig_variant_annotations/download_testing/workflow/valid_gwas_files_snakemake_download.txt"

# # Check that the list file exists
# if [ ! -f "$GWAS_LIST_FILE" ]; then
#     echo "GWAS list file not found: $GWAS_LIST_FILE"
#     exit 1
# fi

# # Read the list of filenames from the file and prepend the directory path
# readarray -t GWAS_FILES < <(sed "s|^|$GWAS_FILES_DIR/|" "$GWAS_LIST_FILE")

# Check if the array is empty
if [ ${#GWAS_FILES[@]} -eq 0 ]; then
    echo "No GWAS files listed in $GWAS_LIST_FILE. Exiting."
    exit 1
fi

# Get the specific GWAS file for this array task
FILE_TO_PROCESS="${GWAS_FILES[$SLURM_ARRAY_TASK_ID]}"

# Check if the file path is valid
if [ -z "$FILE_TO_PROCESS" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID is out of bounds for the GWAS_FILES array."
    exit 1
fi

echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing GWAS file: $FILE_TO_PROCESS"
echo "Variant Info File: $VARIANT_INFO_FILE"
echo "Output Directory: $OUTPUT_DIR"

# Run the Python script
python "$PYTHON_SCRIPT" "$FILE_TO_PROCESS" "$VARIANT_INFO_FILE" "$OUTPUT_DIR"

# Deactivate virtual environment if you activated one
# deactivate # Uncomment if using venv
