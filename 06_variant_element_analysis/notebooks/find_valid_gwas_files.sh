#!/bin/bash
# find_valid_gwas_files.sh

INPUT_DIR="gwas_filtered"
VALID_FILES_DIR="gwas_significant_data" # Directory for files with actual data
INVALID_FILES_DIR="gwas_header_only"    # Directory for files with only headers (or corrupted/incomplete)
VALID_FILE_LIST="valid_gwas_files_snakemake_download.txt"  # List of paths to valid files

mkdir -p "$VALID_FILES_DIR" "$INVALID_FILES_DIR"

# Clear previous list
> "$VALID_FILE_LIST"

echo "Checking files in ${INPUT_DIR} and separating them..."

for gz_file in "$INPUT_DIR"/*.tsv.gz; do
    if [ -f "$gz_file" ]; then # Ensure it's a regular file
        # Step 1: Check gzip integrity
        if gzip -t "$gz_file" 2>/dev/null; then
            # gzip integrity check passed, now check line count
            LINE_COUNT=$(zcat "$gz_file" | wc -l)
            if [[ "$LINE_COUNT" -gt 1 ]]; then
                # File has more than just the header, it's valid and complete
                ln -s "$(realpath "$gz_file")" "$VALID_FILES_DIR/$(basename "$gz_file")"
                echo "$VALID_FILES_DIR/$(basename "$gz_file")" >> "$VALID_FILE_LIST"
                echo "  [OK] Valid and complete: $(basename "$gz_file")"
            else
                # File has only header or is empty, even if gzip is valid
                mv "$gz_file" "$INVALID_FILES_DIR/"
                echo "  [SKIP] Header-only or empty: $(basename "$gz_file") -> moved to $INVALID_FILES_DIR"
            fi
        else
            # gzip integrity check failed (likely incomplete or corrupted file)
            mv "$gz_file" "$INVALID_FILES_DIR/"
            echo "  [FAIL] Incomplete/Corrupted (gzip check failed): $(basename "$gz_file") -> moved to $INVALID_FILES_DIR"
        fi
    fi
done

echo "Separation complete. Valid files are symlinked/copied to ${VALID_FILES_DIR}."
echo "List of valid files saved to ${VALID_FILE_LIST}."
echo "Incomplete/Header-only files moved to ${INVALID_FILES_DIR}."