#!/bin/bash
# Directory containing BLAST result files
RESULTS_DIR="/home/djhs/fly_acetylation_damage/results"
# Number of sample reads used in each BLAST query (adjust if needed)
SAMPLE_SIZE=1000
# Path to the Python summary script
PY_SCRIPT="/home/djhs/fly_acetylation_damage/scripts/BLAST_res_summary.py"
# Master summary output file
OUTPUT_SUMMARY="$RESULTS_DIR/master_blast_summary.txt"

# Empty or create the master summary file
> "$OUTPUT_SUMMARY"

# Loop over each BLAST result file ending with *_blast_results.txt in the RESULTS_DIR
for file in "$RESULTS_DIR"/*_blast_results.txt; do
    # Extract the base filename (e.g., SCC-bulkChIC-UMC-JAN-003_R1_drosophila_blast_results.txt)
    base=$(basename "$file")
    echo "Processing $base ..." | tee -a "$OUTPUT_SUMMARY"
    # Call the Python script, passing the file name (the Python script expects the file to be in RESULTS_DIR)
    python3 "$PY_SCRIPT" "$base" "$SAMPLE_SIZE" >> "$OUTPUT_SUMMARY"
    echo "----------------------------------------" >> "$OUTPUT_SUMMARY"
done

echo "All summaries written to $OUTPUT_SUMMARY"
