#!/bin/bash
# Directory containing BLAST result files
RESULTS_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/BLAST_results"
# Path to the Python summary script
PY_SCRIPT="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/summarise_blast_output.py"

# Number of sample reads used in each BLAST query (adjust if needed)
SAMPLE_SIZE=1000

for folder in "$RESULTS_DIR"; do 

    folder_base = "$(basename "$folder")"
    OUPUT_PATH = "$RESULTS_DIR/$folder_base"
    OUTPUT_SUMMARY = "$OUTPUT_PATH/${folder_base}_summary_blast_results.txt"

    # Empty or create the master summary file
    > "$OUTPUT_SUMMARY"

    # Loop over each BLAST result file ending with *_blast_results.txt in the RESULTS_DIR
    for file in "$RESULTS_DIR"/*_blast_results.txt; do
        # Extract the base filename (e.g., SCC-bulkChIC-UMC-JAN-003_R1_drosophila_blast_results.txt)
        base=$(basename "$file")
        echo "Processing $base ..." | tee -a "$OUTPUT_SUMMARY"
        # Call the Python script, passing the file name (the Python script expects the file to be in RESULTS_DIR)
        python3 "$PY_SCRIPT" "$file" "$SAMPLE_SIZE" "$OUTPUT_PATH">> "$OUTPUT_SUMMARY"
        echo "----------------------------------------" >> "$OUTPUT_SUMMARY"
    done

    echo "All summaries written to $OUTPUT_SUMMARY"
done 