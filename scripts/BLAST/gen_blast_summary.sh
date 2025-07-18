#!/bin/bash
#SBATCH --job-name=gen_blast_summary
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/gen_blast_summary-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/gen_blast_summary-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

# Summarise the BLAST results using the Python script
for folder in "$RES_DIR"/*; do 
    # Check if the item is a directory
    if [ ! -d "$folder" ]; then
        echo "Skipping $folder, not a directory."
        continue
    fi
    folder_base="$(basename "$folder")"
    echo " ----------------------------------------"
    echo "Processing folder: $folder_base"
    BLAST_PATH="$RES_DIR/$folder_base"
    OUTPUT_SUMMARY="$RES_DIR/${folder_base}_summary_blast_results.txt"

    # Empty or create the master summary file
    > "$OUTPUT_SUMMARY"

    # Loop over each BLAST result file ending with *_blast_results.txt in the RES_DIR
    for file in "$BLAST_PATH"/JAN*_blast_results.txt; do
        # Extract the base filename (e.g., SCC-bulkChIC-UMC-JAN-003_R1_drosophila_blast_results.txt)
        base=$(basename "$file")
        echo "Processing file $base ..." | tee -a "$OUTPUT_SUMMARY"
        # Call the Python script, passing the file name (the Python script expects the file to be in RES_DIR)
        python3 "$PY_SCRIPT" "$file" "$SAMPLE_SIZE" >> "$OUTPUT_SUMMARY"
        echo "----------------------------------------" >> "$OUTPUT_SUMMARY"
    done

    echo "All summaries written to $OUTPUT_SUMMARY"
done 