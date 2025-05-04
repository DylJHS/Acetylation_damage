#!/bin/bash
# Directory containing the SCC folders
DATA_DIR="/home/djhs/fly_acetylation_damage/data"
# Directory where FASTA subset files will be saved
FASTA_DIR="$DATA_DIR/FASTA"
# Create FASTA directory if it doesn't exist
mkdir -p "$FASTA_DIR"

# Number of reads to sample from each file
SAMPLE_SIZE=1000

# Loop through each SCC folder (assumed to have names like SCC-bulkChIC-UMC-JAN-*)
for folder in "$DATA_DIR"/SCC-bulkChIC-UMC-JAN-*; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
        # Process each gzipped FASTQ file for R1 and R2 in the folder
        for file in "$folder"/*_R[12]_*.fastq.gz; do
            # Check if file exists (in case pattern doesn't match any files)
            [ -e "$file" ] || continue
            echo "  Processing file: $file"
            # Get the folder's base name (e.g., SCC-bulkChIC-UMC-JAN-003)
            folder_id=$(basename "$folder")
            # Determine the read type (R1 or R2) from the filename
            if [[ "$file" == *"_R1_"* ]]; then
                read_tag="R1"
            elif [[ "$file" == *"_R2_"* ]]; then
                read_tag="R2"
            else
                read_tag="unknown"
            fi
            # Construct the output FASTA file name; e.g., SCC-bulkChIC-UMC-JAN-003_R1_subset.fasta
            OUTPUT_FASTA="$FASTA_DIR/${folder_id}_${read_tag}_subset.fasta"
            echo "    Output FASTA: $OUTPUT_FASTA"
            # Decompress the file, sample reads, and convert to FASTA in one pipeline
            zcat "$file" | seqtk sample -s100 - "$SAMPLE_SIZE" | seqtk seq -a > "$OUTPUT_FASTA"
        done
    fi
done
