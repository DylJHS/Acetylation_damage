#!/bin/bash
#SBATCH --job-name=fasta_subset
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/fasta_subset-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/fasta_subset-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Directory containing the SCC folders
DATA_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage"
# Directory where FASTA subset files will be saved
FASTA_DIR="$DATA_DIR/FASTA_subsets"
# Create FASTA directory if it doesn't exist
mkdir -p "$FASTA_DIR"

# Number of reads to sample from each file
SAMPLE_SIZE=1000

# Loop through each SCC folder
for folder in "$DATA_DIR"/SCC-bulkChIC-UMC-JAN-*; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
        # Process each gzipped FASTQ file for R1 and R2 in the folder
        for file in "$folder"/*_R[12]_*.fastq.gz; do
            # Check if file exists
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
            # Decompress the file, sample reads, and convert to FASTA 
            zcat "$file" | seqtk sample -s100 - "$SAMPLE_SIZE" | seqtk seq -a > "$OUTPUT_FASTA"
        done
    fi
done