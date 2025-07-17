#!/bin/bash
#SBATCH --job-name=generate_fasta_subsets
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/generate_fasta_subsets-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/generate_fasta_subsets-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

# SCC folders
SCC_FOLDERS=("$DATA_DIR"/bulkChIC-UMC-JAN-*)
SCC_FOLDER="${SCC_FOLDERS[$SLURM_ARRAY_TASK_ID]}"

# Check
if [ ! -d "$SCC_FOLDER" ]; then
    echo "SCC folder not found: $SCC_FOLDER"
    exit 1
fi

# Process each gzipped FASTQ file for R1 and R2 in the folder
for file in "$SCC_FOLDER"/*_R[12]_*.fastq.gz; do
    # Check if file exists (in case pattern doesn't match any files)
    [ -e "$file" ] || continue
    echo "  Processing file: $file"

    # Extract FASTQ short ID
    fastq_base=$(basename "$file" .fastq.gz)
    short_id=$(echo "$fastq_base" | sed -E 's/^.*(JAN-[0-9]+)_.*(L[0-9]+)_R([12])_.*/\1_\2_R\3/')
    
    OUTPUT_FASTA="$FASTA_DIR/${short_id}_subset.fasta"
    echo "    Output FASTA: $OUTPUT_FASTA"

    # Decompress the file, sample reads, and convert to FASTA in one pipeline
    zcat "$file" | seqtk sample -s100 - "$SAMPLE_SIZE" | seqtk seq -a > "$OUTPUT_FASTA"
done