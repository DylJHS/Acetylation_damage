#!/bin/bash
#SBATCH --job-name=run_blastn
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/blastn_query-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/blastn_query-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-23
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/BLAST/BLAST_wkflw_config.sh

# Get the FASTA file for this array task
FASTA_FILES=($FASTA_DIR/*_subset.fasta)
FASTA_FILE="${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}"

# Sanity check
if [ ! -f "$FASTA_FILE" ]; then
    echo "FASTA file not found: $FASTA_FILE"
    exit 1
fi

# Loop through each subset FASTA file in the FASTA directory and run BLAST against each database.
base=$(basename "$FASTA_FILE" _subset.fasta)

for species in human drosophila; do
    # Select the appropriate database for each species.
    if [ "$species" = "human" ]; then
        db_name="${HUMAN_REF%.*}_db"
    elif [ "$species" = "drosophila" ]; then
        db_name="${DROS_REF%.*}_db"
    fi
    output_file="${RES_DIR}/${species}/${base}_${species}_blast_results.txt"
    echo "Running BLAST for $base against $species database:"
    echo "Using database: $db_name"
    blastn -query "$FASTA_FILE" \
        -task megablast \
        -perc_identity 95 \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -ungapped \
        -db "$db_name" \
        -num_threads 4 \
        -out "$output_file" \
        -evalue 1e-20 \
        -outfmt 6
done
