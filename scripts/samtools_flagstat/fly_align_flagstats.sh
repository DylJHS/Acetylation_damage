#!/bin/bash
#SBATCH --job-name=fly_align_results
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/fly_align_results-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/fly_align_results-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Fix missing libcrypto issue for samtools
export LD_LIBRARY_PATH=/hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env/lib:$LD_LIBRARY_PATH


# Define paths
RESULTS_DIR="/hpc/shared/onco_janssen/dhaynessimmons/results/fly_acetylation_damage/fly_alignments"
OUTPUT_FILE="${RESULTS_DIR}/alignment_summary.txt"

# Ensure results directory exists
mkdir -p $RESULTS_DIR

# Clear the output file before writing new results
echo "Alignment Statistics Summary - $(date)" > $OUTPUT_FILE
echo "======================================" >> $OUTPUT_FILE

# Loop through BAM files and save stats
for file in ${RESULTS_DIR}/*.bam; do
    if [[ -f $file ]]; then
        echo "Processing: $file"
        echo "Stats for $file:" | tee -a $OUTPUT_FILE
        samtools flagstat "$file" | tee -a $OUTPUT_FILE
        echo "--------------------------------------" >> $OUTPUT_FILE
        echo "" >> "$OUTPUT_FILE"  # Adds a blank line to the output file
        echo 
    else
        echo "Warning: No BAM files found in $RESULTS_DIR" | tee -a $OUTPUT_FILE
    fi
done

echo "Alignment stats saved to: $OUTPUT_FILE"