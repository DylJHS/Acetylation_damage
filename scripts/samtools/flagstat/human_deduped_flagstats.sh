#!/bin/bash
#SBATCH --job-name=flagstat_human_dedup
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/flagstat_human_dedup-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/flagstat_human_dedup-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Define the paths
DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/human_alignments"
INPUT_DIR="${DIR}/dedup"
OUTPUT_FILE="${DIR}/dedup2_flagstats.txt"


# Loop through BAM files and save stats
for file in ${INPUT_DIR}/*.bam; do
    if [[ -f $file ]]; then
        echo "Processing: $file"
        echo "Stats for $file:" | tee -a $OUTPUT_FILE
        samtools flagstat "$file" | tee -a $OUTPUT_FILE
        echo "--------------------------------------" >> $OUTPUT_FILE
        echo "" >> "$OUTPUT_FILE"  # Adds a blank line to the output file
        echo 
    else
        echo "Warning: No BAM files found in $INPUT_DIR" | tee -a $OUTPUT_FILE
    fi
done

echo "Alignment stats saved to: $OUTPUT_FILE"

