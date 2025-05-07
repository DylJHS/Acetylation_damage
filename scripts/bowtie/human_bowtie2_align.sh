#!/bin/bash
#SBATCH --job-name=bowtie2-align
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/human_bowtie2-align-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/human_bowtie2-align-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G  # Adjusted for better efficiency
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Fix missing libcrypto issue for samtools
export LD_LIBRARY_PATH=/hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env/lib:$LD_LIBRARY_PATH

# Define paths
REFERENCE_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage"
INDEX_PREFIX="${REFERENCE_DIR}/GRCh38/GRCh38"
OUTPUT_DIR="/hpc/shared/onco_janssen/dhaynessimmons/results/fly_acetylation_damage/human_alignments"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through all SCC folders
for SCC_DIR in ${REFERENCE_DIR}/SCC-bulkChIC-UMC-JAN-*; do
    if [[ -d $SCC_DIR ]]; then
        SAMPLE_ID=$(basename $SCC_DIR)

        # Identify R1 and R2 files
        R1=$(ls ${SCC_DIR}/*_R1_001.fastq.gz 2>/dev/null | head -n 1)
        R2=$(ls ${SCC_DIR}/*_R2_001.fastq.gz 2>/dev/null | head -n 1)

        # Ensure both files exist before proceeding
        if [[ -z "$R1" || -z "$R2" ]]; then
            echo "Error: Missing R1 or R2 for ${SAMPLE_ID}" >> "${OUTPUT_DIR}/alignment_errors.log"
            continue
        fi

        # Define output file paths
        OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.bam"
        OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_ID}.sam"

        echo "Processing: ${SAMPLE_ID}"

        # Run Bowtie2 alignment
        bowtie2 -x $INDEX_PREFIX -1 $R1 -2 $R2 --threads 8 -S $OUTPUT_SAM

        # Convert SAM to sorted BAM
        samtools view -@ 8 -bS $OUTPUT_SAM | samtools sort -@ 8 -o $OUTPUT_BAM
        samtools index $OUTPUT_BAM


        echo "Finished aligning ${SAMPLE_ID}"
    fi
done

echo "All alignments completed at $(date)"
