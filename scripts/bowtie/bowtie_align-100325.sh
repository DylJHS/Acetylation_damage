#!/bin/bash
#SBATCH --job-name=bowtie2-align        # Job name
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/bowtie2-align-%j.out  # Log file
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/bowtie2-align-%j.err   # Error log
#SBATCH --time=12:00:00                 # Max runtime (adjust if needed)
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=8                # Number of CPU cores
#SBATCH --mem=64G                        # Memory allocation (adjust as needed)
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Define paths
REFERENCE_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage"
INDEX_PREFIX="${REFERENCE_DIR}/BDGP6"  # Prefix of the Bowtie2 index
TRIMMED_DIR="${REFERENCE_DIR}/fastq/trimmed"
OUTPUT_DIR="${REFERENCE_DIR}/aligned"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through all trimmed FASTQ pairs and align them
for R1 in ${TRIMMED_DIR}/*_val_1.fq.gz; do
    SAMPLE=$(basename $R1 _val_1.fq.gz)
    R2="${TRIMMED_DIR}/${SAMPLE}_val_2.fq.gz"
    OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE}.bam"
    OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE}.sam"

    echo "Processing sample: $SAMPLE"

    # Align reads with Bowtie2
    bowtie2 -x $INDEX_PREFIX -1 $R1 -2 $R2 --threads 8 -S $OUTPUT_SAM

    # Convert SAM to sorted BAM
    samtools view -@ 8 -bS $OUTPUT_SAM | samtools sort -@ 8 -o $OUTPUT_BAM

    # Index BAM file
    samtools index $OUTPUT_BAM

    # Remove SAM file to save space
    rm $OUTPUT_SAM

    echo "Finished aligning $SAMPLE"
done

echo "All alignments completed at $(date)"

