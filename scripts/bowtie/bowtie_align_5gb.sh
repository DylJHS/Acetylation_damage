#!/bin/bash
#SBATCH --job-name=bowtie2-align
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/bowtie2-align-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/bowtie2-align-%j.err
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

# Define paths
REFERENCE_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage"
INDEX_PREFIX="${REFERENCE_DIR}/BDGP6"
TRIMMED_DIR="${REFERENCE_DIR}/fastq/trimmed"
OUTPUT_DIR="${REFERENCE_DIR}/aligned"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through R1 reads and find the correct R2 pair
for R1 in ${TRIMMED_DIR}/*_R1_001_val_1.fq.gz; do
    SAMPLE=$(basename $R1 _R1_001_val_1.fq.gz)
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"

    # Check if R2 exists before running Bowtie2
    if [[ ! -f $R2 ]]; then
        echo "Error: R2 file missing for $SAMPLE"
        continue
    fi

    OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE}.bam"
    OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE}.sam"

    echo "Processing: $SAMPLE"

    # Run Bowtie2 alignment
    bowtie2 -x $INDEX_PREFIX -1 $R1 -2 $R2 --threads 8 -S $OUTPUT_SAM

    # Convert and sort BAM
    samtools view -@ 8 -bS $OUTPUT_SAM | samtools sort -@ 8 -o $OUTPUT_BAM
    samtools index $OUTPUT_BAM

    # Remove SAM to save space
    rm $OUTPUT_SAM

    echo "Finished aligning $SAMPLE"
done

echo "All alignments completed at $(date)"
