#!/bin/bash
#SBATCH --job-name=bam_cvrg
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/bam_cvrg-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/bam_cvrg-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

# Set directories based on the alignment type
BAM_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/${2}"
BAM_FILES=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
FILE_NAME=$(basename "${BAM_FILE}" ".bam")

if [[ "$BAM_DIR" == *drosophila* || "$BAM_DIR" == *tagged* ]]; then
    effectivegsize=142573017
elif [[ "$BAM_DIR" == *human* ]]; then 
    effectivegsize=2913022398
fi

if [[ "$1" == "-bedgraph" ]]; then
    BED_DIR="${BAM_DIR}/../bed"
    OUTPUT_FILE="${BED_DIR}/${FILE_NAME}_coverage.bedgraph"
elif [[ "$1" == "-bigwig" ]]; then
    BIGWIG_DIR="${BAM_DIR}/../bigwig"
    OUTPUT_FILE="${BIGWIG_DIR}/${FILE_NAME}_coverage.bw"
else
    echo "Error: Please specigy valid output format -bedgraph or -bigwig."
fi

# Run bamCoverage
echo "Running bamCoverage for $FILE_NAME..."

bamCoverage -b "$BAM_FILE" \
    --outFileName "$OUTPUT_FILE" \
    --outFileFormat "${1#-}" \
    --binSize 25 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "${effectivegsize}" \
    --ignoreForNormalization chrX \
    --ignoreDuplicates \
    --centerReads \
    --extendReads

echo "bamCoverage finished"