#!/bin/bash
#SBATCH --job-name=bam_cvrg
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bam_cvrg-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bam_cvrg-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/deeptools_env_config.sh

# get the arguments for species
if [[ "$1" == "-human" ]]; then
    BAM_DIR="$HUMAN_DEDUP_DIR"
    BED_DIR="$HUMAN_DEDUP_BEDGRAPH"
    BIGWIG_DIR="$HUMAN_DEDUP_BIGWIG"
    effectivegsize=2913022398
elif [[ "$1" == "-drosophila" ]]; then
    BAM_DIR="$DROS_DEDUP_DIR"
    BED_DIR="$DROS_DEDUP_BEDGRAPH"
    BIGWIG_DIR="$DROS_DEDUP_BIGWIG"
    effectivegsize=142573017
elif [[ "$1" == "-tagged" ]]; then
    BAM_DIR="$TAGGED_ALIGNMENT_DIR"
    BED_DIR="$TAGGED_DEDUP_BEDGRAPH"
    BIGWIG_DIR="$TAGGED_DEDUP_BIGWIG"
    effectivegsize=142573017
else
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi

BAM_FILES=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
FILE_NAME=$(basename "${BAM_FILE}")

if [[ "$2" == "-bedgraph" ]]; then
    OUTPUT_FILE="${BED_DIR}/${FILE_NAME}_coverage.bedgraph"
elif [[ "$2" == "-bigwig" ]]; then
    OUTPUT_FILE="${BIGWIG_DIR}/${FILE_NAME}_coverage.bw"
else
    echo "Error: Please specigy valid output format -bedgraph or -bigwig."
fi

# Run bamCoverage
echo "Running bamCoverage for $FILE_NAME..."

bamCoverage -b "$BAM_FILE" \
    --outFileName "$OUTPUT_FILE" \
    --outFileFormat "${2#-}" \
    --binSize 50 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "${effectivegsize}" \
    --ignoreForNormalization chrX \
    --extendReads

echo "bamCoverage finished"