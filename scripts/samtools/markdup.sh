#!/bin/bash
#SBATCH --job-name=dupli_rem
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/dupli_rem-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/dupli_rem-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/wkflw_config.sh

# get the arguments for species
if [[ "$1" == "-human" ]]; then
    BAM_DIR="$HUMAN_ALIGN_BOWTIE_DIR"
    OUTPUT_PATH="$HUMAN_ALIGN_DIR"
elif [[ "$1" == "-drosophila" ]]; then
    BAM_DIR="$DROS_ALIGN_BOWTIE_DIR"
    OUTPUT_PATH="$DROS_ALIGN_DIR"
else
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi

OUTPUT_DIR="${OUTPUT_PATH}"

BAM_LIST=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_LIST[$SLURM_ARRAY_TASK_ID]}"

BASE=$(basename "$BAM_FILE" .bam)

echo "Processing: $BAM_FILE"
echo "-----------------------------------------------------------"

DEDUP_BAM="${OUTPUT_DIR}/dedup/${BASE}_dedup.bam"

# Create directories if they do not exist
mkdir -p "${OUTPUT_DIR}/dedup"

if [[ "$2" == "-r" ]]; then
    echo "Removing duplicates..."
    echo "-----------------------------------------------------------"
    samtools collate -@ 8 -O -u  -T "${TEMP_DIR}/collate_${BASE}" "$BAM_FILE" | \
        samtools fixmate -m -@ 8 -u - - | \
        samtools sort -@ 8 -u - | \
        samtools markdup -r -@ 8 - "$DEDUP_BAM"
else
    echo "Marking duplicates without removal..."
    echo "-----------------------------------------------------------"
    samtools collate -@ 8 -O -u -T "${TEMP_DIR}/collate_${BASE}" "$BAM_FILE" | \
        samtools fixmate -m -@ 8 -u - - | \
        samtools sort -@ 8 -u - | \
        samtools markdup -@ 8 - "$DEDUP_BAM"
fi

echo "Finalizing BAM file..."
echo "-----------------------------------------------------------"

echo "Removing temporary files..."
rm -f ${TEMP_DIR}/*_${BASE}*