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
    DEDUP_PATH="$HUMAN_DEDUP_DIR"
elif [[ "$1" == "-drosophila" ]]; then
    BAM_DIR="$DROS_ALIGN_BOWTIE_DIR"
    DEDUP_PATH="$DROS_DEDUP_DIR"
else
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi
echo "Processing ${1#-} BAM files..."
echo "-----------------------------------------------------------"

BAM_LIST=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_LIST[$SLURM_ARRAY_TASK_ID]}"

BASE=$(basename "$BAM_FILE" .bam)
DEDUP_BAM="${DEDUP_PATH}/${BASE}_dedup.bam"
DEDUP_STATS="${DEDUP_PATH}/deduped_stats"
TEMP_FLDER="${TEMP_DIR}/${BASE}_${1#-}"

echo "Processing: $BASE"
echo "-----------------------------------------------------------"

# Create directories if they do not exist
mkdir -p "${TEMP_FLDER}"

if [[ "$2" == "-r" ]]; then
    echo "Removing duplicates..."
    echo "-----------------------------------------------------------"
    samtools collate -@ 8 -O -u  -T "${TEMP_FLDER}/collate_${BASE}" "$BAM_FILE" | \
        samtools fixmate -m -@ 8 -u - - | \
        samtools sort -@ 8 -u -T "${TEMP_FLDER}/sort_${BASE}" - | \
        samtools markdup -@ 8 -r -f "${DEDUP_STATS}" -T "${TEMP_FLDER}/mrkd_${BASE}" - "$DEDUP_BAM"
else
    echo "Marking duplicates without removal..."
    echo "-----------------------------------------------------------"
    samtools collate -@ 8 -O -u -T "${TEMP_FLDER}/collate_${BASE}" "$BAM_FILE" | \
        samtools fixmate -m -@ 8 -u - - | \
        samtools sort -@ 8 -u -T "${TEMP_FLDER}/sort_${BASE}" - | \
        samtools markdup -@ 8 -f "${DEDUP_STATS}" -T "${TEMP_FLDER}/mrkd_${BASE}" - "$DEDUP_BAM"
fi

echo "Finalizing BAM file..."
echo "-----------------------------------------------------------"

echo "Removing temporary files..."
rm -rf ${TEMP_FLDER}