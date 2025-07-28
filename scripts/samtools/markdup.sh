#!/bin/bash
#SBATCH --job-name=dupli_rem
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/dupli_rem-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/dupli_rem-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/wkflw_config.sh

# get the arguments for species
if [[ "$1" == "human" ]]; then
    BAM_DIR="$HUMAN_ALIGN_BOWTIE_DIR"
    OUTPUT_PATH="$HUMAN_ALIGN_DIR"
elif [[ "$1" == "drosophila" ]]; then
    BAM_DIR="$DROS_ALIGN_BOWTIE_DIR"
    OUTPUT_PATH="$DROS_ALIGN_DIR"
else
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi

OUTPUT_DIR="${OUTPUT_PATH}"

BAM_LIST=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_LIST[$SLURM_ARRAY_TASK_ID]}"

echo "Processing: $BAM_FILE"

BASE=$(basename "$BAM_FILE" .bam)

NAME_SORTED_BAM="${OUTPUT_DIR}/sortd/${BASE}_name_sorted.bam"
FIXMATE_BAM="${OUTPUT_DIR}/fixmate/${BASE}_fixmate.bam"
POS_SORTED_BAM="${OUTPUT_DIR}/sortd/${BASE}_pos_sorted.bam"
DEDUP_BAM="${OUTPUT_DIR}/dedup/${BASE}_dedup.bam"
DEDUP_INDEX="${OUTPUT_DIR}/dedup/${BASE}_dedup.bam.bai"

# Create directories if they do not exist
mkdir -p "${OUTPUT_DIR}/sortd" "${OUTPUT_DIR}/fixmate" "${OUTPUT_DIR}/dedup"

echo "Sorting the input BAM by read name..."
samtools collate -@ 4 -o "$NAME_SORTED_BAM" "$BAM_FILE"

echo "Fixing mate information..."
samtools fixmate -m -@ 4 "$NAME_SORTED_BAM" "$FIXMATE_BAM"

echo "Sorting by coordinate..."
samtools sort -@ 4 -o "$POS_SORTED_BAM" "$FIXMATE_BAM"

echo "Removing duplicates..."
samtools markdup -r -@ 4 "$POS_SORTED_BAM" "$DEDUP_BAM"

echo "Indexing deduplicated BAM..."
samtools index "$DEDUP_BAM" "$DEDUP_INDEX"

echo -e "Duplicate removal completed successfully. \n"

