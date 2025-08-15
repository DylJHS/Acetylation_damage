#!/bin/bash
#SBATCH --job-name=bai_indexing
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/bai_indexing-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/bai_indexing-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/genomics_env_config.sh

# get the arguments for species


for SPECIES in "tagged" "drosophila"; do
    if [[ $SPECIES == "drosophila" ]]; then
        BAM_DIR="$DROS_ALIGN_BOWTIE_DIR"
        DEDUP_PATH="$DROS_DEDUP_DIR"
    elif [[ "$SPECIES" == "human" ]]; then
        BAM_DIR="$HUMAN_ALIGN_BOWTIE_DIR"
        DEDUP_PATH="$HUMAN_DEDUP_DIR"
    elif [[ $SPECIES == "tagged" ]]; then
        BAM_DIR="$TAGGED_BAM_ORI_DIR"
        DEDUP_PATH="$TAGGED_DEDUP_DIR"
    fi
    echo "-----------------------------------------------------------"
    echo "Processing $SPECIES BAM files..."
    echo "-----------------------------------------------------------"

    BAM_LIST=("$DEDUP_PATH"/*.bam)
    BAM_FILE="${BAM_LIST[$SLURM_ARRAY_TASK_ID]}"

    BASE=$(basename "$BAM_FILE" .bam)
    DEDUP_BAI="${DEDUP_PATH}/${BASE}_dedup.bai"


    echo "Processing: $BASE"
    echo "-----------------------------------------------------------"


    samtools index -@ 8 "$BAM_FILE" "$DEDUP_BAI"

    echo "Finalised BAI file: $DEDUP_BAI"
    echo "-----------------------------------------------------------"
done