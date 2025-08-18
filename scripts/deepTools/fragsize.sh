#!/bin/bash
#SBATCH --job-name=fragsize
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fragsize-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fragsize-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

BAM_FOLDER="$TAGGED_RAW_ALIGNMENT_DIR_BAMS"
OUT_FOLDER="$TAGGED_RAW_ALIGNMENT_DIR"/frag_size
mkdir -p "$OUT_FOLDER"

echo "Processing BAM files in $BAM_FOLDER"

bamPEFragmentSize \
    -b "$BAM_FOLDER"/*.bam \
    --histogram "$OUT_FOLDER/frag_size_histogram.pdf" \
    --maxFragmentLength 200 \
    --numberOfProcessors 8 \
    --samplesLabel ctr1 ctr2 ctr3 dsb1 dsb2 dsb3 

echo "Finished processing"
done