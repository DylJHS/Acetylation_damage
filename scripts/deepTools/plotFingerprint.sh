#!/bin/bash
#SBATCH --job-name=plt_fngrprnt
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/plt_fngrprnt-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/plt_fngrprnt-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Define the paths
DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/human_alignments"
OUTPUT_DIR="${DIR}/fingerprints"
mkdir -p "$OUTPUT_DIR"

# Loop over the BAM files in the directory
for INPUT_BAM in "$DIR"/dedup/*JAN-00?_dedup.bam; do
    if [[ -f "$INPUT_BAM" ]]; then
        echo "Processing: $INPUT_BAM"

        BASE=$(basename "$INPUT_BAM" .bam)

        OUTPUT_PNG="${DIR}/fingerprints/${BASE}_fingerprint.png"
        OUTPUT_PDF="${DIR}/fingerprints/${BASE}_fingerprint.pdf"

        echo "Generating fingerprint plot..."
        plotFingerprint -b "$INPUT_BAM" \
            --labels "H3K9ac" \
            --minMappingQuality 30 \
            --outRawCounts "${DIR}/fingerprints/${BASE}_fingerprint.txt" \
            --plotFile "$OUTPUT_PNG" \
            --plotTitle "Fingerprint for $BASE" \
            --dpi 300 \
            --plotHeight 5 \
            --plotWidth 7 \
            --outFileName "$OUTPUT_PDF"
    fi
done
