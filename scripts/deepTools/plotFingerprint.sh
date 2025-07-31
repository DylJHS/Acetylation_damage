#!/bin/bash
#SBATCH --job-name=fingerprint
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/fingerprint-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/fingerprint-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/deeptools_env_config.sh

# Define the paths
# Define paths
if [[ "$1" == "-drosophila" ]]; then
    RESULTS_DIR="$DROS_DEDUP_DIR"
elif [[ "$1" == "-human" ]]; then
    RESULTS_DIR="$HUMAN_DEDUP_DIR"
elif [[ "$1" == "-tagged" ]]; then
    RESULTS_DIR="$TAGGED_ALIGNMENT_DIR"
else
    echo "--------------------------------------------------"
    echo "Error: Invalid alignment type specified."
    exit 1
fi

OUTPUT="${RESULTS_DIR}/fingerprint"
mkdir -p "${OUTPUT}"

echo "--------------------------------------------------"
echo "Generating fingerprint plot for ${1#-} files"

plotFingerprint \
    -b "${RESULTS_DIR}"/*.bam \
    --smartLabels \
    --minMappingQuality 30 \
    --outRawCounts "${OUTPUT}/out_fingerprint.txt" \
    -bs 10 \
    --skipZeros \
    -v \
    -r chr2L:22245050:22245550 \
    --plotFile "${OUTPUT}/fingerprint.png" 

echo "--------------------------------------------------"
echo "finshed running plotFingerprint, output saved to ${OUTPUT}"