#!/bin/bash
#SBATCH --job-name=plt_fngrprnt_tagged
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/plt_fngrprnt_tagged-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/plt_fngrprnt_tagged-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/envs/deeptools_env

# Define the paths
DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/fly_alignments/tagged"
OUTPUT_DIR="${DIR}/fingerprints"
mkdir -p "$OUTPUT_DIR"


OUTPUT_PNG="${DIR}/fingerprints/fingerprint.png"

echo "Generating fingerprint plot..."
plotFingerprint -b /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/fly_alignments/tagged/dedup/*.bam \
    --labels "control_003" "control_005" "H3K9ac_006" "H3K9ac_007" "H3K9ac_008" \
    --minMappingQuality 30 \
    --outRawCounts "${DIR}/fingerprints/out_fingerprint.txt" \
    --skipZeros \
    --plotFile "$OUTPUT_PNG" 
