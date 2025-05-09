#!/bin/bash
#SBATCH --job-name=bam_cvrg
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bam_cvrg-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bam_cvrg-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/envs/deeptools_env

# Define the species
SPECIES="human"
# SPECIES = "fly"

# Define the paths
if [[ "$SPECIES" == "fly" ]]; then
    DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/fly_alignments/tagged"
else
    # Assuming human if not fly
    DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/human_alignments"
fi

OUTPUT_DIR="${DIR}/coverage_files"
mkdir -p "$OUTPUT_DIR"

# Loop over the BAM files in the directory
for INPUT_BAM in "$DIR"/dedup/*JAN-00?_dedup.bam; do
    if [[ -f "$INPUT_BAM" ]]; then
        # Extract the base name of the BAM file
        BASE_NAME=$(basename "$INPUT_BAM" .bam)
        
        # Define the output file name
        OUTPUT_BAM="${OUTPUT_DIR}/${BASE_NAME}_coverage.bw"
        
        # Run bamCoverage
        echo "Running bamCoverage for $INPUT_BAM..."
        bamCoverage -b "$INPUT_BAM" \
            --outFileName "$OUTPUT_BAM" \
            --outFileFormat bigwig \
            --binSize 10 \
            --effectiveGenomeSize 2913022398 \
            --ignoreForNormalization chrX \
            --extendReads
    else
        echo "No BAM files found in $DIR/dedup/"
    fi
done
