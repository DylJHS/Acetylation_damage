#!/bin/bash
#SBATCH --job-name=bam_cvrg
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bam_cvrg-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bam_cvrg-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
<<<<<<< HEAD
#SBATCH --mem=5G
=======
#SBATCH --mem=16G
>>>>>>> 77ff64d29879b3ee7244d73e4d5e2b5e730ced04
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/envs/deeptools_env

# Define the species
<<<<<<< HEAD
# SPECIES="human"
SPECIES="fly"
=======
SPECIES="human"
# SPECIES = "fly"
>>>>>>> 77ff64d29879b3ee7244d73e4d5e2b5e730ced04

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
<<<<<<< HEAD
for INPUT_BAM in "$DIR"/dedup/*00[0-9]*_dedup.bam; do
=======
for INPUT_BAM in "$DIR"/dedup/*JAN-00?_dedup.bam; do
>>>>>>> 77ff64d29879b3ee7244d73e4d5e2b5e730ced04
    if [[ -f "$INPUT_BAM" ]]; then
        # Extract the base name of the BAM file
        BASE_NAME=$(basename "$INPUT_BAM" .bam)
        
        # Define the output file name
        OUTPUT_BAM="${OUTPUT_DIR}/${BASE_NAME}_coverage.bw"
        
        # Run bamCoverage
        echo "Running bamCoverage for $INPUT_BAM..."
<<<<<<< HEAD

        if [[ "$SPECIES" == "fly" ]]; then
            bamCoverage -b "$INPUT_BAM" \
                --outFileName "$OUTPUT_BAM" \
                --outFileFormat bigwig \
                --binSize 10 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize 142573017 \
                --ignoreForNormalization chrX \
                --extendReads
        else
            bamCoverage -b "$INPUT_BAM" \
                --outFileName "$OUTPUT_BAM" \
                --outFileFormat bigwig \
                --binSize 10 \
                --normalizeUsing RPGC \
                --effectiveGenomeSize 2913022398 \
                --ignoreForNormalization chrX \
                --extendReads
        fi
=======
        bamCoverage -b "$INPUT_BAM" \
            --outFileName "$OUTPUT_BAM" \
            --outFileFormat bigwig \
            --binSize 10 \
            --effectiveGenomeSize 2913022398 \
            --ignoreForNormalization chrX \
            --extendReads
>>>>>>> 77ff64d29879b3ee7244d73e4d5e2b5e730ced04
    else
        echo "No BAM files found in $DIR/dedup/"
    fi
done
