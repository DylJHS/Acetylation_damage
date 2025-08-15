#!/bin/bash
#SBATCH --job-name=integration_test
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/inter_test-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/inter_test-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/envs/deeptools_env

# Define the paths
DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/fly_alignments/tagged"

# Define integration window (dm6 coordinates) 600bp wide window
echo -e "2L\t22245250\t22245350\tint_site" > integration.bed
    
# Open output file for writing
OUTFILE="$DIR/integration_counts.txt"
> "$OUTFILE"  # Clear if it exists

# Loop through BAM files
for bam in "$DIR"/dedup/003_tagged_dedup.bam "$DIR"/dedup/005_tagged_dedup.bam \
           "$DIR"/dedup/006_tagged_dedup.bam "$DIR"/dedup/007_tagged_dedup.bam; do

    # Step 1: Fragment stats
    # Sample name without path or extension
    sample=$(basename "$bam" .bam)

    # Step 1: count total primary, paired, non-supplementary alignments
    total=$(samtools view -c -F 0x904 "$bam")

    # Step 2: count reads overlapping the integration region
    count=$(bedtools coverage -a integration.bed -b "$bam" )

    echo -e "$sample\t$total\t$count" >> "$OUTFILE"

    # Step 2: Normalised coverage
    bamCoverage -b "$bam" \
        --outFileName "$DIR/$(basename "$bam" .bam)_RPKM_integration.bw" \
        --outFileFormat bigwig \
        --binSize 50 \
        --normalizeUsing RPKM \
        --effectiveGenomeSize 142573017
done
