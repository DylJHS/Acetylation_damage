#!/bin/bash
#SBATCH --job-name=integration_test
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/inter_test-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/inter_test-%j.err
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
DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/fly_alignments/tagged"


# Define integration window (dm6 coordinates) 600bp wide window
echo -e "chr2L\t22245011\t22245611\tint_site" > integration.bed    

# total mapped fly fragments in each BAM (QC-pass)
for bam in "$DIR"/dedup/003_tagged_dedup.bam "$DIR"/dedup/005_tagged_dedup.bam \
            "$DIR"/dedup/006_tagged_dedup.bam "$DIR"/dedup/007_tagged_dedup.bam; do
    total=$(samtools view -c -F 0x904 $bam)               # paired, primary, non-supplementary
    inSite=$(bedtools coverage -a integration.bed -b $bam -counts \
             | awk '{print $4}')
    echo -e "$bam\t$total\t$inSite"

done > /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/fly_alignments/tagged/integration_counts.txt


# Create scaled bigwig for each BAM
for bam in "$DIR"/dedup/003_tagged_dedup.bam "$DIR"/dedup/005_tagged_dedup.bam \
            "$DIR"/dedup/006_tagged_dedup.bam "$DIR"/dedup/007_tagged_dedup.bam; do
    bamCoverage -b $bam \
        --outFileName "$DIR/$(basename "$bam" .bam)_RPKM_integration.bw" \
        --outFileFormat bigwig \
        --binSize 50 \
        --normalizeUsing RPKM \
        --effectiveGenomeSize 142573017 
    
done