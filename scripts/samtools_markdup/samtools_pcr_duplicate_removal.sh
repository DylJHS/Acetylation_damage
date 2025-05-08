#!/bin/bash
#SBATCH --job-name=dupli_rem
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/dupli_rem-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/dupli_rem-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
eval "$(/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin/conda shell.bash hook)"
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env


# Define base name and output directory
BASE="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/human_alignments/SCC-bulkChIC-UMC-JAN-003"
INPUT_BAM="${BASE}.bam"
NAME_SORTED_BAM="${BASE}_namesort.bam"
FIXMATE_BAM="${BASE}_fixmate.bam"
POS_SORTED_BAM="${BASE}_positionsort.bam"
DEDUP_BAM="${BASE}_dedup.bam"
DEDUP_INDEX="${DEDUP_BAM}.bai"

# Step 1: Sort by read name
echo "Sorting the input BAM by read name..."
samtools sort -n -@ 8 -o "$NAME_SORTED_BAM" "$INPUT_BAM"

# Step 2: Fix mate information
echo "Fixing mate information..."
samtools fixmate -@ 8 "$NAME_SORTED_BAM" "$FIXMATE_BAM"

# Step 3: Sort by genomic coordinate
echo "Sorting by coordinate..."
samtools sort -@ 8 -o "$POS_SORTED_BAM" "$FIXMATE_BAM"

# Step 4: Remove duplicates
echo "Removing duplicates..."
samtools markdup -r -@ 8 "$POS_SORTED_BAM" "$DEDUP_BAM"

# Step 5: Index the final BAM
echo "Indexing deduplicated BAM..."
samtools index "$DEDUP_BAM" "$DEDUP_INDEX"

echo "Duplicate removal completed successfully."
