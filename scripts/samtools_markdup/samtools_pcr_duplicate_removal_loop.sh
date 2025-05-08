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
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Define the paths
DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/human_alignments"

# Loop over the BAM files in the directory
for INPUT_BAM in "$DIR"/orig/*JAN-00?.bam; do
    if [[ -f "$INPUT_BAM" ]]; then
        echo "Processing: $INPUT_BAM"

        BASE=$(basename "$INPUT_BAM" .bam)

        NAME_SORTED_BAM="${DIR}/sortd/${BASE}_name_sorted.bam"
        FIXMATE_BAM="${DIR}/fixmate/${BASE}_fixmate.bam"
        POS_SORTED_BAM="${DIR}/sortd/${BASE}_pos_sorted.bam"
        DEDUP_BAM="${DIR}/dedup/${BASE}_dedup.bam"
        DEDUP_INDEX="${DIR}/dedup/${BASE}_dedup.bam.bai"

        echo "Sorting the input BAM by read name..."
        samtools collate -@ 8 -o "$NAME_SORTED_BAM" "$INPUT_BAM"

        echo "Fixing mate information..."
        samtools fixmate -m -@ 8 "$NAME_SORTED_BAM" "$FIXMATE_BAM"

        echo "Sorting by coordinate..."
        samtools sort -@ 8 -o "$POS_SORTED_BAM" "$FIXMATE_BAM"

        echo "Removing duplicates..."
        samtools markdup -r -@ 8 "$POS_SORTED_BAM" "$DEDUP_BAM"

        echo "Indexing deduplicated BAM..."
        samtools index "$DEDUP_BAM" "$DEDUP_INDEX"

        echo "Duplicate removal completed successfully. \n"
    fi
done
