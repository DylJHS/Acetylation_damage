#!/bin/bash
#SBATCH --job-name=duplicheck
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/duplicheck-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/duplicheck-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/genomics_env_config.sh


# get the arguments for species
if [[ "$1" == "-human" ]]; then
    BAM_DIR="$HUMAN_ALIGN_BOWTIE_DIR"
    DEDUP_BAM_DIR="$HUMAN_DEDUP_DIR"
elif [[ "$1" == "-drosophila" ]]; then
    BAM_DIR="$DROS_ALIGN_BOWTIE_DIR"
    DEDUP_BAM_DIR="$DROS_DEDUP_DIR"
elif [[ "$1" == "-tagged" ]]; then
    BAM_DIR="$TAGGED_ALIGNMENT_DIR"
    DEDUP_BAM_DIR="$TAGGED_DEDUP_DIR"
else
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi

output_file="${DEDUP_BAM_DIR}/duplicate_check.txt"

for bam in $BAM_DIR/*.bam; do
    base=$(basename $bam .bam)
    dedup_bam="${DEDUP_BAM_DIR}/${base}_dedup.bam"

    if [ ! -f $dedup_bam ]; then
        echo "$dedup_bam does not exist"
        continue
    fi

    echo "Processing: $base" >> $output_file
    echo "-----------------------------------------------------------" >> $output_file

    echo " Checking for:" >> $output_file
    echo "for pre duplicate reads:" >> $output_file
    samtools view -c $bam >> $output_file
    echo "for post duplicate reads:" >> $output_file
    samtools view -c $dedup_bam >> $output_file

    echo -e "----------------------------------------------------------- \n\n" >> $output_file
done