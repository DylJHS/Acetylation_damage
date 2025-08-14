#!/bin/bash
#SBATCH --job-name=frip
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/frip-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/frip-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/deeptools_env_config.sh

echo "---------------------------------------------------------------------------------"
echo "Starting FRIP analysis script for $SLURM_ARRAY_TASK_ID"
echo "---------------------------------------------------------------------------------"

# Define paths
if [[ "$1" == "-drosophila" ]]; then
    RESULTS_DIR="$DROS_ALIGNMENT_DIR"
elif [[ "$1" == "-human" ]]; then
    RESULTS_DIR="$HUMAN_ALIGNMENT_DIR"
elif [[ "$1" == "-tagged" ]]; then
    RESULTS_DIR="$TAGGED_ALIGNMENT_DIR"
else
    echo "Error: Invalid species or alignment type specified."
    exit 1
fi

# Has to be run on deduped alignments
if [[ "$2" == "-r" ]]; then
    FRAG_BEDS="$RESULTS_DIR"/deduped_alignments/removed/fragments
    BAM_FLDR="$RESULTS_DIR"/deduped_alignments/removed/bams
elif [[ "$2" == "-m" ]]; then
    FRAG_BEDS="$RESULTS_DIR"/deduped_alignments/marked/fragments
    BAM_FLDR="$RESULTS_DIR"/deduped_alignments/removed/bams
else
    echo "Error: Invalid deduplication type specified."
    exit 1
fi


# Get the bam and the bed peak file
BAM_FILES=("${BAM_FLDR}"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
if [[ "$2" == "-r" ]]; then 
    base_name=$(basename "$BAM_FILE" -r.bam)
else 
    base_name=$(basename "$BAM_FILE" .bam)
fi

peak_file="${FRAG_BEDS}/${base_name}.fragments.bedgraph"

echo -e "trying to find the peak file: \n$peak_file \n and the bam file: \n$BAM_FILE"
echo "-------------------------------------------------------------------------------"

# Create the output file
peaks_in="${FRAG_BEDS}/../bedgraph/${base_name}_frip.txt"

# Check that they both exist
if [[ ! -e "$BAM_FILE" || ! -e "$peak_file" ]]; then
    echo "Error: BAM file or peak file does not exist."
    exit 1
fi

# convert the paired end bam to bed in bedpe format
echo "Conversion of the BAM file to BED format" 
echo "------------------------------------------------------------------------------"
bedtools bamtobed -i "$BAM_FILE" -bedpe | \
awk '$1==$4 {
    s = ()}


# # Run bedtools intersect 
# echo -e "Running bedtools intersect for: \n$BAM_FILE \nand with peak file: \n$peak_file"
# echo "------------------------------------------------------------------------------"

# bedtools intersect -a "$BAM_FILE" -b "$peak_file" -bed > "$peaks_in"
# echo "Finished running bedtools intersect" 
# echo "------------------------------------------------------------------------------"

# Calculate FRIP
