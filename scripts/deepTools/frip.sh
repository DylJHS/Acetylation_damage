#!/bin/bash
#SBATCH --job-name=frip
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/frip-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/frip-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

echo "---------------------------------------------------------------------------------"
echo "Starting FRIP analysis script for $SLURM_ARRAY_TASK_ID"
echo "---------------------------------------------------------------------------------"

FOLDER="/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/${1}"

# Get the bam and the bed peak file
BAM_FILES=("${FOLDER}"/bams/*.bam)
SEACR_BEDS="${FOLDER}/seacr_results"
FRIP_DIR="${FOLDER}/frip_results"

mkdir -p "$SEACR_BEDS"  "$FRIP_DIR"

bam_file="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
base_name=$(basename "$bam_file" .bam)
peak_file="${SEACR_BEDS}/${base_name}_.relaxed.bed"

TEMP="${TEMP_DIR}/${base_name}"
mkdir -p "$TEMP"

echo -e "trying to find the peak file: \n$peak_file \n and the bam file: \n$bam_file"
echo "-------------------------------------------------------------------------------"

# Check for existance
if [[ ! -e "$bam_file" || ! -e "$peak_file" ]]; then
    echo "Error: BAM file or peak file does not exist."
    exit 1
fi

# Get the correctly mapped duplicate free bams
echo "Processing BAM file: $bam_file"
echo "-------------------------------------------------------------------------------"
samtools view -b -f 0x2 -F 0xF04  "$bam_file" > "${TEMP}/${base_name}_deduped.bam"
samtools view "${TEMP}/${base_name}_deduped.bam" | awk '{print $0}' | head -n 10
echo "-------------------------------------------------------------------------------"

# Sort the deduped bam file
echo "Sorting deduped BAM file..."
samtools sort -@ 8 -o "${TEMP}/${base_name}_deduped_sorted.bam" "${TEMP}/${base_name}_deduped.bam"
echo "-------------------------------------------------------------------------------"

# Convert to fragments 
echo "Converting BAM to fragments "
bedtools bamtobed -bedpe -i "${TEMP}/${base_name}_deduped_sorted.bam" | \
   awk 'BEGIN{OFS="\t"} $1==$4{ s=($2<$5?$2:$5); e=($3>$6?$3:$6); print $1,s,e }' | \
   sort -k1,1 -k2,2n -k3,3n \
   > "${TEMP}/${base_name}_fragments.bed"
awk '{print $0}' "${TEMP}/${base_name}_fragments.bed" | head -n 10
echo "-------------------------------------------------------------------------------"

# Compute the total fragments
echo "Calculating total fragments..."
TOTAL=$(wc -l < "${TEMP}/${base_name}_fragments.bed")
echo "Total fragments: $TOTAL"
echo "-------------------------------------------------------------------------------"

# Count the overlapping fragments
echo "Counting overlapping fragments with peaks..."
OVERLAP=$(bedtools intersect -u -a "${TEMP}/${base_name}_fragments.bed" -b "$peak_file" | wc -l)
echo "Overlapping fragments: $OVERLAP"
echo "-------------------------------------------------------------------------------"

# Compute the FRIP score
echo "Calculating FRIP score..."
awk -v o="$OVERLAP" -v t="$TOTAL" 'BEGIN{printf("%.6f\n", (t?o/t:0))}' > "${FRIP_DIR}/${base_name}_deduped_frip_score.txt"
echo "-------------------------------------------------------------------------------"

echo "FRIP score for $base_name: $(cat "${FRIP_DIR}/${base_name}_deduped_frip_score.txt")"
echo "---------------------------------------------------------------------------------"
echo "FRIP analysis completed for $base_name"
echo "---------------------------------------------------------------------------------"

# Clear the temporary folder
rm -r "${TEMP}"