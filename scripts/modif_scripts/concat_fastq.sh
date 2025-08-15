#!/bin/bash
#SBATCH --job-name=fastq_concat
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fastq_concat-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fastq_concat-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

## Concatenates gzipped fastq files in subfolders of a specified directory across lanes.

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/genomics_env_config.sh

# Find the correct directory 
DATA_DIR="${DATA_DIR}/SEND-bulkChIC-UMC-JAN-003-008-2"

# Get the subfolders from the DIR
folders=("${DATA_DIR}"/bulkChIC*)
if [[ ${#folders[@]} -eq 0 ]]; then
    echo "No subfolders found in $DATA_DIR"
    exit 1
fi  

# Loop over the subfolders
for subfolder in "${folders[@]}"; do
    sample_name=$(basename "$subfolder")
    echo -e "\nProcessing subfolder: $sample_name \n"
    # Get the fastq files 
    fastq_files=("${subfolder}"/*fastq.gz)
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        echo "No fastq files found in $subfolder"
        continue
    fi
    R1_fq_files=()
    R2_fq_files=()
   
    for f in "${fastq_files[@]}"; do
        if [[ "$f" == *"_R1_"* ]]; then
            R1_fq_files+=("$f")
        elif [[ "$f" == *"_R2_"* ]]; then
            R2_fq_files+=("$f")
        fi
    done

    # Deduplicate R1 and R2 files
    unique_R1_files=($(printf "%s\n" "${R1_fq_files[@]}" | sort -u))
    echo "Unique R1 files: ${unique_R1_files[@]}"
    unique_R2_files=($(printf "%s\n" "${R2_fq_files[@]}" | sort -u))
    echo "Unique R2 files: ${unique_R2_files[@]}"

    echo "Merging R1 files..."
    cat "${unique_R1_files[@]}" > "${subfolder}/${sample_name}_R1_merged.fastq.gz"
    echo "Merging R2 files..."
    cat "${unique_R2_files[@]}" > "${subfolder}/${sample_name}_R2_merged.fastq.gz"
    echo "Merged files created: ${sample_name}_R1_merged.fastq.gz and ${sample_name}_R2_merged.fastq.gz"

    
done