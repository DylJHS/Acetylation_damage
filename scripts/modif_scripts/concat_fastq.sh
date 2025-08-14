#!/bin/bash
#SBATCH --job-name=fastq_concat
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/fastq_concat-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/fastq_concat-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/genomics_env_config.sh

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
    # Get the fastq files 
    fastq_files=("${subfolder}"/*fastq.gz)
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        echo "No fastq files found in $subfolder"
        continue
    fi
    # Get the unique lanes
    unique_lanes=()
    for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq); do 
        echo "Merging R1"

        cat "$i"_L00*_R1_001.fastq.gz > "$i"_ME_L001_R1_001.fastq.gz

        echo "Merging R2"

        cat "$i"_L00*_R2_001.fastq.gz > "$i"_ME_L001_R2_001.fastq.gz

    done
done