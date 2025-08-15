#!/bin/bash
#SBATCH --job-name=gunzip
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/gunzip-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/gunzip-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/genomics_env_config.sh


FILE="$1"
if [[ ! -f "$FILE" ]]; then
    echo "File not found: $FILE"
    exit 1
fi 

filename=$(basename "$FILE")
multi_ext="${filename#*.}"    # gives 'tar.gz' for file.tar.gz
echo "$multi_ext"

if [[ $multi_ext == "tar.gz" ]]; then
    echo "Compressed file found: ${FILE}"
    echo "Unzipping..."
    tar -xzf "${FILE}" -C /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/data/
else
    echo "No compressed file found or unsupported format."
fi