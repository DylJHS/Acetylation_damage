#!/bin/bash
#SBATCH --job-name=bowtie2-build       # Job name
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/bowtie2-build-%j.out  # Log file (%j inserts job ID)
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/bowtie2-build_%j.err   # Error log
#SBATCH --time=06:00:00                # Max runtime (adjust if needed)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores (adjust as needed)
#SBATCH --mem=32G                       # Memory allocation (adjust as needed)
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/genomics_env_config.sh

for file in ${HUMAN_GEN_FA} ${DROS_GEN_FA}; do
    if [[ ${file} == "*.gz" ]]; then
        echo "Compressed file found: ${file}"
        echo "Unzipping..."
        gunzip "${file}" 
    fi

    if [[ ${file} == ${HUMAN_GEN_FA} ]]; then
        prefix="${DROS_INDEX}/human"
    elif [[ ${file} == ${DROS_GEN_FA} ]]; then
        prefix="${DROS_INDEX}/drosophila"
    fi

    # Run Bowtie2 index build
    echo "Starting Bowtie2 index build at $(date)"
    bowtie2-build --threads 8 $file ${prefix}_bowtie2_index
    echo "Bowtie2 index build completed at $(date)"
done

