#!/bin/bash
#SBATCH --job-name=fastqc_run
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fastqc_run-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fastqc_run-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/BLAST/BLAST_wkflw_config.sh

# SCC folders
SCC_FOLDERS=("$DATA_DIR"/bulkChIC-UMC-JAN-*)
SCC_FOLDER="${SCC_FOLDERS[$SLURM_ARRAY_TASK_ID]}"

# Check
if [ ! -d "$SCC_FOLDER" ]; then
    echo "SCC folder not found: $SCC_FOLDER"
    exit 1
fi

# Process each gzipped FASTQ file 
fastqc "$SCC_FOLDER"/*.fastq.gz -o "$SCC_FOLDER"
# Check if FastQC ran successfully
if [ $? -ne 0 ]; then
    echo "FastQC failed for folder: $SCC_FOLDER"
    exit 1
fi