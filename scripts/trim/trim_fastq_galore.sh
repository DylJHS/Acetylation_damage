#!/bin/bash
#SBATCH --job-name=trim_fastq
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/trim_fastq-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/trim_fastq-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

# Identify SCC folder for this array task
SCC_FOLDERS=("$DATA_DIR"/bulkChIC-UMC-JAN-*)
SCC_FOLDER="${SCC_FOLDERS[$SLURM_ARRAY_TASK_ID]}"

# Output check
if [ ! -d "$SCC_FOLDER" ]; then
    echo "SCC folder not found: $SCC_FOLDER"
    exit 1
fi

echo "Processing folder: $SCC_FOLDER"
echo " -------------------------------------------------------------------"

# Find all R1 files and match them with R2
for R1 in "$SCC_FOLDER"/*_R1_*.fastq.gz; do
    [ -e "$R1" ] || continue
    R2="${R1/_R1_/_R2_}"

    echo "Processing R1: $R1"
    echo "Corresponding R2: $R2"
    
    if [ ! -f "$R2" ]; then
        echo "  No matching R2 for $R1"
        continue
    fi

    echo "------------------------------------------------------------------------"

    # Extract sample ID
    base_name=$(basename "$R1" | sed -E 's/_R1_001\.fastq\.gz$//')

    echo "  Trimming: $base_name"
    
    trim_galore \
        --illumina \
        --output_dir "$TRIMMED_DIR" \
        --no_report_file \
        --paired $R1 $R2
done

echo "Trimming complete for folder: $SCC_FOLDER"
