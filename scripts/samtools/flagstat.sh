#!/bin/bash
#SBATCH --job-name=flagstat
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/flagstat-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/flagstat-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/wkflw_config.sh
echo "---------------------------------------------------------------------------------"
echo "Starting flagstat alignment script for $SLURM_ARRAY_TASK_ID"
echo "---------------------------------------------------------------------------------"

# Define paths
if [[ "$1" == "-drosophila" ]]; then
    RESULTS_DIR="$DROS_ALIGN_DIR"
    if [[ "$2" == "-bowtie" ]]; then
        BAM_FLDR="$DROS_ALIGN_BOWTIE_DIR"
    elif [[ "$2" == "-dedup" ]]; then
        BAM_FLDR="$DROS_DEDUP_DIR"
    fi
elif [[ "$1" == "-human" ]]; then
    RESULTS_DIR="$HUMAN_ALIGN_DIR"
    if [[ "$2" == "-bowtie" ]]; then
        BAM_FLDR="$HUMAN_ALIGN_BOWTIE_DIR"
    elif [[ "$2" == "-dedup" ]]; then
        BAM_FLDR="$HUMAN_DEDUP_DIR"
    fi
fi

if [[ "$2" == "-dedup" ]]; then
    OUTPUT_FILE="$RESULTS_DIR/flagstat_dedup_align_summary.txt"
elif [[ "$2" == "-bowtie" ]]; then
    OUTPUT_FILE="$RESULTS_DIR/flagstat_bowtie_align_summary.txt"
fi

echo "Running $2 for species: $1"
echo "---------------------------------------------------------------------------------"

# Define output file

BAM_FILES=("${BAM_FLDR}"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Define the lock function 
write_with_lock() {
    local msg="$1"
    local tmp_lock="${TEMP_DIR}/lock_${SLURM_JOB_ID}.lock"

    # Wait until we can create the lock
    while ! mkdir "$tmp_lock" 2>/dev/null; do
        sleep 0.1
    done

    # Critical section: write safely
    echo "$msg" >> "$OUTPUT_FILE"

    # Release lock
    rmdir "$tmp_lock"
}

# if first array job, clear the output file and create the file if it doesn't exist
if [[ $SLURM_ARRAY_TASK_ID -eq 0 ]]; then
    if [[ -f $OUTPUT_FILE ]]; then
        rm "$OUTPUT_FILE"
        echo "Cleared existing output file: $OUTPUT_FILE"
        echo "---------------------------------------------------------------------------------"
    fi
    touch "$OUTPUT_FILE"
    echo "Created new output file: $OUTPUT_FILE"
    echo "---------------------------------------------------------------------------------"
fi  

# Loop through BAM files and save stats
if [[ -f $BAM_FILE ]]; then
    name=$(basename "$BAM_FILE")
    echo "Processing: $name"
    intro="Processing BAM file: $name"
    stats=$(samtools flagstat "$BAM_FILE")
    final_output="$intro"$'\n'"$stats"$'\n'"-------------------------------------------------------"$'\n'
    write_with_lock "$final_output"  # Write the stats to the output file
    # Save the stats to a separate file
    echo "$stats" > "${BAM_FLDR}/flagstat_${name}.txt"
    echo "Alignment stats for $name saved."
else
    echo "Warning: No BAM files found in $RESULTS_DIR" 
    write_with_lock "Warning: No BAM files found in $RESULTS_DIR"
fi
echo "---------------------------------------------------------------------------------"
echo "Flagstat alignment script for $SLURM_ARRAY_TASK_ID completed."
echo "---------------------------------------------------------------------------------"