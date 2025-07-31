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
elif [[ "$1" == "-tagged" ]]; then
    RESULTS_DIR="$TAGGED_BAM_ORI_DIR"
    if [[ "$2" == "-bowtie" ]]; then
        BAM_FLDR="$RESULTS_DIR"
    elif [[ "$2" == "-dedup" ]]; then
        BAM_FLDR="$TAGGED_DEDUP_DIR"
    fi
else
    echo "Error: Invalid species or alignment type specified."
    exit 1
fi

FLAGSTAT_FLDR="$BAM_FLDR/flagstat_results"
OUTPUT_FILE="$BAM_FLDR/${1#-}_flagstat_${2#-}_align_summary.txt"

mkdir -p "$FLAGSTAT_FLDR"

echo "Running ${2#-} for : ${1#-}"
echo "---------------------------------------------------------------------------------"

# Define output file

BAM_FILES=("${BAM_FLDR}"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
name=$(basename "$BAM_FILE")

# Define the lock function 
write_with_lock() {
    local msg="$1"
    local tmp_lock="${TEMP_DIR}/lock_${SLURM_JOB_ID}.lock"
    local wait_time=0
    local max_wait=60  # seconds

    # Wait until we can create the lock
    while ! mkdir "$tmp_lock" 2>/dev/null; do
        sleep 1
        ((wait_time++))
        if [[ $wait_time -ge $max_wait ]]; then
            echo "Lock wait timeout exceeded for task $SLURM_ARRAY_TASK_ID" >&2
            return 1
        fi
    done

    # Critical section: write safely
    echo "$msg" >> "$OUTPUT_FILE"

    # Release lock
    rmdir "$tmp_lock"
}

# Check if the output file exists
if [[ -e "$OUTPUT_FILE" ]]; then
    # check if the file contains the name of the BAM file
    if grep -q "$name" "$OUTPUT_FILE"; then
        # clear the output file
        echo "Output file already contains stats for $name. Clearing the file."
        echo "---------------------------------------------------------------------------------"
        > "$OUTPUT_FILE"
    else
        echo "Output file exists but does not contain stats for $name. Appending to the file."
        echo "---------------------------------------------------------------------------------"
    fi
else
    echo "Output file does not exist. Creating a new file."
    echo "---------------------------------------------------------------------------------"
    touch "$OUTPUT_FILE"    
fi

# Loop through BAM files and save stats
if [[ -e $BAM_FILE ]]; then
    echo "Processing: $name"
    echo "---------------------------------------------------------------------------------"
    intro="BAM file: $name"
    stats=$(samtools flagstat "$BAM_FILE")
    final_output="$intro"$'\n'"$stats"$'\n'"-------------------------------------------------------"$'\n'
    echo "$stats" > "${FLAGSTAT_FLDR}/flagstat_${name}.txt"
    # sleep random delay
    sleep $((RANDOM % 30 + 1))  # Sleep for a random time between 1 and 30 seconds
    write_with_lock "$final_output"  # Write the stats to the output file
    # Save the stats to a separate file
    echo "Alignment stats for $name saved."
else
    echo "Warning: No BAM files found in $RESULTS_DIR" 
    write_with_lock "Warning: No BAM files found in $RESULTS_DIR"
fi
echo "---------------------------------------------------------------------------------"
echo "Flagstat alignment script for $SLURM_ARRAY_TASK_ID completed."
echo "---------------------------------------------------------------------------------"