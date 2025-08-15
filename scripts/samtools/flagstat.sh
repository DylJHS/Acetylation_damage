#!/bin/bash
#SBATCH --job-name=flagstat
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/flagstat-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/flagstat-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/genomics_env_config.sh

echo "---------------------------------------------------------------------------------"
echo "Starting flagstat alignment script for $SLURM_ARRAY_TASK_ID"
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

if [[ "$2" == "-raw" ]]; then
    BAM_FLDR="$RESULTS_DIR"/aligned_bams/bams
    FLAGSTAT_FLDR="$RESULTS_DIR"/aligned_bams/flagstat_results
elif [[ "$2" == "-dedup" && "$3" == "-r" ]]; then
    BAM_FLDR="$RESULTS_DIR"/deduped_alignments/removed/bams
    FLAGSTAT_FLDR="$RESULTS_DIR"/deduped_alignments/removed/flagstat_results
elif [[ "$2" == "-dedup" && "$3" == "-m" ]]; then
    BAM_FLDR="$RESULTS_DIR"/deduped_alignments/marked/bams
    FLAGSTAT_FLDR="$RESULTS_DIR"/deduped_alignments/marked/flagstat_results
else
    echo "Error: Invalid duplication or deduplication type specified."
    exit 1
fi
mkdir -p "$FLAGSTAT_FLDR"

# Define output file
OUTPUT_FILE="$FLAGSTAT_FLDR/${1#-}_flagstat_${2#-}_${3#-}_align_summary.txt"

echo -e "Running ${2#-} ${3#-} flagstat for ${1#-} in folder: \n\t $BAM_FLDR"
echo "---------------------------------------------------------------------------------"

# Define input file
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