#!/bin/bash
#SBATCH --job-name=trim_fastq
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/trim_fastq-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/trim_fastq-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
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

# Define the TRIMMED_FLDR
TRIMMED_FLDR="$SCC_FOLDER/trimmed_lane_fastq"
# Clear the directory if it already exists
if [ "$(ls -A "$TRIMMED_FLDR")" ]; then
    echo "Clearing existing files in $TRIMMED_FLDR"
    rm -rf "$TRIMMED_FLDR"/*
else
    echo "Creating new directory: $TRIMMED_FLDR"
    mkdir -p "$TRIMMED_FLDR"
fi

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
        --output_dir "$TRIMMED_FLDR" \
        --no_report_file \
        --paired $R1 $R2
done


# Concatenate trimmed reads across lanes into single files for each sample

echo "Concatenating trimmed lanes for each sample..."

# Find all unique base sample names (without lane and read info)
trimmed_files=("$TRIMMED_FLDR"/*_R1_001_val_1.fq.gz)
sample_prefixes=()

for f in "${trimmed_files[@]}"; do
    name=$(basename "$f" | sed -E 's/_S[0-9]+_L00[0-9]_R1_001_val_1\.fq\.gz//')
    sample_prefixes+=("$name")
done

# Deduplicate sample names
unique_samples=($(printf "%s\n" "${sample_prefixes[@]}" | sort -u))

for sample in "${unique_samples[@]}"; do
    # Define patterns to match all lanes for this sample
    R1_files=("$TRIMMED_FLDR"/${sample}_S*_L00*_R1_001_val_1.fq.gz)
    R2_files=("$TRIMMED_FLDR"/${sample}_S*_L00*_R2_001_val_2.fq.gz)

    # Output filenames
    merged_R1="$TRIMMED_FLDR/${sample}_merged_R1.fq.gz"
    merged_R2="$TRIMMED_FLDR/${sample}_merged_R2.fq.gz"

    echo "  Merging lanes for sample: $sample"
    cat "${R1_files[@]}" > "$merged_R1"
    cat "${R2_files[@]}" > "$merged_R2"

    # Save the merged files to the TRIMMED_DIR
    echo "  Moving merged files to $TRIMMED_DIR"
    mv "$merged_R1" "$TRIMMED_DIR/${sample}_R1.fastq.gz"
    mv "$merged_R2" "$TRIMMED_DIR/${sample}_R2.fastq.gz"
done

echo "Concatenation complete."