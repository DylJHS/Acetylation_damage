#!/bin/bash
#SBATCH --job-name=fastq_concat
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fastq_concat-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/fastq_concat-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/BLAST/BLAST_wkflw_config.sh

echo "      Concatenating trimmed lanes for each sample..."

# Find all unique base sample names (without lane and read info)
trimmed_files=("$TRIMMED_DIR"/*_R1_001_val_1.fq.gz)
echo "trimmed_files: ${trimmed_files[@]}"
sample_prefixes=()

for f in "${trimmed_files[@]}"; do
    name=$(basename "$f" | sed -E 's/_S[0-9]+_L00[0-9]_R1_001_val_1\.fq\.gz//')
    sample_prefixes+=("$name")
done

# Deduplicate sample names
unique_samples=($(printf "%s\n" "${sample_prefixes[@]}" | sort -u))
echo "Found ${#unique_samples[@]} unique samples to process."

echo "----------------------------------------"
for prefix in "${unique_samples[@]}"; do
    echo "Found sample prefix: $prefix"
done
echo "----------------------------------------"

for sample in "${unique_samples[@]}"; do
    echo "Processing sample: $sample"
    # Define patterns to match all lanes for this sample
    R1_files=("$TRIMMED_DIR"/${sample}_S*_L00*_R1_001_val_1.fq.gz)
    R2_files=("$TRIMMED_DIR"/${sample}_S*_L00*_R2_001_val_2.fq.gz)

    echo "----------------------------------------"
    echo "R1 files for sample $sample:"
    for r1 in "${R1_files[@]}"; do
        echo "  Found R1 file: $r1"
    done

    echo "R2 files for sample $sample:"
    for r2 in "${R2_files[@]}"; do
        echo "  Found R2 file: $r2"
    done

    echo "----------------------------------------"
    # Output filenames
    merged_R1="$TRIMMED_MERG_DIR/${sample}_merged_R1.fq.gz"
    merged_R2="$TRIMMED_MERG_DIR/${sample}_merged_R2.fq.gz"

    echo "  Merging lanes for sample $sample to:"
    echo "    Merged R1: $merged_R1"
    echo "    Merged R2: $merged_R2"
    cat "${R1_files[@]}" > "$merged_R1"
    cat "${R2_files[@]}" > "$merged_R2"
    echo "----------------------------------------"
done

echo "          Concatenation complete."