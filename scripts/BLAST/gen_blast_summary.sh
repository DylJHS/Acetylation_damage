#!/bin/bash
#SBATCH --job-name=gen_blast_summary
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/gen_blast_summary-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/gen_blast_summary-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env
export LD_LIBRARY_PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env/lib:$LD_LIBRARY_PATH"


# Directory where result files will be saved
RES_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/BLAST_results"

# Path to the Python summary script
PY_SCRIPT="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/summarise_blast_output.py"

# Number of reads to sample from each file
SAMPLE_SIZE=1000

# Summarise the BLAST results using the Python script
for folder in "$RES_DIR"/*; do 
    folder_base="$(basename "$folder")"
    echo "Processing folder: $folder_base"
    BLAST_PATH="$RES_DIR/$folder_base"
    OUTPUT_SUMMARY="$RES_DIR/${folder_base}_summary_blast_results.txt"

    # Empty or create the master summary file
    > "$OUTPUT_SUMMARY"

    # Loop over each BLAST result file ending with *_blast_results.txt in the RES_DIR
    for file in "$BLAST_PATH"/JAN*_blast_results.txt; do
        # Extract the base filename (e.g., SCC-bulkChIC-UMC-JAN-003_R1_drosophila_blast_results.txt)
        base=$(basename "$file")
        echo "Processing file $base ..." | tee -a "$OUTPUT_SUMMARY"
        # Call the Python script, passing the file name (the Python script expects the file to be in RES_DIR)
        python3 "$PY_SCRIPT" "$file" "$SAMPLE_SIZE" >> "$OUTPUT_SUMMARY"
        echo "----------------------------------------" >> "$OUTPUT_SUMMARY"
    done

    echo "All summaries written to $OUTPUT_SUMMARY"
done 