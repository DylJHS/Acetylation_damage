#!/bin/bash
#SBATCH --job-name=run_multiqc
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/run_multiqc-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/run_multiqc-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/wkflw_config.sh

# ARGUMENTS
# if the first argument is "human", use the human alignments directory
if [[ "$1" == "human" ]]; then
    # Define the folder
    FOLDER="$HUMAN_ALIGN_BOWTIE_DIR"  # Change to HUMAN_ALIGN_DIR for human alignments
elif [[ "$1" == "drosophila" ]]; then
    # Define the folder
    FOLDER="$DROS_ALIGN_BOWTIE_DIR"  # Change to DROS_ALIGN_DIR for Drosophila melanogaster alignments
else
    # if nothing is supplied, exit with an error
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi

# Check if the folder exists
if [[ ! -d "$FOLDER" ]]; then
    echo "Error: The specified folder '$FOLDER' does not exist."
    exit 1
fi

# Go to the specified folder
cd "$FOLDER" 

# if second argument is "fastqc", run fastqc
if [[ "$2" == "fastqc" ]]; then
    echo "Running fastqc in directory: $PWD"
    # run the fastqc command
    fastqc -t 6 -o "$FOLDER" *.fq.gz
    echo "FastQC completed in directory: $PWD"
fi

echo "-----------------------------------------------------"

echo "Running MultiQC in directory: $PWD"

# run the multiqc command
multiqc .
echo "MultiQC completed in directory: $PWD"
echo "-----------------------------------------------------"
echo "-----------------------------------------------------"