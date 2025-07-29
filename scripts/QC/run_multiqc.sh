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
if [[ "$1" == "-human" ]]; then
    # Define the folder
    if [[ "$2" == "-bowtie" ]]; then
        FOLDER="$HUMAN_ALIGN_BOWTIE_DIR"  # Change to HUMAN_ALIGN_DIR for human alignments
    elif [[ "$2" == "-dedup" ]]; then
        FOLDER="$HUMAN_DEDUP_DIR"  # Change to HUMAN_DEDUP_DIR for deduplicated human alignments
    fi
    FOLDER="$HUMAN_ALIGN_BOWTIE_DIR"  # Change to HUMAN_ALIGN_DIR for human alignments
elif [[ "$1" == "-drosophila" ]]; then
    # Define the folder
    if [[ "$2" == "-bowtie" ]]; then
        FOLDER="$DROS_ALIGN_BOWTIE_DIR"  # Change to DROS_ALIGN_DIR for drosophila alignments
    elif [[ "$2" == "-dedup" ]]; then
        FOLDER="$DROS_DEDUP_DIR"  # Change to DROS_DEDUP_DIR for deduplicated drosophila alignments
    fi
else
    echo "Usage: $0 -human|-drosophila [-bowtie|-dedup]"
    echo "Please specify either -human or -drosophila followed by -bowtie or -dedup if needed."
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
if [[ "$3" == "-fastqc" ]]; then
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