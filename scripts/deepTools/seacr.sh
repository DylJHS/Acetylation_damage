#!/bin/bash
#SBATCH --job-name=seacr
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/seacr-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/seacr-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/genomics_env_config.sh

echo "---------------------------------------------------------------------------------"
echo "Starting SEACR script for $SLURM_ARRAY_TASK_ID"
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

# Has to be run on deduped alignments
if [[ "$2" == "-r" ]]; then
    FRAG_BEDS="$RESULTS_DIR"/deduped_alignments/removed/fragments
elif [[ "$2" == "-m" ]]; then
    FRAG_BEDS="$RESULTS_DIR"/deduped_alignments/marked/fragments
else
    echo "Error: Invalid deduplication type specified."
    exit 1
fi

# Create the output folder
SEACR_RES_DIR="$FRAG_BEDS/../seacr_results"
mkdir -p "$SEACR_RES_DIR"

# Set the fragment for the current job
fragments=("${FRAG_BEDS}"/*.fragments.bedgraph)
fragment="${fragments[$SLURM_ARRAY_TASK_ID]}"

output_prefix="${SEACR_RES_DIR}/(basename "$fragment" .fragments.bedgraph)_"


# Run the SEACR analysis
echo "Runnign the SEACR analysis using frament: $fragment"
echo "---------------------------------------------------------------------------------"
SEACR_1.3.sh $fragment 0.05 non stringent output_prefix
eccho "Finished running the SEACR analysis" 
echo "---------------------------------------------------------------------------------"