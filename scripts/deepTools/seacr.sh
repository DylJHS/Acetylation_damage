#!/bin/bash
#SBATCH --job-name=seacr
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/seacr/seacr/seacr-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/seacr/seacr-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5  # Adjust based on the number of BAM files
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

echo "---------------------------------------------------------------------------------"
echo "Starting SEACR script for $SLURM_ARRAY_TASK_ID"
echo "---------------------------------------------------------------------------------"

# Define paths
FRAG_BEDS="/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/${1}"

echo "running the analysis on the files in folder: $FRAG_BEDS"
echo "---------------------------------------------------------------------------------"

# Create the output folder
SEACR_RES_DIR="$FRAG_BEDS/../seacr_results"
mkdir -p "$SEACR_RES_DIR"

# Set the fragment for the current job
fragments=("${FRAG_BEDS}"/*.fragments.bedgraph)
if [[ ${#fragments[@]} -eq 0 ]]; then
    echo "No .fragments.bedgraph files found in: $FRAG_BEDS"
    exit 1
else
    fragment="${fragments[$SLURM_ARRAY_TASK_ID]}"
fi

output_prefix="${SEACR_RES_DIR}/$(basename "$fragment" .fragments.bedgraph)_"


# Run the SEACR analysis with settigns of non relaxed and 0.05 FDR
echo "Running the SEACR analysis using fragment: $(basename "$fragment" .fragments.bedgraph)"
echo "---------------------------------------------------------------------------------"
SEACR_1.3.sh $fragment 0.05 non relaxed $output_prefix
echo "---------------------------------------------------------------------------------"
echo "Finished running the SEACR analysis" 
echo "---------------------------------------------------------------------------------" 

