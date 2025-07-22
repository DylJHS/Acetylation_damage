#!/bin/bash
#SBATCH --job-name=bowtie2-build       # Job name
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bowtie2-build-%j.out  # Log file (%j inserts job ID)
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bowtie2-build_%j.err   # Error log
#SBATCH --time=06:00:00                # Max runtime (adjust if needed)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores (adjust as needed)
#SBATCH --mem=32G                       # Memory allocation (adjust as needed)
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh


# Define paths
REFERENCE_DIR="${PROJ_DIR}/data/ref_genomes"
GENOME_FILE="${REFERENCE_DIR}/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"
BUILD_PATH="${REFERENCE_DIR}/BDGP6"
INDEX_PREFIX="${BUILD_PATH}/BDGP6"

# Create output directory if it doesn't exist
mkdir -p "$BUILD_PATH"

# Navigate to the reference genome location
cd $REFERENCE_DIR

# Run Bowtie2 index build
echo "Starting Bowtie2 index build at $(date)"
bowtie2-build --threads 8 $GENOME_FILE $INDEX_PREFIX
echo "Bowtie2 index build completed at $(date)"

