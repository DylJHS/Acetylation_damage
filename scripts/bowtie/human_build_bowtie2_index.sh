#!/bin/bash
#SBATCH --job-name=bowtie2-build       # Job name
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/bowtie2-build-%j.out  # Log file (%j inserts job ID)
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/bowtie2-build_%j.err   # Error log
#SBATCH --time=06:00:00                # Max runtime (adjust if needed)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPU cores (adjust as needed)
#SBATCH --mem=32G                       # Memory allocation (adjust as needed)
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
# Manually initialize Conda in a batch job
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env

# Define paths
REFERENCE_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage/reference_genomes"
GENOME_FILE="${REFERENCE_DIR}/Homo_sapiens.GRCh38.dna.toplevel.fa"
INDEX_PREFIX="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage/GRCh38/GRCh38"


# Navigate to the reference genome location
cd $REFERENCE_DIR

# Run Bowtie2 index build
echo "Starting Bowtie2 index build at $(date)"
bowtie2-build --threads 8 $GENOME_FILE $INDEX_PREFIX
echo "Bowtie2 index build completed at $(date)"

