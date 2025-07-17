#!/bin/bash
#SBATCH --job-name=run_blastn
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/blastn_query-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/blastn_query-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-23
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env
export LD_LIBRARY_PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env/lib:$LD_LIBRARY_PATH"

# Directory where FASTA subset files are saved
DATA_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/data/SCC-bulkChIC-2"

FASTA_DIR="$DATA_DIR/FASTA_subsets"
# Create FASTA directory if it doesn't exist
mkdir -p "$FASTA_DIR"

# Directory where result files will be saved
RES_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results/BLAST_results"
mkdir -p "$RES_DIR/human" "$RES_DIR/drosophila"

# Reference FASTA files for the three species
DB_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/data/ref_genomes"
HUMAN_REF="$DB_DIR/Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa" # Human reference
DROS_REF="$DB_DIR/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"  # Drosophila reference

# Get the FASTA file for this array task
FASTA_FILES=($FASTA_DIR/*_subset.fasta)
FASTA_FILE="${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}"

# Sanity check
if [ ! -f "$FASTA_FILE" ]; then
    echo "FASTA file not found: $FASTA_FILE"
    exit 1
fi

# Loop through each subset FASTA file in the FASTA directory and run BLAST against each database.
base=$(basename "$FASTA_FILE" _subset.fasta)

for species in human drosophila; do
    # Select the appropriate database for each species.
    if [ "$species" = "human" ]; then
        db_name="${HUMAN_REF%.*}_db"
    elif [ "$species" = "drosophila" ]; then
        db_name="${DROS_REF%.*}_db"
    fi
    output_file="${RES_DIR}/${species}/${base}_${species}_blast_results.txt"
    echo "Running BLAST for $base against $species database..."
    blastn -query "$FASTA_FILE" \
        -task megablast \
        -perc_identity 95 \
        -max_target_seqs 1 \
        -max_hsps 1 \
        -ungapped \
        -db "$db_name" \
        -num_threads 4 \
        -out "$output_file" \
        -evalue 1e-20 \
        -outfmt 6
done
