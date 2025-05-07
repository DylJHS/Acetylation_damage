#!/bin/bash
#SBATCH --job-name=BLAST_query
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/blast_query-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/blast_query-%j.err
#SBATCH --time=12:00:00
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


# Directory containing reference genomes 
DB_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage/reference_genomes"
# Directory where FASTA subset files are saved
FASTA_DIR="/hpc/shared/onco_janssen/dhaynessimmons/data/fly_acetylation_damage/FASTA_subsets"
# Directory where result files will be saved
RES_DIR="/hpc/shared/onco_janssen/dhaynessimmons/results/fly_acetylation"

# Reference FASTA files for the three species
HUMAN_REF="$DB_DIR/Homo_sapiens.GRCh38.dna.toplevel.fa" # Human reference
DROS_REF="$DB_DIR/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"  # Drosophila reference

# Create BLAST databases if they do not already exist.
# The script checks for the presence of a .nhr file (for nucleotide databases).
for ref in "$HUMAN_REF" "$DROS_REF"; do
    db_name="${ref%.*}_db"  # e.g., human_db, BDGP6_db
    if [ ! -f "${db_name}.nhr" ]; then
        echo "Creating BLAST database for $ref"
        makeblastdb -in "$ref" -dbtype nucl -out "$db_name"
    else
        echo "BLAST database for $ref already exists."
    fi
done

# Loop through each subset FASTA file in the FASTA directory and run BLAST against each database.
for subset in "$FASTA_DIR"/*_subset.fasta; do
    base=$(basename "$subset" _subset.fasta)
    for species in human drosophila; do
        # Select the appropriate database for each species.
        if [ "$species" = "human" ]; then
            db_name="${HUMAN_REF%.*}_db"
        elif [ "$species" = "drosophila" ]; then
            db_name="${DROS_REF%.*}_db"
        fi
        output_file="${RES_DIR}/${base}_${species}_blast_results.txt"
        echo "Running BLAST for $subset against $species database..."
        blastn -query "$subset" -db "$db_name" -num_threads 4 -out "$output_file" -evalue 1e-5 -outfmt 6
    done
done