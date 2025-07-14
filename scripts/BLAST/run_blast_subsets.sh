#!/bin/bash
#SBATCH --job-name=run_blast_subsets
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/blast_query-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/blast_query-%j.err
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


# Directory containing reference genomes 
DB_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/data/ref_genomes"
# Directory where FASTA subset files are saved
DATA_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/data/SCC-bulkChIC-2"

FASTA_DIR="$DATA_DIR/FASTA_subsets"
# Create FASTA directory if it doesn't exist
mkdir -p "$FASTA_DIR"
# Directory where result files will be saved
RES_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/results"
mkdir -p "$RES_DIR"

# Number of reads to sample from each file
SAMPLE_SIZE=1000

# Loop through each SCC folder (assumed to have names like SCC-bulkChIC-UMC-JAN-*)
for folder in "$DATA_DIR"/bulkChIC-UMC-JAN-*; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"

        # Process each gzipped FASTQ file for R1 and R2 in the folder
        for file in "$folder"/*_R[12]_*.fastq.gz; do
            # Check if file exists (in case pattern doesn't match any files)
            [ -e "$file" ] || continue
            echo "  Processing file: $file"
            # Get the folder's base name (e.g., SCC-bulkChIC-UMC-JAN-003)
            folder_id=$(basename "$folder")

            # Extract FASTQ short ID
            fastq_base=$(basename "$file" .fastq.gz)
            short_id=$(echo "$fastq_base" | sed -E 's/^.*(JAN-[0-9]+)_.*(L[0-9]+)_R([12])_.*/\1_\2_R\3/')
            
            OUTPUT_FASTA="$FASTA_DIR/${short_id}_subset.fasta"
            echo "    Output FASTA: $OUTPUT_FASTA"

            # Decompress the file, sample reads, and convert to FASTA in one pipeline
            zcat "$file" | seqtk sample -s100 - "$SAMPLE_SIZE" | seqtk seq -a > "$OUTPUT_FASTA"
        done
    fi
done

# Reference FASTA files for the three species
HUMAN_REF="$DB_DIR/Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa" # Human reference
DROS_REF="$DB_DIR/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"  # Drosophila reference

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