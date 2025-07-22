#!/bin/bash

# Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env
export LD_LIBRARY_PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env/lib:$LD_LIBRARY_PATH"

# Project root
PROJ_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage"

# Data
DATA_DIR="$PROJ_DIR/data/SCC-bulkChIC-2"
FASTA_DIR="$DATA_DIR/FASTA_subsets"
DB_DIR="$PROJ_DIR/data/ref_genomes"
TRIMMED_DIR="$DATA_DIR/trimmed_fastq"

# Results
RES_DIR="$PROJ_DIR/results/BLAST_results"
mkdir -p "$FASTA_DIR" "$RES_DIR/human" "$RES_DIR/drosophila"

# Reference genomes
HUMAN_REF="$DB_DIR/Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa"
DROS_REF="$DB_DIR/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"

# Sample size
SAMPLE_SIZE=5000

# Python script for summary
PY_SCRIPT="$PROJ_DIR/scripts/BLAST/summarise_blast_output.py"
