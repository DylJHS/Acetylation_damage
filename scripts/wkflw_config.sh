#!/bin/bash
# This is the workflow configuration script for the Fly Acetylation Damage project.

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
REF_DIR="$PROJ_DIR/data/ref_genomes"
TRIMMED_MERG_DIR="$DATA_DIR/trimmed_merged_fastq"
TRIMMED_DIR="$DATA_DIR/trimmed_lane_fastq"
TAGGED_BAM_ORI_DIR="$DATA_DIR/tagged_bam_files"

# Results
RES_DIR="$PROJ_DIR/results"

BLAST_DIR="$RES_DIR/BLAST_results"

DROS_ALIGN_DIR="$RES_DIR/fly_alignments"
HUMAN_ALIGN_DIR="$RES_DIR/human_alignments"
TAGGED_ALIGN_DIR="$RES_DIR/tagged_alignments"

DROS_ALIGN_BOWTIE_DIR="$DROS_ALIGN_DIR/bowtie2_alignments"
HUMAN_ALIGN_BOWTIE_DIR="$HUMAN_ALIGN_DIR/bowtie2_alignments"

DROS_DEDUP_DIR="$DROS_ALIGN_DIR/deduped_alignments"
HUMAN_DEDUP_DIR="$HUMAN_ALIGN_DIR/deduped_alignments"
TAGGED_DEDUP_DIR="$TAGGED_ALIGN_DIR/deduped_alignments"

DROS_DEDUP_STATS="$DROS_DEDUP_DIR/deduped_stats"
HUMAN_DEDUP_STATS="$HUMAN_DEDUP_DIR/deduped_stats"
TAGGED_DEDUP_STATS="$TAGGED_DEDUP_DIR/deduped_stats"

TEMP_DIR="$RES_DIR/temp"

# Reference genomes
HUMAN_REF="$REF_DIR/Homo_sapiens.GRCh38.dna.toplevel.fa/Homo_sapiens.GRCh38.dna.toplevel.fa"
DROS_REF="$REF_DIR/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa"
DROS_INDEX="$REF_DIR/BDGP6"
HUMAN_INDEX="$REF_DIR/GRCh38"

# Create necessary directories
mkdir -p "$FASTA_DIR" "$BLAST_DIR/human" "$BLAST_DIR/drosophila" \
    "$TRIMMED_MERG_DIR" "$TRIMMED_DIR" "$DROS_ALIGN_DIR" \
    "$HUMAN_ALIGN_DIR" "$TEMP_DIR" "$DROS_DEDUP_DIR" \
    "$HUMAN_DEDUP_DIR" "$DROS_DEDUP_STATS" "$HUMAN_DEDUP_STATS" \
    "$TAGGED_ALIGN_DIR" "$TAGGED_DEDUP_DIR" "$TAGGED_DEDUP_STATS"

# BLAST read subsample #
SAMPLE_SIZE=5000

# Python script for summary
PY_SCRIPT="$PROJ_DIR/scripts/BLAST/summarise_blast_output.py"