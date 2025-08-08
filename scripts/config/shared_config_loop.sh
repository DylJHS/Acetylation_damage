#!/bin/bash

# Project root
PROJ_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage"

# Data
DATA_DIR="$PROJ_DIR/data/SCC-bulkChIC-2"
FASTA_DIR="$DATA_DIR/FASTA_subsets"
REF_DIR="$PROJ_DIR/data/ref_genomes"
TRIMMED_MERG_DIR="$DATA_DIR/trimmed_merged_fastq"
TRIMMED_DIR="$DATA_DIR/trimmed_lane_fastq"

# Results
RES_DIR="$PROJ_DIR/results" # Main results folder

# Temporary files
TEMP_DIR="$RES_DIR/temp"

# BLAST folders
BLAST_DIR="$RES_DIR/BLAST_results"

# Alignment folders
DROS_ALIGN_DIR="$RES_DIR/fly_alignments"
HUMAN_ALIGN_DIR="$RES_DIR/human_alignments"
TAGGED_ALIGN_DIR="$RES_DIR/tagged_alignments"

# Create the folders if they don't exist
mkdir -p "$DROS_ALIGN_DIR" "$HUMAN_ALIGN_DIR" "$TAGGED_ALIGN_DIR" \
    "$FASTA_DIR" "$BLAST_DIR/human" "$BLAST_DIR/drosophila" \
    "$TRIMMED_MERG_DIR" "$TRIMMED_DIR" \
    "$BLAST_DIR" "$TEMP_DIR"

for type in "DROS" "HUMAN" "TAGGED"; do

    # Original alignment folders
    eval "{type}_ALIGN_BOWTIE_DIR=\"${type}_ALIGN_DIR/bowtie2_alignments\""

    # Deduplicated alignment folders
    eval "{type}_DEDUP_DIR=\"${type}_ALIGN_DIR/deduped_alignments\""

    # Deduplicated alignment bam folders
    eval "{type}_DEDUP_BAM_DIR=\"${type}_DEDUP_DIR/bams\""

    # Removed deduplicated alignment bam folders
    eval "{type}_DEDUP_REM_BAMS=\"${type}_DEDUP_BAM_DIR/removed\""

    # Marked deduplicated alignment bam folders
    eval "{type}_DEDUP_MARK_BAMS=\"${type}_DEDUP_BAM_DIR/marked\""

    # Deduplicated alignment stats
    eval "{type}_DEDUP_STATS=\"${type}_DEDUP_DIR/deduped_stats\""

    # Removed deduplicated alignment stats
    eval "{type}_DEDUP_STATS_REM=\"${type}_DEDUP_STATS/removed\""

    # Marked deduplicated alignment stats
    eval "{type}_DEDUP_STATS_MARK=\"${type}_DEDUP_STATS/marked\""

    # Deduplicated alignment bigwigs
    eval "{type}_DEDUP_BIGWIG=\"${type}_DEDUP_DIR/bigwig\""

    # Deduplicated alignment bedgraphs
    eval "{type}_DEDUP_BEDGRAPH=\"${type}_DEDUP_DIR/bedgraph\""

    # Create necessary directories
    eval "mkdir -p \
        \${${type}_ALIGN_BOWTIE_DIR} \${${type}_DEDUP_BAM_DIR} \${${type}_DEDUP_REM_BAMS} \
        \${${type}_DEDUP_MARK_BAMS} \${${type}_DEDUP_STATS_REM} \${${type}_DEDUP_STATS_MARK} \
        \${${type}_DEDUP_BIGWIG} \${${type}_DEDUP_BEDGRAPH}"
done

# Reference genomes
HUMAN_REF="$REF_DIR/human"
HUMAN_GEN_FA="$HUMAN_REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DROS_REF="$REF_DIR/drosophila"
DROS_GEN_FA="$DROS_REF/dmel-all-chromosome-current.fasta.gz"
DROS_INDEX="$DROS_REF/dmel_index"
HUMAN_INDEX="$HUMAN_REF/human_index"

# BLAST read subsample #
SAMPLE_SIZE=5000

# Python script for summary
PY_SCRIPT="$PROJ_DIR/scripts/BLAST/summarise_blast_output.py"