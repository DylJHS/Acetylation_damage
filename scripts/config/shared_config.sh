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

# Original alignment folders
DROS_ALIGNMENT_DIR="$DROS_ALIGN_DIR/aligned_bams"
HUMAN_ALIGNMENT_DIR="$HUMAN_ALIGN_DIR/aligned_bams"
TAGGED_ALIGNMENT_DIR="$TAGGED_ALIGN_DIR/aligned_bams"

# Deduplicated alignment folders
DROS_DEDUP_DIR="$DROS_ALIGN_DIR/deduped_alignments"
HUMAN_DEDUP_DIR="$HUMAN_ALIGN_DIR/deduped_alignments"
TAGGED_DEDUP_DIR="$TAGGED_ALIGN_DIR/deduped_alignments"

# Deduplicated alignment bam folders
DROS_DEDUP_BAM_DIR="$DROS_DEDUP_DIR/bams"
HUMAN_DEDUP_BAM_DIR="$HUMAN_DEDUP_DIR/bams"
TAGGED_DEDUP_BAM_DIR="$TAGGED_DEDUP_DIR/bams"

# Removed deduplicated alignment bam folders
DROS_DEDUP_REM_BAMS="$DROS_DEDUP_BAM_DIR/removed"
HUMAN_DEDUP_REM_BAMS="$HUMAN_DEDUP_BAM_DIR/removed"
TAGGED_DEDUP_REM_BAMS="$TAGGED_DEDUP_BAM_DIR/removed"

# Marked deduplicated alignment bam folders
DROS_DEDUP_MARK_BAMS="$DROS_DEDUP_BAM_DIR/marked"
HUMAN_DEDUP_MARK_BAMS="$HUMAN_DEDUP_BAM_DIR/marked"
TAGGED_DEDUP_MARK_BAMS="$TAGGED_DEDUP_BAM_DIR/marked"

# Deduplicated alignment stats
DROS_DEDUP_STATS="$DROS_DEDUP_DIR/deduped_stats"
HUMAN_DEDUP_STATS="$HUMAN_DEDUP_DIR/deduped_stats"
TAGGED_DEDUP_STATS="$TAGGED_DEDUP_DIR/deduped_stats"

# Removed deduplicated alignment stats
DROS_DEDUP_STATS_REM="$DROS_DEDUP_STATS/removed"
HUMAN_DEDUP_STATS_REM="$HUMAN_DEDUP_STATS/removed"
TAGGED_DEDUP_STATS_REM="$TAGGED_DEDUP_STATS/removed"

# Marked deduplicated alignment stats
DROS_DEDUP_STATS_MARK="$DROS_DEDUP_STATS/marked"
HUMAN_DEDUP_STATS_MARK="$HUMAN_DEDUP_STATS/marked"
TAGGED_DEDUP_STATS_MARK="$TAGGED_DEDUP_STATS/marked"

# Deduplicated alignment bigwigs
DROS_DEDUP_BIGWIG="$DROS_DEDUP_DIR/bigwig"
HUMAN_DEDUP_BIGWIG="$HUMAN_DEDUP_DIR/bigwig"
TAGGED_DEDUP_BIGWIG="$TAGGED_DEDUP_DIR/bigwig"

# Deduplicated alignment bedgraphs
DROS_DEDUP_BEDGRAPH="$DROS_DEDUP_DIR/bedgraph"
HUMAN_DEDUP_BEDGRAPH="$HUMAN_DEDUP_DIR/bedgraph"
TAGGED_DEDUP_BEDGRAPH="$TAGGED_DEDUP_DIR/bedgraph"

# Reference genomes
HUMAN_REF="$REF_DIR/human"
HUMAN_GEN_FA="$HUMAN_REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DROS_REF="$REF_DIR/drosophila"
DROS_GEN_FA="$DROS_REF/dmel-all-chromosome-current.fasta.gz"
DROS_INDEX="$DROS_REF/dmel_index"
HUMAN_INDEX="$HUMAN_REF/human_index"

# Create necessary directories
mkdir -p "$FASTA_DIR" "$BLAST_DIR/human" "$BLAST_DIR/drosophila" \
    "$TRIMMED_MERG_DIR" "$TRIMMED_DIR" "$DROS_ALIGN_DIR" \
    "$HUMAN_ALIGN_DIR" "$TEMP_DIR" "$DROS_DEDUP_DIR" \
    "$HUMAN_DEDUP_DIR" "$DROS_DEDUP_STATS" "$HUMAN_DEDUP_STATS" \
    "$TAGGED_ALIGN_DIR" "$TAGGED_DEDUP_DIR" "$TAGGED_DEDUP_STATS" \
    "$DROS_DEDUP_BIGWIG" "$HUMAN_DEDUP_BIGWIG" "$TAGGED_DEDUP_BIGWIG" \
    "$DROS_DEDUP_BEDGRAPH" "$HUMAN_DEDUP_BEDGRAPH" "$TAGGED_DEDUP_BEDGRAPH" \
    "$DROS_DEDUP_BAM_DIR" "$HUMAN_DEDUP_BAM_DIR" "$TAGGED_DEDUP_BAM_DIR" \
    "$DROS_DEDUP_REM_BAMS" "$HUMAN_DEDUP_REM_BAMS" "$TAGGED_DEDUP_REM_BAMS" \
    "$DROS_DEDUP_STATS_REM" "$HUMAN_DEDUP_STATS_REM" "$TAGGED_DEDUP_STATS_REM" \
    "$DROS_DEDUP_MARK_BAMS" "$HUMAN_DEDUP_MARK_BAMS" "$TAGGED_DEDUP_MARK_BAMS" \
    "$DROS_DEDUP_STATS_MARK" "$HUMAN_DEDUP_STATS_MARK" "$TAGGED_DEDUP_STATS_MARK"

# BLAST read subsample #
SAMPLE_SIZE=5000

# Python script for summary
PY_SCRIPT="$PROJ_DIR/scripts/BLAST/summarise_blast_output.py"