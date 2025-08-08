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

for align_type in "DROS" "HUMAN" "TAGGED"; do

    # Base folders
    eval "${align_type}_ALIGN_BOWTIE_DIR=\"\${${align_type}_ALIGN_DIR}/aligned_bams\""
    eval "${align_type}_DEDUP_DIR=\"\${${align_type}_ALIGN_DIR}/deduped_alignments\""
    eval "${align_type}_DEDUP_BAM_DIR=\"\${${align_type}_DEDUP_DIR}/bams\""
    eval "${align_type}_DEDUP_STATS=\"\${${align_type}_DEDUP_DIR}/deduped_stats\""
    eval "${align_type}_DEDUP_BIGWIG=\"\${${align_type}_DEDUP_DIR}/bigwig\""
    eval "${align_type}_DEDUP_BEDGRAPH=\"\${${align_type}_DEDUP_DIR}/bedgraph\""

    for dedup_type in "REM" "MARK"; do
        if [[ "$dedup_type" == "REM" ]]; then
            tag="removed"
        else
            tag="marked"
        fi

        # Define dedup subpaths
        eval "${align_type}_DEDUP_${dedup_type}_BAMS=\"\${${align_type}_DEDUP_BAM_DIR}/${tag}\""
        eval "${align_type}_DEDUP_STATS_${dedup_type}=\"\${${align_type}_DEDUP_STATS}/${tag}\""
        eval "${align_type}_DEDUP_BIGWIG_${dedup_type}=\"\${${align_type}_DEDUP_BIGWIG}/${tag}\""
        eval "${align_type}_DEDUP_BEDGRAPH_${dedup_type}=\"\${${align_type}_DEDUP_BEDGRAPH}/${tag}\""

        # Create the folders
        eval "mkdir -p \
            \${${align_type}_DEDUP_${dedup_type}_BAMS} \
            \${${align_type}_DEDUP_STATS_${dedup_type}} \
            \${${align_type}_DEDUP_BIGWIG_${dedup_type}} \
            \${${align_type}_DEDUP_BEDGRAPH_${dedup_type}}"
    done
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