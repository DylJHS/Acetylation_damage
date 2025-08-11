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

for alignment_type in "DROS" "HUMAN" "TAGGED"; do
    echo "creating the alignment type folder"
    if [[ "$alignment_type" == "DROS" ]]; then
        alignment_type_tag="drosophila"
    elif [[ "$alignment_type" == "HUMAN" ]]; then
        alignment_type_tag="human"
    else
        alignment_type_tag="tagged"
    fi

    # Alignment folders
    eval "${alignment_type}_ALIGN_DIR=\"\${RES_DIR}/${alignment_type_tag}_alignments\""

    # Original alignment folders
    eval "${alignment_type}_ALIGNMENT_DIR=\"\${${alignment_type}_ALIGN_DIR}/aligned_bams\""
    eval "${alignment_type}_ALIGN_DIR_BAMS=\"\${${alignment_type}_ALIGNMENT_DIR}/bams\""

    # Deduplicated alignment folders
    eval "${alignment_type}_DEDUP_DIR=\"\${${alignment_type}_ALIGN_DIR}/deduped_alignments\""

    for dedup_type in "REM" "MARK"; do
        if [[ "$dedup_type" == "REM" ]]; then
            dedup_type_tag="removed"
        else
            dedup_type_tag="marked"
        fi

        # Dedup subfolder base
        eval "${alignment_type}_DEDUP_${dedup_type}_DIR=\"\${${alignment_type}_DEDUP_DIR}/${dedup_type_tag}\""

        # Define subfolders within each dedup mode
        eval "${alignment_type}_DEDUP_${dedup_type}_BAMS=\"\${${alignment_type}_DEDUP_${dedup_type}_DIR}/bams\""
        eval "${alignment_type}_DEDUP_STATS_${dedup_type}=\"\${${alignment_type}_DEDUP_${dedup_type}_DIR}/deduped_stats\""
        eval "${alignment_type}_DEDUP_BIGWIG_${dedup_type}=\"\${${alignment_type}_DEDUP_${dedup_type}_DIR}/bigwig\""
        eval "${alignment_type}_DEDUP_BEDGRAPH_${dedup_type}=\"\${${alignment_type}_DEDUP_${dedup_type}_DIR}/bedgraph\""

        # Create the folders
        eval "mkdir -p \
            \${${alignment_type}_DEDUP_${dedup_type}_BAMS} \
            \${${alignment_type}_DEDUP_STATS_${dedup_type}} \
            \${${alignment_type}_DEDUP_BIGWIG_${dedup_type}} \
            \${${alignment_type}_DEDUP_BEDGRAPH_${dedup_type}}"

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