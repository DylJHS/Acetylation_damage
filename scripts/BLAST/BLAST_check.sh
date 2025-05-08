#!/bin/bash
# Set file names and parameters
READS="/home/djhs/fly_acetylation_damage/data/SCC-bulkChIC-UMC-JAN-004/SCC-bulkChIC-UMC-JAN-004_AAG7VNLM5_S8_L001_R2_001.fastq.gz"  # Your original gzipped FASTQ file
SUBSET_FA="/home/djhs/fly_acetylation_damage/data/FASTA/JAN004R2_subset.fasta"      # FASTA version of the subset
SAMPLE_SIZE=500               # Number of reads to sample

# Step 1: Sample a subset of reads using seqtk (for gzipped file, use zcat)
zcat "$READS" | seqtk sample -s100 - "$SAMPLE_SIZE" | seqtk seq -a > "$SUBSET_FA"

