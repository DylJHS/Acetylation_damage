#!/bin/bash
# Directory containing the trimmed FASTQ files
FASTQ_DIR="/home/djhs/fly_acetylation_damage/data/fastq/trimmed"
FASTA_DIR="/home/djhs/fly_acetylation_damage/data/trimmed_fasta"

# Number of reads to sample from each file
SAMPLE_SIZE=100

# Loop over all gzipped FASTQ files that include "_val_" in the filename (to avoid reports)
for READS in ${FASTQ_DIR}/*_val_*.f*q.gz; do
    # Get the base name (strip the path)
    base=$(basename "$READS")
    # Remove the .fq.gz or .fastq.gz extension
    base_name=${base%%.fq.gz}
    if [ "$base_name" = "$base" ]; then
        base_name=${base%%.fastq.gz}
    fi

    # Define output file names for the subset FASTQ and FASTA files
    SUBSET_FA="${FASTA_DIR}/${base_name}_subset.fasta"

    echo "Processing file: $READS"
    echo "  -> Subset FASTA: $SUBSET_FA"

    # Sample a subset of reads (decompress on the fly with zcat) and write to subset FASTQ
    zcat "$READS" | seqtk sample -s100 - "$SAMPLE_SIZE" | seqtk seq -a > "$SUBSET_FA"
done
