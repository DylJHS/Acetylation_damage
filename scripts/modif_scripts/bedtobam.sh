#!/bin/bash
#SBATCH --job-name=bamtobed
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bamtobed-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bamtobed-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda or module environment with samtools, bedtools, UCSC

source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/deeptools_env_config.sh

# Parse species flag  (-human / -drosophila / -tagged)
if   [[ "$1" == "-human" ]]; then
    BAM_DIR="$HUMAN_DEDUP_DIR"
    GENOME_FILE="$HUMAN_GENOME_SIZES"
elif [[ "$1" == "-drosophila" ]]; then
    BAM_DIR="$DROS_DEDUP_DIR"
    GENOME_FILE="$DROS_GENOME_SIZES" 
elif [[ "$1" == "-tagged" ]]; then
    BAM_DIR="$TAGGED_ALIGNMENT_DIR"
    GENOME_FILE="$DROS_GENOME_SIZES"
else
    echo "Error: first argument must be -human, -drosophila, or -tagged"
    exit 1
fi

#  Determine the BAM file for this array index
BAM_FILES=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
FILE_BASENAME=$(basename "${BAM_FILE%%.bam}")

# Set output paths

OUTPUT_PATH="${BAM_DIR}/coverage"
mkdir -p "${OUTPUT_PATH}"
BEDGRAPH_FILE="${OUTPUT_PATH}/${FILE_BASENAME}.bedgraph"
BIGWIG_FILE="${OUTPUT_PATH}/${FILE_BASENAME}.bw"

echo "[$(date)] Processing ${FILE_BASENAME}.bam → fragment coverage"

# ------------------------------------------------------------------
# 5. OPTIONAL: light BAM filtering to drop low-MAPQ & secondary reads
# ------------------------------------------------------------------
TMP_BAM=$(mktemp --suffix=.bam)
samtools view -@ "$SLURM_CPUS_PER_TASK" -b \
    -F 1804   \     # remove unmapped/secondary/duplicates
    -q 30     \
    "$BAM_FILE" -o "$TMP_BAM"
samtools index "$TMP_BAM"

# ------------------------------------------------------------------
# 6. BEDPE  ➜  fragment BED3  ➜  bedGraph
# ------------------------------------------------------------------
TMP_BEDPE=$(mktemp)
TMP_FRAG=$(mktemp)

# 6a. paired-end BAM → BEDPE (one line/read-pair)
bedtools bamtobed -bedpe -i "$TMP_BAM" > "$TMP_BEDPE"

# 6b. keep pairs on same chrom and ≤1 kb apart, collapse to fragment
awk '$1==$4 && $6-$2<1000 {print $1"\t"$2"\t"$6}' "$TMP_BEDPE" \
  | sort -k1,1 -k2,2n -k3,3n > "$TMP_FRAG"

# 6c. fragment BED3 → coverage bedGraph (raw counts, no zeros)
bedtools genomecov -bg -i "$TMP_FRAG" -g "$GENOME_FILE" > "$BEDGRAPH_FILE"

# ------------------------------------------------------------------
# 7. If requested, convert to BigWig
# ------------------------------------------------------------------
if [[ "$OUTPUT_EXT" == "bw" ]]; then
    sort -k1,1 -k2,2n "$BEDGRAPH_FILE" -o "$BEDGRAPH_FILE"
    bedGraphToBigWig "$BEDGRAPH_FILE" "$GENOME_FILE" "$BIGWIG_FILE"
    echo "[$(date)] Created BigWig: $BIGWIG_FILE"
    # Keep or delete the intermediate bedGraph as you prefer
    # rm "$BEDGRAPH_FILE"
fi

# ------------------------------------------------------------------
# 8. Clean-up
# ------------------------------------------------------------------
rm -f "$TMP_BAM" "$TMP_BAM.bai" "$TMP_BEDPE" "$TMP_FRAG"
echo "[$(date)] Coverage job done for ${FILE_BASENAME}"
