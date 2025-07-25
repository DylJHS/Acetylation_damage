#!/bin/bash
#SBATCH --job-name=fly_bowtie2-align
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/fly_bowtie2-align-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/fly_bowtie2-align-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G  # Adjusted for better efficiency
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl
echo "-----------------------------------------------------------------------"
echo "Starting Bowtie2 alignment at $(date)"

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

# Define the files and directories
FILE_SET=("$TRIMMED_MERG_DIR"/*_merged_R1.fq.gz)
FQ_FILE= "${FILE_SET[$SLURM_ARRAY_TASK_ID]}"
RES_PATH=${DROS_ALIGN_DIR}     # Change for species alignment results path
INDEX_PREFIX="${DROS_INDEX}"   # Change for species genome index build
OUTPUT_DIR="${RES_PATH}/bowtie2_alignments"
SAMPLE_ID=$(basename "$FQ_FILE" _merged_R1.fq.gz)

# Define output file paths
OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.bam"
OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_ID}.sam"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Check if the sample exists
if [[ ! -e "$FQ_FILE" ]]; then
    echo "No FASTQ file found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "-----------------------------------------------------------------------"
echo "Processing: ${FQ_FILE}"

# Define the corresopnding R1 and R2 files
R1="${FQ_FILE}"
R2="${FQ_FILE/_R1/_R2}"

# Check that both R1 and R2 files exist
if [[ ! -e "$R1" || ! -e "$R2" ]]; then
    echo "Missing R1 or R2 for sample $SAMPLE_ID"
    exit 1
fi

# Run Bowtie2 alignment
echo "Aligning ${SAMPLE_ID} with Bowtie2"
bowtie2 -x $INDEX_PREFIX -1 $R1 -2 $R2 --threads 8 -S $OUTPUT_SAM

# Convert SAM to sorted BAM
echo "-----------------------------------------------------------------------"
echo "Converting SAM to sorted BAM for ${SAMPLE_ID}"
samtools view -@ 8 -bS $OUTPUT_SAM | samtools sort -@ 8 -o $OUTPUT_BAM
samtools index $OUTPUT_BAM
rm $OUTPUT_SAM  # Clean up SAM file to save space

echo "-----------------------------------------------------------------------"
echo "Finished aligning ${SAMPLE_ID}"
echo "-----------------------------------------------------------------------"