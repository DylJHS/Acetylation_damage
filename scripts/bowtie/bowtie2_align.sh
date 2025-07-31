#!/bin/bash
#SBATCH --job-name=bowtie2-align
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bowtie2-align-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/bowtie2-align-%j.err
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
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/wkflw_config.sh

# Define the files and directories
FILE_SET=("$TRIMMED_MERG_DIR"/*_merged_R1.fq.gz)
echo "-----------------------------------------------------------------------"
echo "Processing file set: ${FILE_SET[*]}"
FQ_FILE="${FILE_SET[$SLURM_ARRAY_TASK_ID]}"
echo "-----------------------------------------------------------------------"
echo "Selected FASTQ file: $FQ_FILE"
SAMPLE_ID=$(basename "$FQ_FILE" _merged_R1.fq.gz)

if [[ "$1" == "-drosophila" ]]; then
    INDEX_PREFIX="${DROS_INDEX}/BDGP6"  # Change for species genome index
    OUTPUT_DIR="${DROS_ALIGN_DIR}"  # Change for species output directory
elif  [[ "$1" == "-human" ]]; then
    INDEX_PREFIX="${HUMAN_INDEX}/GRCh38"  # Change for species genome index
    OUTPUT_DIR="${HUMAN_ALIGN_DIR}"  # Change for species output directory
fi
echo "-----------------------------------------------------------------------"
echo "Running bowtie with the ${1#-} ref genome"

# Define output file paths
OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.bam"
OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_ID}.sam"

# Check if the sample exists
if [[ ! -e "$FQ_FILE" ]]; then
    echo "-----------------------------------------------------------------------"
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
    echo "-----------------------------------------------------------------------"
    echo "Missing R1 or R2 for sample $SAMPLE_ID"
    exit 1
fi

# Run Bowtie2 alignment
echo "Aligning ${SAMPLE_ID} with Bowtie2"
bowtie2 -p 8 -q --no-unal --local --sensitive-local -N 1 -x $INDEX_PREFIX -1 $R1 -2 $R2 -S $OUTPUT_SAM

# Convert SAM to sorted BAM
echo "-----------------------------------------------------------------------"
echo "Converting SAM to sorted BAM for ${SAMPLE_ID}"
samtools view -@ 8 -bS $OUTPUT_SAM | samtools sort -@ 8 -o $OUTPUT_BAM
samtools index $OUTPUT_BAM
rm $OUTPUT_SAM  # Clean up SAM file to save space

echo "-----------------------------------------------------------------------"
echo "Finished aligning ${SAMPLE_ID}"
echo "-----------------------------------------------------------------------"