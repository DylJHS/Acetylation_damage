#!/bin/bash
#SBATCH --job-name=fly_bowtie2-align
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/logs/fly_bowtie2-align-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/logs/fly_bowtie2-align-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G  # Adjusted for better efficiency
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

FILE_SET=$(ls $TRIMMED_MERG_DIR/*_merged_R1.fq.gz)
FQ_FILE= "${FILE_SET[$SLURM_ARRAY_TASK_ID]}"

# For each sample 
# Check if the sample exists
if [[ ! -e "$FQ_FILE" ]]; then
    echo "No FASTQ file found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi


# if [[ -d $SCC_DIR ]]; then
#     SAMPLE_ID=$(basename $SCC_DIR)

#     # Identify R1 and R2 files
#     R1=$(ls ${SCC_DIR}/*_R1_001.fastq.gz 2>/dev/null | head -n 1)
#     R2=$(ls ${SCC_DIR}/*_R2_001.fastq.gz 2>/dev/null | head -n 1)

#     # Ensure both files exist before proceeding
#     if [[ -z "$R1" || -z "$R2" ]]; then
#         echo "Error: Missing R1 or R2 for ${SAMPLE_ID}" >> "${OUTPUT_DIR}/alignment_errors.log"
#         continue
#     fi

#     # Define output file paths
#     OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.bam"
#     OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_ID}.sam"

#     echo "Processing: ${SAMPLE_ID}"

#     # Run Bowtie2 alignment
#     bowtie2 -x $INDEX_PREFIX -1 $R1 -2 $R2 --threads 8 -S $OUTPUT_SAM

#     # Convert SAM to sorted BAM
#     samtools view -@ 8 -bS $OUTPUT_SAM | samtools sort -@ 8 -o $OUTPUT_BAM
#     samtools index $OUTPUT_BAM


#     echo "Finished aligning ${SAMPLE_ID}"
# fi

# echo "All alignments completed at $(date)"
