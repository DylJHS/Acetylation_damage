#!/bin/bash
#SBATCH --job-name=frag_bed
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/frag_bed-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/frag_bed-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

# Get the arguments for species
if [[ "$1" == "-human" ]]; then
    DEDUP_PATH="$HUMAN_DEDUP_DIR"
    my_genome="${HUMAN_REF}/human.genome"
elif [[ "$1" == "-drosophila" ]]; then
    DEDUP_PATH="$DROS_DEDUP_DIR"
    my_genome="${DROS_REF}/drosophila.genome"
elif [[ "$1" == "-tagged" ]]; then
    DEDUP_PATH="$TAGGED_DEDUP_DIR"
    my_genome="${DROS_REF}/drosophila.genome"
else
    echo "Error: Please specify 'human' or 'drosophila' as the first argument."
    exit 1
fi

echo "Processing ${1#-} BAM files..."
echo "-----------------------------------------------------------"

# Get the list of the 
BAM_LIST=("$DEDUP_PATH"/*.bam)
BAM_FILE="${BAM_LIST[$SLURM_ARRAY_TASK_ID]}"
name=$(basename "$BAM_FILE" .bam)

# Create the output folder
OUTPUT_FLD="${DEDUP_PATH}/fragments"
mkdir -p $OUTPUT_FLD

output_file="$OUTPUT_FLD/$name"

echo "Generating the fragment coverage for $name"
echo "-----------------------------------------------------------"

bedtools bamtobed -bedpe -i $BAM_FILE | \
awk '$1==$4 && ((($3 > $6) ? $3 : $6) - (($2 < $5) ? $2 : $5) < 1000)'  | \
cut -f 1,2,6 | \
sort -k1,1 -k2,2n -k3,3n | \
bedtools genomecov -bg -i - -g $my_genome > $output_file.fragments.bedgraph

echo -e "Finished generating the fragment coverage and saved to: \n\t $output_file"
echo "-----------------------------------------------------------"