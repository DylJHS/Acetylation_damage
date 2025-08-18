#!/bin/bash
#SBATCH --job-name=frag_bed
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/frag_bed-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/frag_bed-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

# Set directories based on the alignment type
BAM_DIR="/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/results/${1}"
BAM_FILES=("$BAM_DIR"/*.bam)
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"
name=$(basename "$BAM_FILE" .bam)

# Create the output folder
OUTPUT_FLD="${BAM_DIR}/../fragments"
mkdir -p $OUTPUT_FLD
output_file="$OUTPUT_FLD/$name"

# Determine the genome based on the alignment type
if [[ "$BAM_DIR" == *drosophila* || "$BAM_DIR" == *tagged* ]]; then
    my_genome=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/data/ref_genomes/drosophila/drosophila.genome
elif [[ "$BAM_DIR" == *human* ]]; then 
    my_genome=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/data/ref_genomes/human/human.genome
fi

echo "Processing ${1#-} BAM files..."
echo "-----------------------------------------------------------"

echo "Generating the fragment coverage for $name"
echo "-----------------------------------------------------------"

bedtools bamtobed -bedpe -i $BAM_FILE | \
awk '$1==$4 && ((($3 > $6) ? $3 : $6) - (($2 < $5) ? $2 : $5) < 1000)'  | \
cut -f 1,2,6 | \
sort -k1,1 -k2,2n -k3,3n | \
bedtools genomecov -bg -i - -g $my_genome > $output_file.fragments.bedgraph

echo -e "Finished generating the fragment coverage and saved to: \n\t $output_file"
echo "-----------------------------------------------------------"