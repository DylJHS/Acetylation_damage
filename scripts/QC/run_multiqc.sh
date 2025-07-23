#!/bin/bash
#SBATCH --job-name=run_multiqc
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/run_multiqc-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/run_multiqc-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration (defines $DATA_DIR , etc.)
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

cd "$TRIMMED_DIR" 

echo "Running fastqc in directory: $PWD"
# run the fastqc command
fastqc -t 5 -o "$TRIMMED_DIR" *.fq.gz
echo "FastQC completed in directory: $PWD"

echo "-----------------------------------------------------"

echo "Running MultiQC in directory: $PWD"

# run the multiqc command
multiqc .
echo "MultiQC completed in directory: $PWD"