#!/bin/bash
#SBATCH --job-name=faidx_gen
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/faidx_gen-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/faidx_gen-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/genomics_env_config.sh

# get the arguments for species


for SPECIES in "human" "drosophila"; do
    if [[ $SPECIES == "drosophila" ]]; then
        REFERENCE="$DROS_GEN_FA"
        OUTPUT_DIR="${DROS_REF}"
    elif [[ "$SPECIES" == "human" ]]; then
        REFERENCE="$HUMAN_GEN_FA"
        OUTPUT_DIR="${HUMAN_REF}"
    fi
    
    OUTPUT="${OUTPUT_DIR}/${SPECIES}.fai"
    echo "-----------------------------------------------------------"
    echo "Processing $REFERENCE reference genome"
    echo "-----------------------------------------------------------"

    samtools faidx "${REFERENCE}" -o "${OUTPUT}"
    cut -f1,2 "${OUTPUT}" > "${OUTPUT_DIR}/${SPECIES}.genome"

    echo "Done and saved to: ${OUTPUT_DIR}/${SPECIES}.genome " 
    echo "-----------------------------------------------------------"
done