#!/bin/bash
#SBATCH --job-name=tagged_transformation
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/tagged_transformation-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/tagged_transformation-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Transform the names of the taggged bam files and save them to the same directory

# Load Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/wkflw_config.sh

FODLERS=("$DATA_DIR"/bulkChIC*)
TAGGED_FLDR="$DATA_DIR/tagged_bam_files"

mkdir -p "$TAGGED_FLDR"

for FOLDER in "${FODLERS[@]}"; do
    FOLD_NAME=$(basename "$FOLDER")
    echo "Processing folder: $FOLD_NAME"
    for BAM_FILE in "$FOLDER"/tagged*; do
        if [[ -e "$BAM_FILE" ]]; then
            BASENAME=$(basename "$BAM_FILE")
            NEW_NAME="${FOLD_NAME}_${BASENAME}"
            NEW_PATH="$TAGGED_FLDR/$NEW_NAME"
            echo "Renaming $BAM_FILE to $NEW_PATH"
            mv "$BAM_FILE" "$NEW_PATH"
        else
            echo "No tagged files found in $FOLDER"
        fi
    done
done