#!/bin/bash
#SBATCH --job-name=diffbind
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/diffbind-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/diffbind-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/deeptools_env_config.sh

echo "Running DiffBind analysis script ..." 
echo "-------------------------------------------------------"

# Set defaut args
outdir="${PROJ_DIR}/results/tagged_alignments/aligned_bams/diffbind"
sample_sheet="${PROJ_DIR}/data/inputs/diffbind_samplesheet_3.csv"
cofactor=""

# parse the args
while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outdir)
        outdir="${PROJ_DIR}/results/$2"
        shift 2
        ;;
    -s|--sample-sheet)
        sample_sheet="${PROJ_DIR}/data/inputs/$2"
        shift 2
        ;;
    -c|--cofactor)
        cofactor="$2"
        shift 2
        ;;
    -h|--help|--?)
        echo "Valid argumetns are:"
        echo "-o|--outdir <output directory> (default: ${outdir})"
        echo "-s|--sample-sheet <sample sheet> (default: ${sample_sheet})"
        echo "-c|--cofactor <cofactor> (default: ${cofactor})"
        exit 1
        ;;
    --)
        shift
        break
        ;;
    *)
        echo "Unknown argument: $1"
        echo "For help, run: $0 --help or $0 -h or $0 --?"
        exit 1
        ;;
    esac
done

echo "Output directory: $outdir"
echo "Sample sheet: $sample_sheet"
echo "Cofactor: $cofactor"
echo "-------------------------------------------------------"

# Run DiffBind
R script /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/deepTools/run_diffbind.R \
    "$outdir" \
    "$sample_sheet" \
    "$cofactor"

# Check if the DiffBind script ran successfully
if [[ $? -ne 0 ]]; then
    echo "DiffBind analysis failed. Please check the logs for details."
    exit 1
else
    echo "DiffBind analysis completed successfully."
fi

