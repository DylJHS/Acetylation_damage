#!/bin/bash
#SBATCH --job-name=run_blast_wkflw
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/run_blast_wkflw-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/logs/run_blast_wkflw-%j.err
#SBATCH --time=18:00:00
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/BLAST/BLAST_wkflw_config.sh

# Submit FASTA subset generation (array: 1 job per JAN folder)
jid_fasta=$(sbatch --parsable --array=0-5 "$PROJ_DIR/scripts/BLAST/gen_fasta_subsets_arr.sh")
echo "Submitted FASTA subset generation as job $jid_fasta"

# Wait until FASTA files are created (optional: add a wait loop here)

# Count how many *_subset.fasta files were generated
sleep 5
count=$(ls "$FASTA_DIR"/*_subset.fasta 2>/dev/null | wc -l)
echo "Detected $count FASTA subset files"

# Submit BLASTn query job as array
jid_blast=$(sbatch --parsable --dependency=afterok:$jid_fasta --array=0-$((count - 1)) "$PROJ_DIR/scripts/BLAST/BLASTn_arr.sh")
echo "Submitted BLASTn job as array job $jid_blast"

# Submit summary generation
jid_summary=$(sbatch --parsable --dependency=afterok:$jid_blast "$PROJ_DIR/scripts/BLAST/gen_blast_summary.sh")
echo "Submitted summary generation as job $jid_summary"
