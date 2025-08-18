#!/bin/bash
#SBATCH --job-name=header
#SBATCH --output=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/header-%j.out
#SBATCH --error=/hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/logs/header-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --mail-type=all
#SBATCH --mail-user=d.j.haynes-simmons@umcutrecht.nl

# Load the workflow configuration
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/genomics_env_config.sh

FOLDER="$TAGGED_RAW_ALIGNMENT_DIR_BAMS"

for bam_file in "$FOLDER"/*.bam; do
    echo -e "\nProcessing BAM file: $(basename "$bam_file")"
    samtools view "$bam_file" | awk '{print $0}' | head -n 10
    echo -e "\n"
    echo "----------------------------------------"
    samtools view "$bam_file" \
    | awk '{
        mi="";
        for(i=12;i<=NF;i++){
            if($i ~ /^MI:Z:/){ split($i,a,":"); mi=a[3]; break }
        }
        if(mi!="") c[mi]++
    }
    END{
        dup=0;
        for(mi in c) if(c[mi]>1) dup += (c[mi]-1);
        print dup
    }'
    echo -e "\nDone processing BAM file: $(basename "$bam_file")\n"
done
