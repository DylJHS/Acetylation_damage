#!/bin/bash
# This is the workflow configuration script for the Fly Acetylation Damage project.

# Conda environment
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/envs/deeptools_env

# Load shared config
source /hpc/shared/onco_janssen/dhaynessimmons/projects/Dros_H3K9ac_bulkChIC_Analysis/scripts/config/shared_config.sh