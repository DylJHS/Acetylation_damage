#!/bin/bash
# This is the workflow configuration script for the Fly Acetylation Damage project.

# Conda environment
export PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/bin:$PATH"
source /hpc/shared/onco_janssen/dhaynessimmons/envs/miniconda3/etc/profile.d/conda.sh
conda activate /hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env
export LD_LIBRARY_PATH="/hpc/shared/onco_janssen/dhaynessimmons/envs/genomics_env/lib:$LD_LIBRARY_PATH"

# Load shared config
source /hpc/shared/onco_janssen/dhaynessimmons/projects/fly_acetylation_damage/scripts/config/shared_config.sh