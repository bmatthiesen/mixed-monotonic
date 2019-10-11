#!/bin/bash

#SBATCH -J "gee"
#SBATCH --array=0-99
#SBATCH --account=p_mwrc
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2583
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/log/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/bench

ml load h5py

srun python3 run_gee.py 5
