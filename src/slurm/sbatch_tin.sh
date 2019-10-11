#!/bin/bash

#SBATCH -J "tin"
#SBATCH --array=0-899
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=2583
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/logPA/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/benchPA

ml load h5py

srun python3 run_tin.py
