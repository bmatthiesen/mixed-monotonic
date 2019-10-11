#!/bin/bash

#SBATCH -J "tin"
#SBATCH --array=607,627,658,662,713,715,719,722,730-732,737,740-743,750,751,756,759,765,769,771,772,774,776,778,780,783,794,796,864,898
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=10583
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/logPA/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/benchPA

ml load h5py

srun python3 run_tin2.py
