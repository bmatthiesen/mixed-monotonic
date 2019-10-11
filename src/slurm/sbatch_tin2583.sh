#!/bin/bash

#SBATCH -J "tin"
#SBATCH --array=0-516,518-599,602-606,608,610-616,618-626,628,630-638,640-650,652-657,659-661,663-674,676-680,682-684,686,687,689-692,694-699
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=2583
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/logPA/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/benchPA

ml load h5py

srun python3 run_tin2.py
