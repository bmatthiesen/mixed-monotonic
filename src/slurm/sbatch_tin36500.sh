#!/bin/bash

#SBATCH -J "tin"
#SBATCH --array=601,617,651,704-708,710,712,718,721,723-726,733,734,736,738,744,746-749,752,754,758,760,763,766,768,770,775,777,781,782,784-788,790-792,795,797,799,803,814,820,845
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=36500
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/logPA/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/benchPA

ml load h5py

srun python3 run_tin2.py
