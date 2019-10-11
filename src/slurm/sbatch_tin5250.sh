#!/bin/bash

#SBATCH -J "tin"
#SBATCH --array=517,600,609,629,639,675,681,685,688,693,702,703,711,714,720,735,753,755,757,761,773,789,816
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=5250
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/logPA/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/benchPA

ml load h5py

srun python3 run_tin2.py
