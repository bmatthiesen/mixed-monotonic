#!/bin/bash

#SBATCH -J "tin"
#SBATCH --array=700,701,709,716,717,727-729,739,745,762,764,767,779,793,798,800-802,804-813,815,817-819,821-844,846-863,865-897,899
#SBATCH --account=p_mwrc
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=109500
#SBATCH --partition=haswell
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/p_mwrc/diss/logPA/%x-wp_%a-job_%A.out

export JOB_HPC_SAVEDIR=/scratch/p_mwrc/diss/benchPA

ml load h5py

srun python3 run_tin2.py
