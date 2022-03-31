#!/bin/bash

#SBATCH --job-name=AlphaTesting
#SBATCH --mincpus=1
#SBATCH --ntasks=1
#SBATCH --mem=5GB
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=guangyao
#SBATCH --array 1-100

module load R

Rscript ArrayAlphaTesting.R ${SLURM_ARRAY_TASK_ID}