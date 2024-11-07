#!/bin/bash

#SBATCH --mem=20G
#SBATCH --time=10-00:00:00
#SBATCH --job-name=sir
#SBATCH -q long
##SBATCH --cpus-per-task=10
#SBATCH --cpus-per-task=1

#SBATCH --job-name=vorcands
#SBATCH --output=%x-%A-%3a.out
#SBATCH --error=%x-%A-%3a.err

## Array range doesn't have to start at 0.
#SBATCH --array=1-800

echo "Hello from $(hostname)."

echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

## module load Mamba/23.1.0-4
module load conda/latest
module load r/4.4.0

#conda activate nebs05
conda activate vorcands

Rscript R/generic_sim.R $SLURM_ARRAY_TASK_ID
