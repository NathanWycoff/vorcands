#!/bin/bash
##
## Minimal working example of passing the task ID to your program in a job array.
## This is useful when the input to your code can be an integer or the parameter
## you're sweeping over can be mapped to / indexed by an integer.
##

##SBATCH --account=group-rci
#SBATCH --partition=unsafe
##SBATCH --partition=priority
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=4-00:00:00

#SBATCH --job-name=task-id-as-input
## %A is the slurm job number, and %a is the array number within the job.
## %x is the job name.
#SBATCH --output=%x-%A-%3a.out
#SBATCH --error=%x-%A-%3a.err

## Array range doesn't have to start at 0.
#SBATCH --array=0-800

# Print out which node we are running on. Each job in an array might run on
# a different node.
echo "Hello from $(hostname)."

# Print the task ID
# ${SLURM_ARRAY_TASK_ID} is the job array number set by slurm (0, 1, or 2, in this example).
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
module load Mamba/23.1.0-4
. ~/.bashrc
mamba activate nebs05

Rscript R/generic_sim.R $SLURM_ARRAY_TASK_ID
