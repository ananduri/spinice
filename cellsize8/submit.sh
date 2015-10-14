#!/bin/bash
#SBATCH -J spinice8
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 23:59:00

cd /home/ananduri/kmc/cellsize8

prod/a $S $step $T $SLURM_ARRAY_TASK_ID
