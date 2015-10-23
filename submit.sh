#!/bin/bash
#SBATCH -J spinice
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 03:59:00

cd /home/ananduri/kmc

prod/a $cellsize $S $step $T $Jdis $Ddis $SLURM_ARRAY_TASK_ID
