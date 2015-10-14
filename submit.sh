#!/bin/bash
#SBATCH -J spinice
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 23:59:00

cd /home/ananduri/kmc

anneal/a $cellsize $S $step $T $SLURM_ARRAY_TASK_ID
