#!/bin/bash
#SBATCH -J spinice8
#SBATCH -p New
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 48:59:00

cd /home/ananduri/kmc/cellsize8

prod/a $S $step $T $label
