#!/bin/bash
#SBATCH -J spinice
#SBATCH -p New
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 23:59:00

cd /home/ananduri/kmc

anneal/a $cellsize $S $step $T $label
