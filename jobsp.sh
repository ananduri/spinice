for x in 102 # {1..40}
do
sbatch --export=cellsize=6,S=20,step=0,T=2,label=$x submit.sh
done
