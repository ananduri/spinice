for x in {1..40}
do
sbatch --export=S=10,step=0,T=10,label=$x submit.sh
done
