#for x in 0 # {0..40}
#do
#sbatch --export=cellsize=6,S=40,step=0,T=1.5,label=$x submit.sh
#done

sbatch --export=cellsize=4,S=1,step=0,T=3 --array=0 submit.sh
