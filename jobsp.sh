#for x in 0 # {0..40}
#do
#sbatch --export=cellsize=6,S=40,step=0,T=1.5,label=$x submit.sh
#done

sbatch --export=cellsize=6,S=40,step=0,T=1.5 --array=3-40 submit.sh
