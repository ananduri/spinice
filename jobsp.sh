#for x in 0 # {0..40}
#do
#sbatch --export=cellsize=6,S=40,step=0,T=1.5,label=$x submit.sh
#done

sbatch --export=cellsize=4,S=40,step=0,T=2.0 --array=0-9 submit.sh
