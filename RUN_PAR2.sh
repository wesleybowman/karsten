#$ -S /bin/bash
#$ -cwd
#$ -N dngridoctre
#$ -R y
#$ -j y
#$ -l h_rt=47:59:00
#$ -pe ompi 64
gridname=$gridName
time mpirun ./fvcom2d_fast  --casename=$gridname>run_output
echo "Application ends at `date`"
