import os

gridName = os.environ['gridName']
processors = os.environ['processors']
exe = os.environ['exe']

runpar = '''#$ -S /bin/bash
#$ -cwd
#$ -N {0}
#$ -R y
#$ -j y
#$ -l h_rt=47:59:00
#$ -pe ompi {1}
gridname={0}
time mpirun ./{2}  --casename={0}>run_output
echo "Application ends at `date`"
'''.format(gridName, processors, exe)

outputFile = 'RUN_PAR2.sh'

runpar = runpar.split('\n')

with open(outputFile, 'w') as f:
    for line in runpar:
        print >> f, line
