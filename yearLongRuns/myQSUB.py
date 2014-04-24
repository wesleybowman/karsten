import subprocess
import re
import time

commandQsub = subprocess.check_output(['qsub', 'RUN_PAR2.sh'])
#print commandQsub

pattern = r'\b\d{1,9}\b'
match = re.search(pattern,commandQsub)

jobID = match.group(0)
#print jobID

match=True
pattern = r'\b{}\b'.format(jobID)

while match:
    commandQstat = subprocess.check_output(["qstat"])
    match = re.search(pattern,commandQstat)
    #print 'Still running'
    time.sleep(120)
