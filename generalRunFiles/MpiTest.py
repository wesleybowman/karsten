from mpi4py import MPI
import numpy as np


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

temp = np.zeros((30,1))

rows = [rank + size * i for i in range(int(30/size)+1) if comm.rank + comm.size*i < 30]

for i in rows:
    print i
    temp[i]=i

comm.Barrier()
print temp

print 'before gather'
temp = comm.gather(temp,root=0)
print temp
