from mpi4py import MPI
import numpy as np
import pandas as pd


def fixDataframe(array):
    array = comm.gather(array, root=0)

    array2 = np.sum(array, axis=0)
    return array2



comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

temp = np.zeros((30,5))
data = pd.DataFrame(temp)

rows = [rank + size * i for i in range(int(30/size)+1) if comm.rank + comm.size*i < 30]

for i in rows:
    print i
    j = np.ones((1,5))
    temp[i,:]=i*j
    data.iloc[i] = i*j

comm.Barrier()
#print temp
#print data
#data = comm.gather(data.values, root=0)
#data2 = np.sum(data, axis=0)
#print data2

temp2 = fixDataframe(data)
print temp2
#print 'before gather'
#temp = comm.gather(temp,root=0)
#
#print temp
#
#if rank ==0:
#    temp2 = np.sum(temp,axis=0)
#    print 'temp2'
#    print temp2
