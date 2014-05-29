from mpi4py import MPI
import numpy as np
import pandas as pd


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
    data = np.ones((8,1))
    data = np.arange(20)
    #data = pd.DataFrame(data)
    print data
else:
   data = None
data = comm.bcast(data, root=0)
#data = comm.scatter(data, root=0)
#data[0] = (rank+1)
data = (rank+1)
comm.Barrier()
data = comm.gather(data, root=0)
if rank==0:
    print data
