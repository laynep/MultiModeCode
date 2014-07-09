#!/usr/bin/python

import numpy as np

"""Module that defines methods for handling data use with mpi4py."""

def init_parallel():
    """Parallelize with MPI.  Returns comm (the set of all processes); size (total # of processes); rank (the # of a given process); and name (the processor's name)."""

    try:
        from mpi4py import MPI
    except:
        print "Couldn't find mpi4py.  Continuing without parallelization."
        return None, None, None, None

    #Processor data
    comm=MPI.COMM_WORLD #Set of all processes.
    size=comm.Get_size()
    rank=comm.Get_rank()
    name=MPI.Get_processor_name()
    print "Process %d of %d on %s. \n" %(rank, size-1,name)
    return comm, size, rank, name



def chunk(data,processes,group=True):
    """Breaks a numpy array 'data' into list of smaller arrays with 'processes' elements for scattering with MPI.  Defaults to assigning chunks in blocks over array 'group=True'."""
    #largely from http://stackoverflow.com/questions/12812422/how-can-i-send-part-of-an-array-with-scatter

    #Initialize list with processes numb of elmnts
    chunks = [[] for _ in range(processes)]

    #Convert numpy array to list.
    nlist = data.tolist()

    for i, row in enumerate(nlist):
        if not group:
            #Assign elements of nlist to chunks iteratively
            #i.e. first element to first chunk, second to second, etc...
            chunks[i % processes].append(row)

        else:
            #Build the chunks by assigning part of nlist to first chunk,
            #then next part to second chunk, etc.
            #Last chunk may be smaller.
            chunksize = len(nlist)/processes+1
            chunks[i / chunksize].append(row)

    #Convert each element back to numpy array separately.
    chunks=[np.array(chunks[i]) for i in range(processes)]

    return chunks
