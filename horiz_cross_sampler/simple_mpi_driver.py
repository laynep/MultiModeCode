#!/usr/bin/python

"""Program to sample the horizon crossing surface."""

import numpy as np
import cosmology as cosmo
import pandas as pd
import mpi_routines as par
import sys


nfields= np.arange(2.0, 60.0, 1.0)
p=1.5

low=-14
high=-12

nsamples=2

obs_to_calc = ['n_t','r']


other_params = ['nfields', 'p', 'low', 'high']



#MPI parallelized
mpi_comm, mpi_size, mpi_rank, mpi_name = par.init_parallel()


if mpi_rank==0:

    #Uniform weighting of dimensions for ICs
    ic_weight = [np.ones(f_numb) for f_numb in nfields]

    #Iterate over the rows in this list
    loop_params = zip(ic_weight, nfields)

else:

    loop_params = None



if mpi_size>1:
    if mpi_rank==0:
        #Chunk the loop_params into pieces so each process can loop over a subset.
        loop_params = par.chunk(np.array(loop_params),mpi_size,group=False)

    #Halt processors until run parameters are chunked.
    mpi_comm.barrier()

    #Scatter loop_params to all processes
    loop_params = mpi_comm.scatter(loop_params,root=0)


output = {val:[] for val in obs_to_calc+other_params}
for weight, f_numb in loop_params:

    run = cosmo.Nmono_universe(sampler="log",HC_approx=True,
            model="Nmono", nfields=f_numb)
    radius = np.sqrt(2.0*p*run.N_pivot)

    sample_total = run.sample_Nmono_uniform(obs_to_calc, nsamples,
            low, high, radius, weight, p)



    for row in sample_total:
        for val in obs_to_calc:
            output[val].append(row[val])
        output['nfields'].append(int(f_numb))
        output['p'].append(p)
        output['low'].append(low)
        output['high'].append(high)


output = pd.DataFrame(output, columns=obs_to_calc+other_params)

output.save('output'+str(mpi_rank)+'.pkl')

print output
