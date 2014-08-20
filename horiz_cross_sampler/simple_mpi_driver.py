#!/usr/bin/python

"""Program to sample the horizon crossing surface."""

import numpy as np
import cosmology as cosmo
import pandas as pd
import mpi_routines as par
import sys


#nfields= np.arange(2.0, 502.0, 100.0)
nfields= np.arange(5.0, 500.0, 5.0)
#nfields= np.arange(2.0, 3.0, 1.0)


nfields = list(nfields)
nfields.append(2.0)
nfields = np.sort(np.array(nfields))


p=2.0
#p=2.0/3.0
#p=1.5
#p=1.0
#p=4.0
#p=8.0

#sampler='log'
#sampler='uniform'
sampler='MP_and_horizcross'

#low=-14
#high=-12

low=1e-14
high=1e-13

beta = 0.5
m_avg=5e-7

nsamples=1000

#obs_to_calc = ['n_t','r']
obs_to_calc = ['PR','alpha_s','n_s','n_t','r']


#other_params = ['nfields', 'p', 'low', 'high']
other_params = ['nfields', 'beta', 'm_avg']



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

    run = cosmo.Nmono_universe(sampler=sampler,HC_approx=True,
            model="Nmono", nfields=f_numb)
    radius = np.sqrt(2.0*p*run.N_pivot)

    if sampler=='MP_and_horizcross':
        nmoduli = f_numb/beta
        sample_total = run.sample_Nmono_MP(obs_to_calc, nsamples,
            nmoduli, radius, m_avg, weight, p)
    else:
    	sample_total = run.sample_Nmono_uniform(obs_to_calc, nsamples,
    	        low, high, radius, weight, p)



    for row in sample_total:
        for val in obs_to_calc:
            output[val].append(row[val])
        output['nfields'].append(int(f_numb))
        #output['p'].append(p)
        #output['low'].append(low)
        #output['high'].append(high)
        output['beta'].append(beta)
        output['m_avg'].append(m_avg)


output = pd.DataFrame(output, columns=obs_to_calc+other_params)

output.save('output'+str(mpi_rank)+'.pkl')

print output
