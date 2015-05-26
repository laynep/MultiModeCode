#!/usr/bin/python

"""Program to sample the horizon crossing surface."""

import numpy as np
import cosmology as cosmo
import pandas as pd

nfields=100
p=2.0

nsamples=100
ic_weight = np.ones(nfields)

low=1e-14
high=1e-13

obs_to_calc = ['n_s','r']

run = cosmo.Nmono_universe(sampler="uniform",HC_approx=True,
        model="Nmono", nfields=nfields)
radius = np.sqrt(2.0*p*run.N_pivot)

sample_total = run.sample_Nmono_uniform(obs_to_calc, nsamples,
        low, high, radius, ic_weight, p)

output = {val:[] for val in obs_to_calc}
for row in sample_total:
    for val in output:
        output[val].append(row[val])

output = pd.DataFrame(output, columns=obs_to_calc)

print output

output.save('simple_output'+str(nfields)+'.pkl')
