#!/usr/bin/python

"""
Simple driver function.
"""

import numpy as np
from classes import *
import sys
import cProfile

nfields=100
y = universe(sampler="MP_and_uniformsphere",HC_approx=True, model="Nquad",
        nfields=nfields )
radius = 2.0*np.sqrt(y.N_pivot)

obs_to_calc = ('n_t','r', 'n_s', 'alpha_s')

profile=False
if profile:
    cProfile.run(
        'y.get_new_params( nmoduli=1*nfields, radius=radius, m_avg=5.0e-7)'
    )
    cProfile.run(
        "y.calc_observs(y.phi_hc,obs_to_calc=obs_to_calc)"
    )
else:
    for i in xrange(10):

        # beta = naxions/nmoduli
        # nmoduli = axions + dilaton + heavy moduli
        nmoduli = nfields + 1 + 1.0*nfields

        y.get_new_params( nmoduli=nmoduli, radius=radius, m_avg=5.0e-7)

        y.calc_observs(y.phi_hc,obs_to_calc=obs_to_calc)

        #print y.observ
        print y.observ['n_t']/( -y.observ['r']/8.0)
