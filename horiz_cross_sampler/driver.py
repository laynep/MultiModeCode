#!/usr/bin/python

"""
Simple driver function.
"""

import numpy as np
from classes import *
import sys

def main():

    nfields=1
    y = universe(sampler="MP_and_uniformsphere",HC_approx=True,
            model="Nquad", nfields=nfields )
    radius = 2.0*np.sqrt(y.N_pivot)

    obs_to_calc = ('n_t', 'n_s', 'alpha_s')

    # beta = naxions/nmoduli
    # nmoduli = axions + dilaton + heavy moduli
    nmoduli = nfields + 1 + 1.0*nfields

    # <m_avg^2> = sigma^2 for GRM w/entries of std sigma
    m_avg = 5e-7

    # How many samples to build PDF from
    nsamples = 100

    for i in xrange(nsamples):

        y.get_new_params( nmoduli=nmoduli, radius=radius, m_avg=m_avg)

        y.xi_i(y.phi_hc)

        y.calc_observs(y.phi_hc,obs_to_calc=obs_to_calc)

        print y.observ
        #print y.params["Nquad"]["m2"]
        #print y.params/np.min(y.params["Nquad"]["m2"])
        #print y.params["Nquad"]["m2"]/np.min(y.params["Nquad"]["m2"])
        #print y.observ['n_t']/( -y.observ['r']/8.0)


if __name__=="__main__":

    profile=False
    if profile:
        import cProfile
        cProfile.run('main()')
    else:
        main()
