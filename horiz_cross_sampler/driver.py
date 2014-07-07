#!/usr/bin/python

import numpy as np
from classes import *
import sys
import processing as proc


def main():
    """Simple driver function."""


    nfields=100
    run = SR_universe(sampler="MP_and_uniformsphere",HC_approx=True,
            model="Nquad", nfields=nfields )
    radius = 2.0*np.sqrt(run.N_pivot)

    obs_to_calc = ('n_t', 'n_s', 'alpha_s')

    # beta = naxions/nmoduli
    # nmoduli = axions + dilaton + heavy moduli
    nmoduli = nfields + 1 + 1.0*nfields

    # <m_avg^2> = sigma^2 for GRM w/entries of std sigma
    m_avg = 5e-7

    # How many samples to build PDF from
    nsamples = 10

    sample = run.sample_Nquad(obs_to_calc, nsamples, nmoduli, radius, m_avg)


    print sample

    proc.build_pdf_hist(sample)




if __name__=="__main__":

    profile=False
    if profile:
        import cProfile
        cProfile.run('main()')
    else:
        main()
