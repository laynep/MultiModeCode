#!/usr/bin/python

import numpy as np
from classes import *
import sys
import processing as proc


def main():
    """Simple driver function."""


    nfields=100
    run = SR_universe(sampler="MP_and_uniformsphere",HC_approx=True,
            model="Nquad", nfields=nfields)
    radius = 2.0*np.sqrt(run.N_pivot)

    obs_to_calc = ['n_s', 'alpha_s']

    # beta = naxions/nmoduli
    # nmoduli = axions + dilaton + heavy moduli
    nmoduli = nfields + 1 + 1.0*nfields

    # <m_avg^2> = sigma^2 for GRM w/entries of std sigma
    m_avg = 5e-7

    # How many samples to build PDF from
    nsamples = 100

    sample = run.sample_Nquad(obs_to_calc, nsamples, nmoduli, radius, m_avg)


    counts, edges = proc.build_pdf_hist(sample,['n_s','alpha_s'],normed=False)
    print "these are the counts"
    print counts
    print "these are the edges"
    print edges



if __name__=="__main__":

    profile=False
    if profile:
        import cProfile
        cProfile.run('main()')
    else:
        main()
