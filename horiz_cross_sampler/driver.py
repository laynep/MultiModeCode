#!/usr/bin/python

import numpy as np
from classes import *
import sys
import processing

def main():
    """Simple driver function."""

    #List of possible observables
    poss_observables = ['PR', 'n_s', 'alpha_s',
            'r', 'n_t',
            'f_NL', 'tau_NL']

    #Which observables do we want?
    obs_to_calc = ['n_s']

    #List of possible hyperparameters
    #[ LP: ] Should probably also include N_piv...
    hyperparams = ['nfields', 'beta',
            'm_avg', 'dimn_weight']

    ######
    #Cycle over grid of hyperparameters
    ######


    nfields_max = 10
    nfields_min = 2
    nfields_unit = 1
    nfields_list = list(np.arange(nfields_min,nfields_max+1,nfields_unit))

    #For Marcenko-Pastur distribution
    #beta = naxions/nmoduli
    beta_ratio_max = 0.6
    beta_ratio_min = 0.4
    beta_ratio_numb = 1
    beta_list = np.linspace(beta_ratio_min,beta_ratio_max,beta_ratio_numb)

    #<m_avg^2> = sigma^2 for GRM w/entries of std sigma
    m_avg = 5e-7

    #For initial conditions:
    #Uniform weighting of dimensions for ICs
    ic_weight = [np.ones(f_numb) for f_numb in nfields_list]

    def load_hist_dictionary():
        """Loads the hyper parameters into a dictionary for each run."""
        hist_total.append({})
        hist_total[-1]['beta'] = beta
        hist_total[-1]['nfields'] = nfields
        hist_total[-1]['m_avg'] = m_avg
        hist_total[-1]['dimn_weight'] = dimn_weight
        hist_total[-1]['observs'] = obs_to_calc


    #Iterate over all desired combinations of hyperparameters
    hist_total = []
    for dimn_weight, nfields in zip(ic_weight, nfields_list):
        for beta in beta_list:
            print "New run:"
            print "beta=", beta
            print "nfields=", nfields

            load_hist_dictionary()

            #How many samples to build PDF from
            samp_ref=2e7
            if nfields>=10:
                nsamples = int(np.ceil(samp_ref/nfields**2))
            else:
                nsamples = int(np.ceil(samp_ref/10**2))

            nsamples=1000

            print "nsamples=", nsamples

            #nmoduli = axions + dilaton + heavy moduli
            #nmoduli = nfields + 1 + 1.0*nfields
            nmoduli = nfields/beta

            run = SR_universe(sampler="MP_and_horizcross",HC_approx=True,
                    model="Nquad", nfields=nfields)
            radius = 2.0*np.sqrt(run.N_pivot)

            sample = run.sample_Nquad(obs_to_calc, nsamples, nmoduli, radius, m_avg, dimn_weight)


            hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                    processing.hist_estimate_pdf(sample,normed=False, bin_method=processing.scott_rule)


    for item in hist_total:
        print "nfields:"
        print item['nfields']
        print "these are the counts:"
        print item['counts']
        print "these are the edges:"
        print item['edges']




if __name__=="__main__":

    profile=False
    if profile:
        import cProfile
        cProfile.run('main()')
    else:
        main()
