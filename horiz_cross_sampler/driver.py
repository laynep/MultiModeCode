#!/usr/bin/python

import numpy as np
import cosmology as cosmo
import sys
import processing
import cPickle
import mpi_routines as par

def main():
    """Simple driver function."""

    #MPI parallelized
    mpi_comm, mpi_size, mpi_rank, mpi_name = par.init_parallel()
    parallel = True if mpi_comm!=None else False

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


    #Set up grid of hyperparameters

    #For Marcenko-Pastur distribution
    #beta = naxions/nmoduli
    beta_ratio_max = 0.6
    beta_ratio_min = 0.4
    beta_ratio_numb = 3
    beta_list = np.linspace(beta_ratio_min,beta_ratio_max,beta_ratio_numb)

    #<m_avg^2> = sigma^2 for GRM w/entries of std sigma
    m_avg = 5e-7

    if not parallel or mpi_rank==0:

        nfields_max = 5
        nfields_min = 2
        nfields_unit = 1
        nfields_list = np.arange(nfields_min,nfields_max+1,nfields_unit)

        #For initial conditions:
        #Uniform weighting of dimensions for ICs
        ic_weight = [np.ones(f_numb) for f_numb in nfields_list]

        #Iterate over the rows in this list
        loop_params = zip(ic_weight, nfields_list)

    else:

        loop_params = None


    if parallel and mpi_rank==0:
        #Chunk the loop_params into pieces so each process can loop over a subset.
        loop_params = par.chunk(np.array(loop_params),mpi_size,group=False)

    #Halt processors until run parameters are chunked.
    if parallel: mpi_comm.barrier()


    if parallel and mpi_size>1:
        #Scatter loop_params to all processes
        loop_params = mpi_comm.scatter(loop_params,root=0)


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
    for dimn_weight, nfields in loop_params:
        for beta in beta_list:
            print "\nNew run:"
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

            run = cosmo.SR_universe(sampler="MP_and_horizcross",HC_approx=True,
                    model="Nquad", nfields=nfields)
            radius = 2.0*np.sqrt(run.N_pivot)

            sample = run.sample_Nquad(obs_to_calc, nsamples, nmoduli, radius, m_avg, dimn_weight)


            hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                    processing.hist_estimate_pdf(sample,normed=False, bin_method=processing.scott_rule)


    #for item in hist_total:
    #    print "nfields:"
    #    print item['nfields']
    #    print "these are the counts:"
    #    print item['counts']
    #    print "these are the edges:"
    #    print item['edges']

    if parallel:
        myfile = open("outdata"+str(mpi_rank)+".dat","w")
    else:
        myfile = open("outdata.dat","w")
    cPickle.dump(hist_total, myfile)
    myfile.close()





if __name__=="__main__":

    profile=False
    if profile:
        import cProfile
        cProfile.run('main()')
    else:
        main()
