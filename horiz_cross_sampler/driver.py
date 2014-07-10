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
    parallel = not mpi_comm==None

    #List of possible observables
    poss_observables = ['PR', 'n_s', 'alpha_s',
            'r', 'n_t',
            'f_NL', 'tau_NL']

    #Which observables do we want?
    obs_to_calc = ['n_s', 'alpha_s', 'f_NL']

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

    #Use to force the histogram to give same number of bins over some pre-defined
    #region in observable space
    fixed_bins=True
    if fixed_bins:
        obs_range = {'n_s': [0.88, 0.965],
                'alpha_s': [-1.0e-2,-1.0e-3],
                'f_NL': [-1.0e-2,-5.0e-3]
                }
        fixed_range = [obs_range[obs] for obs in sorted(obs_to_calc)] #Sort bc in alphab order later
        nbins = 10
        print fixed_range


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


    if parallel and mpi_size>1:
        if mpi_rank==0:
            #Chunk the loop_params into pieces so each process can loop over a subset.
            loop_params = par.chunk(np.array(loop_params),mpi_size,group=False)

        #Halt processors until run parameters are chunked.
        mpi_comm.barrier()

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

            nsamples=5000

            print "nsamples=", nsamples

            #nmoduli = axions + dilaton + heavy moduli
            #nmoduli = nfields + 1 + 1.0*nfields
            nmoduli = nfields/beta

            run = cosmo.SR_universe(sampler="MP_and_horizcross",HC_approx=True,
                    model="Nquad", nfields=nfields)
            radius = 2.0*np.sqrt(run.N_pivot)

            sample = run.sample_Nquad(obs_to_calc, nsamples, nmoduli, radius, m_avg, dimn_weight)

            if fixed_bins:
                hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                        processing.hist_estimate_pdf(sample,normed=False,
                                datarange=fixed_range, nbins=nbins)
            else:
                hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                        processing.hist_estimate_pdf(sample,normed=False,
                                bin_method=processing.scott_rule)

    #for i in hist_total:
    #    print "this is count:", i['counts']
    #    print "this is edges:", i['edges']

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
