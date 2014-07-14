#!/usr/bin/python

"""Parameters that would change between runs of the horizon crossing sampler.  Variables that should be available for all the threads.  For serial runs or only one MPI process, all parameters are read by the same thread."""

#List of possible observables
poss_observables = ['PR', 'n_s', 'alpha_s',
        'r', 'n_t',
        'f_NL', 'tau_NL']

#Which observables do we want?
obs_to_calc = ['n_s', 'alpha_s']


#List of possible hyperparameters
#[ LP: ] Should probably also include N_piv...
hyperparams = ['nfields', 'beta',
        'm_avg', 'dimn_weight']


#Range of field number to iterate over
#Will create an array between min-max with stepsize of nfields_unit
nfields_max = 200
nfields_min = 2
nfields_unit = 1


#Range of betas from the Marcenko-Pastur distribution to iterate over
#Will create an array between min-max with numb of grid points
beta_ratio_max = 0.6
beta_ratio_min = 0.4
beta_ratio_numb = 10
#beta_ratio_max = 0.5
#beta_ratio_min = 0.5
#beta_ratio_numb = 1

#<m_avg^2> = sigma^2 for GRM w/entries of std sigma
m_avg = 5e-7


#Number of sample points to get for each set of hyperparameters
#nsamples=2e7
nsamples=200

#Should we get less samples with more fields?
scale_nsamples = False

#Output file name "root"
#Will create file called root#.dat where #=mpi_rank
fileroot = "data/outdata"
