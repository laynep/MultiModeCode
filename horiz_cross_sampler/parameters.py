#!/usr/bin/python

"""Parameters that would change between runs of the horizon crossing sampler.  Variables that should be available for all the threads.  For serial runs or only one MPI process, all parameters are read by the same thread."""

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


#Range of field number to iterate over
#Will create an array between min-max with stepsize of nfields_unit
nfields_max = 100
nfields_min = 50
nfields_unit = 10


#Range of betas from the Marcenko-Pastur distribution to iterate over
#Will create an array between min-max with numb of grid points
beta_ratio_max = 0.6
beta_ratio_min = 0.4
beta_ratio_numb = 10

#<m_avg^2> = sigma^2 for GRM w/entries of std sigma
m_avg = 5e-7

#Use to force the histogram to give same number of bins over some pre-defined
#region in observable space
fixed_bins=True
obs_range = {'n_s': [0.88, 0.965],
        'alpha_s': [-1.0e-2,-1.0e-3],
        'f_NL': [-1.0e-2,-5.0e-3],
        'r': [0e0,0.5e0]
        }
fixed_range = [obs_range[obs] for obs in sorted(obs_to_calc)] #Sort bc in alphab order later
nbins = 20

#Number of sample points to get for each set of hyperparameters
#nsamples=2e7
nsamples=1e4

#Should we get less samples with more fields?
scale_nsamples = False

#Output file name "root"
#Will create file called root#.dat where #=mpi_rank
fileroot = "outdata"
