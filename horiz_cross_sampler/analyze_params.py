"""Parameters that would change between runs of analysis program analyze.py"""

#Which observables do we care about?
observs_tostudy = ['n_s']

#Which parameters didn't vary
fixed_params =['dimn_weight','m_avg']

#Things that aren't parameters or observables
aux_params = ['sample']

#To marginalize or not to marginalize?
#    marginalize=True ==> If there is more than one hyperparameter, then
#       compress the marginalized data and overplot everything else.
#       Make one plot for each combination of observables and hyperparameters.
#    marginalize=False ==> If there is more than one hyperparameter, then
#       make a new plot for each iteration of the other hyperparameters.
#       Make one plot for each individual data run.

#params_to_marginalize = []
params_to_marginalize = ['beta']
#params_to_marginalize = ['nfields']

#Use to force the histogram to give same number of bins over some pre-defined
#region in observable space
fixed_bins=True
obs_range = {'n_s': [0.88, 0.965],
        'alpha_s': [-1.0e-2,-1.0e-3],
        'f_NL': [-1.0e-2,-5.0e-3],
        'r': [0e0,0.5e0]
        }
nbins = 5
