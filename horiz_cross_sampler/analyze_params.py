"""Parameters that would change between runs of analysis program analyze.py"""

#Which observable do we care about?
observs_tostudy = ['n_t']
#observs_tostudy = ['PR']

#Which parameters didn't vary
#fixed_params =['dimn_weight','m_avg', 'p']
#fixed_params =['dimn_weight','m_avg','beta','p']

fixed_params =['dimn_weight','low', 'high', 'p']

#Things that aren't parameters or observables
#You probably don't need to change this.
aux_params = ['sample']

#To marginalize or not to marginalize?
#    marginalize=True ==> If there is more than one hyperparameter, then
#       compress the marginalized data and overplot everything else.
#       Make one plot for each combination of observables and hyperparameters.
#    marginalize=False ==> If there is more than one hyperparameter, then
#       make a new plot for each iteration of the other hyperparameters.
#       Make one plot for each individual data run.

params_to_marginalize = []
#params_to_marginalize = ['beta']
#params_to_marginalize = ['nfields']

#Use to force the histogram to give same number of bins over some pre-defined
#region in observable space
fixed_bins=True
norm_PDF=False
obs_range = {
        'PR': [1e-10, 1e-9],
        #'n_s': [0.94, 0.965],
        'n_s': [0.93, 0.95],
        'alpha_s': [-0.005,-0.0015],
        'f_NL': [-1.0e-2,-5.0e-3],
        'r': [0e0,0.5e0],
        #'n_t': [-0.04, -0.015],
        #'n_t': [-0.15, -0.015],
        #'n_t': [-0.027, -0.017],
        'n_t': [-0.027, -0.022],
        }
nbins = 30

#Names of observables
obs_name = {'PR':r'$\mathcal{P}_\mathcal{R}$',
        'n_s':r'$n_s$',
        'alpha_s':r'$\alpha_s$',
        'r':r'$r$',
        'n_t':r'$n_T$',
        'f_NL':r'$f_\mathrm{NL}$',
        'tau_NL':r'$\tau_\mathrm{NL}$'}

#Names of parameters in plots
param_name = {'nfields':r'$N_f$',
        'beta':r'$\beta$',
        'dimn_weight':r'$\lambda$',
        'm_avg':r'$\bar{m}_\mathrm{avg}$'}
