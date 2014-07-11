#!/usr/bin/python

"""Program to analyze the data output from driver.py."""

import numpy as np
import cPickle
import sys
import argparse

import processing

def parse_commandline():
    """Parse command line arguments to find the data files."""

    data_files = None

    if not sys.argv[1:]:
        #Empty list.
        raise TypeError("Please specify the parameter files at the command line " \
                "with the -d flag")
    else:
        parser = argparse.ArgumentParser(description="Analyze data that was output "\
                "from the horizon crossing sampler driver.py (cPickle dump). "\
                "Give the data files on the command line. ")
        parser.add_argument("-d", help="Data file(s) in the format (cPickle dump) of driver.py. "\
                "Default = fileroot.dat for fileroot from parameter.py.",
                nargs="+")

        data_files= parser.parse_args().d

    return data_files


def main():
    """Program that loads the data, which was output from the horizon crossing sampler driver.py, from a file and analyzes it."""

    #Get the data file names from the command line (or default).
    data_files = parse_commandline()
    if data_files==None:
        raise Exception("Error in processing data_files given on the command line.")


    #Load the histogram output and run data
    data = []
    print "Reading from data files..."
    def readin(fname):
        #print "Reading from %s" %files
        myfile = open(fname,'r')
        data.append(cPickle.load(myfile))
        myfile.close()

    map(readin, data_files)

    #Combine data to one list
    #data = [inner for outer in data for inner in outer]

    #Find observables and parameters that were iterated over
    #Get default observable and (hyper)parameter names from parameters.py
    #Which observables do we care about?
    aux_params = ['sample']
    fixed_params =['dimn_weight','m_avg']
    observs_tostudy = ['n_s']

    observs = data[0]['observs']
    var_params = [key for key in data[0].keys()
            if key not in ['observs'] + aux_params + fixed_params]


    print "%s are the observables that were calculated." %observs
    print "%s are the observables that we will study here." %observs_tostudy
    print "%s are the parameters that were iterated over." %var_params
    print "%s are fixed." %fixed_params
    #print "%s are auxiliary." %aux_params

    for obs in observs_tostudy:
        if obs not in observs:
            #observs_tostudy doesn't match any of the computed observables
            raise TypeError("The observable we are studying %s does not match "\
                    "any of the observables in what were computed: %s." \
                    %(obs,observs))


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

    params_nomarg = [params for params in var_params
            if params not in params_to_marginalize]

    print "To marginalize:", params_to_marginalize
    print "To plot:", params_nomarg

    if not params_nomarg:
        raise Exception("You have asked to marginalize over all parameters, "\
                "which isn't supported.")

    #Collapse the dimensions of fixed params or margin. params from each sample
    for sample in data:
        map(sample.pop, fixed_params + params_to_marginalize + ['observs'])

    #Make combined datasets based off non-marg. and non-fixed params
    #Results in dictionary with:
    #    key=tuple(unique non-marg parameters), value=(total sample)
    #but we remove the unwanted observables from the total sample
    proc_data = {}
    for sample in data:
        key = tuple([sample[p] for p in params_nomarg])
        temp_samp = sample.pop('sample')
        temp_samp = [ {obs:x[obs]  for obs in observs_tostudy} for x in temp_samp]
        if key in proc_data:
            proc_data[key] += temp_samp
        else:
            proc_data[key] = temp_samp

        map(sample.pop, sample.keys())


    #Make histograms for each observable's PDF versus each hyperparameter

    #Use to force the histogram to give same number of bins over some pre-defined
    #region in observable space
    fixed_bins=True
    obs_range = {'n_s': [0.88, 0.965],
            'alpha_s': [-1.0e-2,-1.0e-3],
            'f_NL': [-1.0e-2,-5.0e-3],
            'r': [0e0,0.5e0]
            }
    fixed_range = [obs_range[obs] for obs in sorted(observs_tostudy)] #Sort bc in alphab order later
    nbins = 5

    hist_total = processing.make_histograms(data=proc_data, params=params_nomarg,
            fixed_bins=fixed_bins, normed=False,
            fixed_range=fixed_range, nbins=nbins)

    for i in hist_total:
        print i['counts']
        print i['edges']


if __name__=="__main__":
    main()
