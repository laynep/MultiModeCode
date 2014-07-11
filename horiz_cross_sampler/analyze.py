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
    for files in data_files:
        #print "Reading from %s" %files
        myfile = open(files,'r')
        data.append(cPickle.load(myfile))
        myfile.close()
    #Combine data to one list
    #data = [inner for outer in data for inner in outer]

    #Find observables and parameters that were iterated over
    #Get default observable and (hyper)parameter names from parameters.py
    aux_params = ['sample']
    fixed_params =['dimn_weight','m_avg']
    observs = data[0]['observs']
    var_params = [key for key in data[0].keys()
            if key not in ['observs'] + aux_params + fixed_params]


    print "Params %s are the observables that were calculated." %observs
    print "Params %s are the parameters that were iterated over." %var_params
    print "Params %s are fixed." %fixed_params
    print "Params %s are auxiliary." %aux_params


    #To marginalize or not to marginalize?
    #    marginalize=True ==> If there is more than one hyperparameter, then
    #       compress the marginalized data and overplot everything else.
    #       Make one plot for each combination of observables and hyperparameters.
    #    marginalize=False ==> If there is more than one hyperparameter, then
    #       make a new plot for each iteration of the other hyperparameters.
    #       Make one plot for each individual data run.

    params_to_marginalize = ['beta']

    params_nomarg = [params for params in var_params
            if params not in params_to_marginalize]

    print "To marginalize:", params_to_marginalize
    print "To plot:", params_nomarg

    #Collapse the dimensions of fixed params or margin. params from each sample
    for sample in data:
        map(sample.pop, fixed_params + params_to_marginalize)

    #Make combined datasets based off non-marg. and non-fixed params
    for params in params_nomarg :
        for sample in data:
            #test = [sample['sample']
            sys.exit()



    #Plot histograms for each observable's PDF versus each hyperparameter

    #Make the histograms

    #Use to force the histogram to give same number of bins over some pre-defined
    #region in observable space
    fixed_bins=True
    obs_range = {'n_s': [0.88, 0.965],
            'alpha_s': [-1.0e-2,-1.0e-3],
            'f_NL': [-1.0e-2,-5.0e-3],
            'r': [0e0,0.5e0]
            }
    fixed_range = [obs_range[obs] for obs in sorted(observs)] #Sort bc in alphab order later
    nbins = 20

    #if p1.fixed_bins:
    #    hist_total[-1]['counts'], hist_total[-1]['edges'] = \
    #            processing.hist_estimate_pdf(sample,normed=False,
    #                    datarange=p1.fixed_range, nbins=p1.nbins)
    #else:
    #    hist_total[-1]['counts'], hist_total[-1]['edges'] = \
    #            processing.hist_estimate_pdf(sample,normed=False,
    #                    bin_method=processing.scott_rule)



if __name__=="__main__":
    main()
