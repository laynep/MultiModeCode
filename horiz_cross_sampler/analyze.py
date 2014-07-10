#!/usr/bin/python

"""Program to analyze the data output from driver.py."""

import numpy as np
import cPickle
import sys
import argparse

def parse_commandline():
    """Parse command line arguments to find the data files."""

    data_files = None

    if not sys.argv[1:]:
        #Empty list.  Use default parameter file names.
        try:
            params = __import__("parameters",fromlist=['fileroot'])
            data_files = [params.fileroot+".dat"]
        except:
            raise ImportError("Couldn't find the default parameter file " \
            "parameters.py or the fileroot parameter isn't defined.")
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
    for files in data_files:
        print "Reading from %s" %files
        myfile = open(files,'r')
        data.append(cPickle.load(myfile))
        myfile.close()
    #Combine data to one list
    data = [inner for outer in data for inner in outer]

    #Find observables and parameters that were iterated over
    #Get default observable and (hyper)parameter names from parameters.py
    hist_params = ['counts','edges']
    observs = data[0]['observs']
    hyperparams = [key for key in data[0].keys()
            if key not in hist_params + ['observs',
                                            'dimn_weight','m_avg'] ] #Ad hoc


    print "These %s are the observables that were calculated." %observs
    print "These %s are the parameters that were iterated over." %hyperparams


    #Plot histograms for each observable's PDF versus each hyperparameter

    #To marginalize or not to marginalize?
    #    marginalize=True ==> If there is more than one hyperparameter, then
    #       compress the marginalized data and overplot everything else.
    #       Make one plot for each combination of observables and hyperparameters.
    #    marginalize=False ==> If there is more than one hyperparameter, then
    #       make a new plot for each iteration of the other hyperparameters.
    #       Make one plot for each individual data run.

    params_to_marginalize = ['beta']
    marginalize = {}
    for param in hyperparams:
        marginalize[param] = param in params_to_marginalize

    print marginalize




if __name__=="__main__":
    main()
