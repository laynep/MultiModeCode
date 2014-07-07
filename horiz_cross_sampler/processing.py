#!/usr/bin/python

"""Module for functions that are used to process a sample in order to build things like PDF estimators."""

import numpy as np
import sys

def bins_scott_rule(sample):
    """Returns the bin counts per dimension when using a modified Scott binning technique.  Sample is an N-dimensional Numpy array."""

    bins = np.ceil( [ (max(colmn)-min(colmn))/((3.5/len(colmn)**(1.0/3.0))*np.std(colmn))
            for colmn in sample.T])

    return bins


def build_pdf_hist(sample, observables=None, equal_bins=False, normed=True):
    """Estimates the probability density function (PDF) of an N-dimensional sample by using a histogram.  Returns the bin counts and edges of the histogram bins.

    The sample is expected to be in the form of a list of dictionaries, where each dictionary corresponds to one sample and has keys that are the observables' names and values that are their value.  This will probably not work if any two dictionaries in sample report different observables.

    If observables is present, then it is a list of length N that describes which keys in the dictionaries of the sample are to be used in the estimation.  To determine how many bins to give the histogram the default is (1) Bayesian Blocks and (2) the Scott rule, if there is no astroML module.  Bayesian Blocks only works in 1D.  If equal_bins=True, then will not use Bayesian Blocks, which in general gives a non-uniform bin size.  Defaults to PDF normalized bin counts."""

    #Name the dimensions
    if observables==None:
        observables = sample[1].keys()

    #Focus on the dimensions that coincide with observables
    #Matches ordering of observables
    try:
        if len(observables)==1:
            sample = np.array([entries[key] for key in observables for entries in sample])
        else:
            sample = np.array([[entries[key] for key in observables] for entries in sample])

    except:
        raise TypeError("The observables %s have given an error." %observables)


    #Make the histogram
    try:
        #If the dependency exists, then use AstroML,
        #which has a neat Bayesian Blocks implementation in 1D

        if len(observables)==1:
            import astroML.density_estimation as aML

            if equal_bins:
                counts, edges = aML.histogram(sample,bins='scotts',normed=normed)
            else:
                counts, edges = aML.histogram(sample,normed=normed)
        else:
            raise TypeError()

    except:
        #Otherwise, just use the Numpy histogram feature

        bins = bins_scott_rule(sample)
        counts, edges = np.histogramdd(sample, bins=bins, normed=normed)

    return counts, edges
