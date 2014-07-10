#!/usr/bin/python

"""Module for functions that are used to process a sample in order to build things like PDF estimators."""

import numpy as np
import sys

def scott_rule(sample):
    """Returns the bin positions per dimension when using a modified Scott binning technique.  The sample is an N-dimensional Numpy array."""

    if len(sample.shape) >1 :
        nbins = np.ceil( [ (max(colmn)-min(colmn))/((3.5/len(colmn)**(1.0/3.0))*np.std(colmn))
                for colmn in sample.T])
    elif len(sample.shape) ==1 :
        nbins = np.ceil( [ (max(sample)-min(sample))/((3.5/len(sample)**(1.0/3.0))*np.std(sample)) ])
    else:
        raise TypeError("Something wrong with the sample, since sample.shape = %s." %sample.shape)

    return nbins


def hist_estimate_pdf(sample, observables=None, normed=True,
        nbins=None, bin_method=scott_rule, datarange=None):
    """Estimates the probability density function (PDF) of an N-dimensional sample by using a histogram.  Returns the bin counts and edges of the histogram bins.

    The sample is expected to be in the form of a list of dictionaries, where each dictionary corresponds to one sample and has keys that are the observables' names and values that are their value.  This will probably not work if any two dictionaries in sample report different observables.

    If observables is present, then it is a list of length N that describes which keys in the dictionaries of the sample are to be used in the estimation.  Use nbins to force the histogram to use this many bins; to use an automated "principled" approach, make sure nbins=None.  Set bin_method to determine how many bins to give the histogram; the default is scott_rule.  To give histogramdd a range, specify datarange for each dimension."""

    #Name the dimensions
    if observables==None:
        observables = sorted(sample[1].keys())

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
    #try:
    #    #If the dependency exists, then use AstroML
    #    if len(observables)==1:
    #        import astroML.density_estimation as aML

    #        counts, edges = aML.histogram(sample,bins=bin_method,normed=normed)
    #    else:
    #        raise TypeError()


    #Just use the damn Numpy histogram
    if nbins==None: nbins = bin_method(sample)

    if datarange==None:
        counts, edges = np.histogramdd(sample, bins=nbins, normed=normed)
    else:
        counts, edges = np.histogramdd(sample, bins=nbins, normed=normed, range=datarange)

    return counts, edges
