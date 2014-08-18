#!/usr/bin/python

"""Module for functions that are used to process a sample in order to build things like PDF estimators."""

import numpy as np
import sys

def scott_rule(sample):
    """Returns the bin positions per dimension when using a modified Scott binning technique.  The sample is an N-dimensional Numpy array."""

    if len(sample.shape) >1 :
        nbins = np.ceil(  (max(colmn)-min(colmn))/((3.5/len(colmn)**(1.0/3.0))*np.std(colmn))
                for colmn in sample.T)
    elif len(sample.shape) ==1 :
        nbins = np.ceil(  (max(sample)-min(sample))/((3.5/len(sample)**(1.0/3.0))*np.std(sample)) )
    else:
        raise TypeError("Something wrong with the sample, since sample.shape = %s." %sample.shape)

    return nbins


def hist_estimate_pdf(sample, preproccess=False, observables=None, normed=True,
        nbins=None, bin_method=scott_rule, datarange=None):
    """Estimates the probability density function (PDF) of an N-dimensional sample by using a histogram.  Returns the bin counts and edges of the histogram bins.

    The sample is expected to be in the form of a list of dictionaries, where each dictionary corresponds to one sample and has keys that are the observables' names and values that are their value.  This will probably not work if any two dictionaries in sample report different observables.

    If observables is present, then it is a list of length N that describes which keys in the dictionaries of the sample are to be used in the estimation.  Use nbins to force the histogram to use this many bins; to use an automated "principled" approach, make sure nbins=None.  Set bin_method to determine how many bins to give the histogram; the default is scott_rule.  To give histogramdd a range, specify datarange for each dimension."""

    if preproccess==False:

        if observables==None:
            #Name the dimensions.  Alphabetize to maintain some ordering.
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

def make_histograms(data, params,
        fixed_bins=False, normed=False, fixed_range=None, nbins=None):
    """Routine to return histogram bin counts and edges from a dataset data that is expressed as a dictionary with unique keys that describe unique parameter choices and values corresponding to a sample made by driver.py."""

    hist_total=[]
    for key in data:
        hist_total.append({})

        for entry, value in zip(params,key):
            hist_total[-1][entry]=value

        if fixed_bins:
            if nbins==None:
                raise TypeError("Set nbins when calling make_histograms")
            elif fixed_range==None:
                raise TypeError("Set fixed_range when calling make_histograms")

            hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                    hist_estimate_pdf(data[key],
                            normed=normed,
                            datarange=fixed_range,
                            nbins=nbins)
        else:
            hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                    hist_estimate_pdf(data[key],
                            normed=normed,
                            bin_method=scott_rule)
    return hist_total

def arrange_for_imshow(hist_total):
    """Routine to arrange the histogram data into 2D histograms based off the bin count and edge results from make_histograms().  Plots the observable vs each hyperparameter for the range over which both are defined.  Needs fixed bin sizes in both dimensions (observs and params).  Expects a list of dictionaries hist_total that have keys that correspond to unique values of hyperparameter combinations and 'counts' and 'edges'."""


    #Check the bin counts are in the correct form
    hist_params = ['edges', 'counts']
    params = [keys for keys in hist_total[0].keys() if keys not in hist_params]

    if len(hist_total[0]['counts'].shape) > 1:
        raise TypeError("To use plot_histogram, the bin counts need to be "\
                "one-dimensional for each set of unique parameters.")
    elif len(params)>1:
        raise TypeError("Plotting more than one non-marginalized parameter is" \
                "not currently implemented in plot_histogram.")


    #Plotting range and parameter lists
    obs_range = [hist_total[0]['edges'][0][0],hist_total[0]['edges'][-1][-1]]
    extent = [min( hist[p] for p in params for hist in hist_total),
                max( hist[p] for p in params for hist in hist_total),
                obs_range[0], obs_range[1]]

    #Make sure this is in numerical order for the parameter
    #params_data = np.array([[counts for counts in elts['counts']] for elts in hist_total])

    for p in params:
        params_data = np.array([counts for name, counts in
            sorted([elts[p],elts['counts']] for elts in hist_total)])


    return params_data, extent


