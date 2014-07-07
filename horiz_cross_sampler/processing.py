#!/usr/bin/python

"""Module for functions that are used to process a sample in order to build things like PDF estimators."""

def build_pdf_hist(sample):
    """Estimates the probability density function (PDF) of a multidimensional sample by using a histogram.

    The sample is expected to be in the form of a list of dictionaries, where each dictionary corresponds to one sample and has keys that are the observables' names and values that are their value.  This will probably not work if any two dictionaries in sample report different observables."""

    print "boob"
