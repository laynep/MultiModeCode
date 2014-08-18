#!/usr/bin/python

"""Program to analyze the data output from driver.py."""

import numpy as np
import cPickle
import sys
import argparse

import processing
import matplotlib.pyplot as plt
import pyplotsetup

import pandas as pd


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


def process():
    """Program that loads the data, which was output from the horizon crossing sampler driver.py, from a file and analyzes it."""

    #Get the data file names from the command line (or default).
    data_files = parse_commandline()
    if data_files==None:
        raise Exception("Error in processing data_files given on the command line.")


    #Load the histogram output and run data
    data = []
    print "Reading from data files..."
    for fname in data_files:
        #print "Reading from %s" %files
        myfile = open(fname,'r')
        data.append(cPickle.load(myfile))
        myfile.close()


    output = {val:[] for val in data[0].keys()}
    for row in data:
        for val in output:
            output[val].append(row[val])


    for obs in output['observs'][0]:
        output[obs]=[]


    print output['observs'][0]
    print output.keys()
    for samp in output['sample']:
        for item in samp:
            for obs in output['observs'][0]:
                #print item[obs]
                #print obs
                output[obs].append(item[obs])
    #print output['sample'][0][0]
    #print output['n_t']

    #output = pd.DataFrame(output, columns=output.keys())

    #output.save('simple_output'+str(nfields)+'.pkl')



    sys.exit()


    #Find observables and parameters that were iterated over
    #Get default observable and (hyper)parameter names from analyze_params.py

    fromlist = [ "observs_tostudy", "fixed_params", "aux_params", \
            "params_to_marginalize", \
            "fixed_bins","norm_PDF", "obs_range", "nbins", \
            "obs_name","param_name"]
    p1 = __import__("analyze_params", fromlist=fromlist)

    observs = data[0]['observs']
    var_params = [key for key in data[0].keys()
            if key not in ['observs'] + p1.aux_params + p1.fixed_params]

    print "%s are the observables that were calculated." %observs
    print "%s are the observables that we will study here." %p1.observs_tostudy
    print "%s are the parameters that were iterated over." %var_params
    print "%s are fixed." %p1.fixed_params
    #print "%s are auxiliary." %aux_params

    for params in p1.fixed_params:
        if params in p1.params_to_marginalize:
            raise ImportError("Parameter %s is declared to be both " \
                    "fixed and a parameter that should be marginalized. "\
                    "Please fix." %params)

    for obs in p1.observs_tostudy:
        if obs not in observs:
            #observs_tostudy doesn't match any of the computed observables
            raise TypeError("The observable we are studying %s does not match "\
                    "any of the observables in what were computed: %s." \
                    %(obs,observs))

    #Set up parameters that aren't being marginalized over
    params_nomarg = [params for params in var_params
            if params not in p1.params_to_marginalize]

    print "To marginalize:", p1.params_to_marginalize
    print "To plot:", params_nomarg

    if not params_nomarg:
        raise Exception("You have asked to marginalize over all parameters, "\
                "which isn't supported.")

    #Collapse the dimensions of fixed params or margin. params from each sample
    for sample in data:
        map(sample.pop, p1.fixed_params + p1.params_to_marginalize + ['observs'])



    #Make combined datasets based off non-marg. and non-fixed params
    #Results in dictionary with:
    #    key=tuple(unique non-marg parameters), value=(total sample)
    #but we remove the unwanted observables from the total sample
    proc_data = {}
    for sample in data:
        key = tuple([sample[p] for p in params_nomarg])
        temp_samp = sample.pop('sample')
        temp_samp = [ {obs:x[obs]  for obs in p1.observs_tostudy} for x in temp_samp]
        if key in proc_data:
            proc_data[key] += temp_samp
        else:
            proc_data[key] = temp_samp

        map(sample.pop, sample.keys())




    #Make histograms for each observable's PDF versus each hyperparameter


    #Sort bc in alphab order later
    fixed_range = [p1.obs_range[obs] for obs in sorted(p1.observs_tostudy)]
    hist_total = processing.make_histograms(data=proc_data, params=params_nomarg,
            fixed_bins=p1.fixed_bins, normed=p1.norm_PDF,
            fixed_range=fixed_range, nbins=p1.nbins)

    #Make the plots

    plot_data, extent = processing.arrange_for_imshow(hist_total)


    #Figure options

    fig = plt.figure(**pyplotsetup.figprops)
    fig.subplots_adjust(**pyplotsetup.adjustprops)

    color_map = plt.get_cmap('RdYlBu_r')
    #color_map = plt.get_cmap('Blues')

    ax1 = fig.add_subplot('111')

    #Prediction
    N_piv=55.0
    r=8.0/N_piv

    #MP
    #x = np.linspace(2,200,100)
    #y= (-1.0/8.0)*np.ones(100)*0.14545
    #y2= (-1.0/8.0)*(1.0+0.5)*np.ones(100)*0.14545
    #sf_pred, = plt.plot(x,y, 'g--', linewidth=3.0, alpha=0.75)
    #mf_pred, = plt.plot(x,y2, 'k-', linewidth=3.0, alpha=0.95)

    #Log
    x = np.linspace(2,200,100)
    y= (-1.0/8.0)*np.ones(100)*0.14545
    b = -12.0
    a = -14.0
    c = (b-a)/np.log(10.0)
    y2= (-1.0*c/16.0)*((10.0**b + 10.0**a)/(10.0**b-10.0**a))*np.ones(100)*0.14545
    sf_pred, = plt.plot(x,y, 'g--', linewidth=3.0, alpha=0.75)
    mf_pred, = plt.plot(x,y2, 'k-', linewidth=3.0, alpha=0.95)

    #Unif
    #x = np.linspace(2,10000,100)
    #y= (-1.0/8.0)*np.ones(100)*r
    #b = 1.0e-13
    ##a = 1.0e-14
    #a = 1.0e-16
    #y2= (-1.0/6.0)*((b**3-a**3)/(b-a)/(b+a)**2)*np.ones(100)*r
    #sf_pred, = plt.plot(x,y, 'g--', linewidth=3.0, alpha=0.75)
    #mf_pred, = plt.plot(x,y2, 'k-', linewidth=3.0, alpha=0.95)

    import scipy.special as sp
    #sig = r*0.148/np.sqrt(x)
    sig = r*0.11/np.sqrt(x)

    #pval = 1.0 - 0.0455 #2 sig
    pval = 1.0 - 0.0027 #3 sig
    #pval = 1.0 - 6e-5 #4 sig
    a = np.sqrt(2.0)*sig*sp.erfinv(pval)
    y3= y2+a
    y4= y2-a
    siglow_pred, = plt.plot(x,y3, 'k--', linewidth=1.0, alpha=0.75)
    sighigh_pred, = plt.plot(x,y4,'k--', linewidth=1.0, alpha=0.75)

    import seaborn as sns

    sns.set_style('darkgrid')



    plt.rc('legend',**{'fontsize':8})

    #plt.legend([sf_pred,mf_pred],[r'$N_f=1$', r'$N_f \to \infty$'],
    #        loc='upper right')
    plt.legend([mf_pred],[r'$N_f \to \infty$'],
            loc='upper right')


    #plt.imshow(plot_data.T, origin='lower', extent=extent, aspect='auto',
    #        cmap=color_map, interpolation='nearest')

    plt.imshow(plot_data.T, origin='lower', extent=extent, aspect='auto')

    ax1.set_xlabel(r'$N_f$')
    ax1.set_ylabel(r'$n_t$')

    plt.xlim(2,1100)


    save=True
    if save:
        name = 'manyfield_log'
        direct='/home/lpri691/LaTex/multifield_modecode/consistency_rel/plots/'
        plt.savefig(direct+name+'.pdf')



    plt.show()


if __name__=="__main__":
    process()
