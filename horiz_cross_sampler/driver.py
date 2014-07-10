#!/usr/bin/python

import numpy as np
import cosmology as cosmo
import processing
import mpi_routines as par

import cPickle
import sys
import getopt
import argparse


def parse_commandline():
    """Parse command line arguments to find the parameter files."""

    param_file = None

    if not sys.argv[1:]:
        #Empty list.  Use default parameter file names.
        param_file = "parameters"
    else:
        parser = argparse.ArgumentParser(description="Sample the horizon crossing surface "\
                "for an inflation model to build PDFs for given (hyper)parameters. "\
                "Give the parameter files on the command line ")
        parser.add_argument("-p", help="Parameter file "\
                "that is seen by all the threads. "\
                "Default = parameters.py")

        def remv_last3(string):
            #Take off the ".py" ending
            return string.replace(' ', '')[:-3]

        param_file= remv_last3(parser.parse_args().p)

    return param_file


def main():
    """Driver function for sampling N-quadratic inflation."""

    #MPI parallelized
    mpi_comm, mpi_size, mpi_rank, mpi_name = par.init_parallel()
    parallel = not mpi_comm==None


    #Find parameter files
    param_file = parse_commandline()


    #Get the run parameters that all the threads see
    fromlist = [ "obs_to_calc", "hyperparams",
        "beta_ratio_max", "beta_ratio_min", "beta_ratio_numb",
        "m_avg",
        "fixed_bins", "obs_range", "fixed_range", "nbins",
        "samp_ref", "scale_nsamples",
        "fileroot",
        "nfields_max", "nfields_min", "nfields_unit"]

    p1 = __import__(param_file, fromlist=fromlist)

    #Set up grid of hyperparameters

    beta_list = np.linspace(p1.beta_ratio_min, p1.beta_ratio_max, p1.beta_ratio_numb)


    if not parallel or mpi_rank==0:

        nfields_list = np.arange(p1.nfields_min,p1.nfields_max+1,p1.nfields_unit)

        #For initial conditions:
        #Uniform weighting of dimensions for ICs
        ic_weight = [np.ones(f_numb) for f_numb in nfields_list]

        #Iterate over the rows in this list
        loop_params = zip(ic_weight, nfields_list)

    else:

        loop_params = None


    if parallel and mpi_size>1:
        if mpi_rank==0:
            #Chunk the loop_params into pieces so each process can loop over a subset.
            loop_params = par.chunk(np.array(loop_params),mpi_size,group=False)

        #Halt processors until run parameters are chunked.
        mpi_comm.barrier()

        #Scatter loop_params to all processes
        loop_params = mpi_comm.scatter(loop_params,root=0)


    def load_hist_dictionary():
        """Loads the hyper parameters into a dictionary for each run."""
        hist_total.append({})
        hist_total[-1]['beta'] = beta
        hist_total[-1]['nfields'] = nfields
        hist_total[-1]['m_avg'] = p1.m_avg
        hist_total[-1]['dimn_weight'] = dimn_weight
        hist_total[-1]['observs'] = p1.obs_to_calc


    #(Parallelized) iteration over all desired combinations of hyperparameters
    hist_total = []
    for dimn_weight, nfields in loop_params:
        for beta in beta_list:
            print "\nNew run:"
            print "beta=", beta
            print "nfields=", nfields

            load_hist_dictionary()

            #How many samples to build PDF from
            if p1.scale_nsamples:
                if nfields>=10:
                    nsamples = int(np.ceil(p1.samp_ref/nfields**2))
                else:
                    nsamples = int(np.ceil(p1.samp_ref/10**2))
            else:
                nsamples = int(p1.samp_ref)

            print "nsamples=", nsamples

            #nmoduli = axions + dilaton + heavy moduli
            #nmoduli = nfields + 1 + 1.0*nfields
            nmoduli = nfields/beta

            run = cosmo.SR_universe(sampler="MP_and_horizcross",HC_approx=True,
                    model="Nquad", nfields=nfields)
            radius = 2.0*np.sqrt(run.N_pivot)

            sample = run.sample_Nquad(p1.obs_to_calc, nsamples, nmoduli, radius, p1.m_avg, dimn_weight)

            if p1.fixed_bins:
                hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                        processing.hist_estimate_pdf(sample,normed=False,
                                datarange=p1.fixed_range, nbins=p1.nbins)
            else:
                hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                        processing.hist_estimate_pdf(sample,normed=False,
                                bin_method=processing.scott_rule)

    #for i in hist_total:
    #    print "this is count:", i['counts']
    #    print "this is edges:", i['edges']

    #Write histogram output to file
    if parallel:
        myfile = open(p1.fileroot+str(mpi_rank)+".dat","w")
    else:
        myfile = open(p1.fileroot+".dat","w")
    cPickle.dump(hist_total, myfile)
    myfile.close()




if __name__=="__main__":

    profile=False
    if not profile:
        main()
    else:
        import cProfile
        cProfile.run('main()')
