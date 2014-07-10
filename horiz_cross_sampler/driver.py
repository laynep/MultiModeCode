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
    """Parse command line arguments to find the parameter files.  argv should be sys.argv[1:]."""

    param_files_threads = None
    param_files_master = None

    if not sys.argv[1:]:
        #Empty list.  Use default parameter file names.
        param_files_threads = "params_allthreads"
        param_files_master = "params_master"
    else:
        parser = argparse.ArgumentParser(description="Parameter files "\
                "to find the all threads variables " \
                "and the master-only variables.")
        parser.add_argument("-p1", help="Parameter file "\
                "that is seen by all the threads.")
        parser.add_argument("-p2", help="Parameter file "\
                "that is seen by only the master thread.")

        def remv_last3(string):
            #Take off the ".py" ending
            return string.replace(' ', '')[:-3]

        param_files_threads= remv_last3(parser.parse_args().p1)
        param_files_master = remv_last3(parser.parse_args().p2)

    return param_files_threads, param_files_master


def main():
    """Driver function for sampling N-quadratic inflation."""

    #Find parameter files
    param_files_threads, param_files_master = parse_commandline()

    #MPI parallelized
    mpi_comm, mpi_size, mpi_rank, mpi_name = par.init_parallel()
    parallel = not mpi_comm==None

    #Get the run parameters that all the threads see
    fromlist = [ "obs_to_calc", "hyperparams",
        "beta_ratio_max", "beta_ratio_min", "beta_ratio_numb",
        "m_avg",
        "fixed_bins", "obs_range", "fixed_range", "nbins",
        "samp_ref", "scale_nsamples",
        "fileroot"]

    p1 = __import__(param_files_threads, fromlist=fromlist)


    #Set up grid of hyperparameters

    beta_list = np.linspace(p1.beta_ratio_min, p1.beta_ratio_max, p1.beta_ratio_numb)


    if not parallel or mpi_rank==0:

        #Get the run parameters that only the master thread needs to see
        fromlist = ["nfields_max", "nfields_min", "nfields_unit"]
        try:
            p2 = __import__(param_files_master, fromlist=fromlist)
        except:
            raise ImportError("TEST ERROR2")

        nfields_list = np.arange(p2.nfields_min,p2.nfields_max+1,p2.nfields_unit)

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
                                datarange=fixed_range, nbins=nbins)
            else:
                hist_total[-1]['counts'], hist_total[-1]['edges'] = \
                        processing.hist_estimate_pdf(sample,normed=False,
                                bin_method=processing.scott_rule)

    for i in hist_total:
        print "this is count:", i['counts']
        print "this is edges:", i['edges']

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
