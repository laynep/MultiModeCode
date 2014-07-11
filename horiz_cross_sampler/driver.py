#!/usr/bin/python

"""Program to sample the horizon crossing surface."""

import numpy as np
import cosmology as cosmo
import mpi_routines as par

import cPickle
import sys
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
                "Give the parameter files on the command line.")
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
        "nsamples", "scale_nsamples",
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


    def load_sample_dictionary():
        """Loads the hyper parameters into a dictionary for each run."""
        sample_total['beta'] = beta
        sample_total['nfields'] = nfields
        sample_total['m_avg'] = p1.m_avg
        sample_total['dimn_weight'] = dimn_weight
        sample_total['observs'] = p1.obs_to_calc


    #(Parallelized) iteration over all desired combinations of hyperparameters
    counter=0
    for dimn_weight, nfields in loop_params:
        for beta in beta_list:

            print "\nNew run:"
            print "beta=", beta
            print "nfields=", nfields

            sample_total = {}

            load_sample_dictionary()

            #How many samples to build PDF from
            if p1.scale_nsamples:
                if nfields>=10:
                    this_nsamples = int(np.ceil(p1.nsamples/nfields**2))
                else:
                    this_nsamples = int(np.ceil(p1.nsamples/10**2))
            else:
                this_nsamples = int(p1.nsamples)

            print "nsamples=", this_nsamples

            #nmoduli = axions + dilaton + heavy moduli
            #nmoduli = nfields + 1 + 1.0*nfields
            nmoduli = nfields/beta

            run = cosmo.SR_universe(sampler="MP_and_horizcross",HC_approx=True,
                    model="Nquad", nfields=nfields)
            radius = 2.0*np.sqrt(run.N_pivot)

            #Keep the sample
            sample_total['sample'] = run.sample_Nquad(p1.obs_to_calc, this_nsamples,
                    nmoduli, radius, p1.m_avg, dimn_weight)

            counter += 1

            #Write sample output to file
            if parallel:
                myfile = open(p1.fileroot+'_rank'+str(mpi_rank)+'_samp'+str(counter)+".dat","w")
            else:
                myfile = open(p1.fileroot+'_samp'+str(counter)+".dat","w")
            cPickle.dump(sample_total, myfile)
            myfile.close()




if __name__=="__main__":

    profile=False
    if not profile:
        main()
    else:
        import cProfile
        cProfile.run('main()')
