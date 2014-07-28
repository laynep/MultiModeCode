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
    fromlist = [ "obs_to_calc", "hyperparams", "sampler",
        "nsamples", "scale_nsamples",
        "fileroot"]

    p1 = __import__(param_file, fromlist=fromlist)

    #Set up grid of hyperparameters
    hgrid = {}
    for param,values in p1.hyperparams.iteritems():
        if isinstance(values[0],float):
            hgrid[param] = np.linspace(values[0], values[1], values[2])
        elif isinstance(values[0],int):
            hgrid[param] = np.arange(values[0], values[1]+1, values[2])
        else:
            raise TypeError("Hyperparameter %s does not have appropriate ranges." %param)



    if not parallel or mpi_rank==0:

        #For initial conditions:
        #Uniform weighting of dimensions for ICs
        ic_weight = [np.ones(f_numb) for f_numb in hgrid['nfields']]

        #Iterate over the rows in this list
        loop_params = zip(ic_weight, hgrid['nfields'])

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


    def get_nsamples(nsamples, nfields, scale=False):
        """How many samples to build PDF from."""
        if scale:
            if nfields>=10:
                this_nsamples = int(np.ceil(nsamples/nfields**2))
            else:
                this_nsamples = int(np.ceil(nsamples/10**2))
        else:
            this_nsamples = int(nsamples)

        print "nsamples=", this_nsamples
        return this_nsamples

    def write_to_file():
        """Write sample output to file."""

        if parallel:
            myfile = open(p1.fileroot+'_rank'+str(mpi_rank)+'_samp'+str(counter)+".dat","w")
        else:
            myfile = open(p1.fileroot+'_samp'+str(counter)+".dat","w")
        cPickle.dump(sample_total, myfile)
        myfile.close()


    #(Parallelized) iteration over all desired combinations of hyperparameters
    counter=0
    for dimn_weight, nfields in loop_params:


        if p1.sampler == "MP":
            #Use HCA, sample horiz cross surf uniformly, use Marcenko-Pastur for couplings
            #Designed for p=2, ie, N-quadratic

            def load_sample_dictionary():
                """Loads the hyper parameters into a dictionary for each run."""
                sample_total['observs'] = p1.obs_to_calc
                sample_total['beta'] = beta
                sample_total['nfields'] = nfields
                sample_total['m_avg'] = m_avg
                sample_total['dimn_weight'] = dimn_weight
                sample_total['p'] = p

            for beta, m_avg, p in [ [beta, m_avg, p] for beta in hgrid['beta']
                    for m_avg in hgrid['m_avg'] for p in hgrid['p']]:


                print "\nNew run:"
                print "beta=", beta
                print "nfields=", nfields
                print "p=", p

                sample_total = {}

                load_sample_dictionary()

                this_nsamples = get_nsamples(p1.nsamples, nfields, p1.scale_nsamples)

                #nmoduli = axions + dilaton + heavy moduli
                #nmoduli = nfields + 1 + 1.0*nfields
                nmoduli = nfields/beta

                run = cosmo.Nmono_universe(sampler="MP_and_horizcross",HC_approx=True,
                        model="Nmono", nfields=nfields)
                radius = np.sqrt(2.0*p*run.N_pivot)

                sample_total['sample'] = run.sample_Nmono_MP(p1.obs_to_calc, this_nsamples,
                        nmoduli, radius, m_avg, dimn_weight, p)


                counter += 1

                write_to_file()

        elif p1.sampler == "uniform" or p1.sampler=="log":
            #Use HCA, sample the horiz cross surf uniformly
            #Use the uniform distrib for either the couplings or the logs

            def load_sample_dictionary():
                """Loads the hyper parameters into a dictionary for each run."""
                sample_total['observs'] = p1.obs_to_calc
                sample_total['low'] = low
                sample_total['high'] = high
                sample_total['nfields'] = nfields
                sample_total['dimn_weight'] = dimn_weight
                sample_total['p'] = p

            for low, high, p in [ [low, high, p] for low in hgrid['low']
                    for high in hgrid['high'] for p in hgrid['p']]:


                print "\nNew run:"
                print "low=", low
                print "high=", high
                print "nfields=", nfields
                print "p=", p

                sample_total = {}

                load_sample_dictionary()

                this_nsamples = get_nsamples(p1.nsamples, nfields, p1.scale_nsamples)

                run = cosmo.Nmono_universe(sampler=p1.sampler,HC_approx=True,
                        model="Nmono", nfields=nfields)
                radius = np.sqrt(2.0*p*run.N_pivot)

                sample_total['sample'] = run.sample_Nmono_uniform(p1.obs_to_calc, this_nsamples,
                        low, high, radius, dimn_weight, p)


                counter += 1

                write_to_file()


        else:
            raise TypeError("The sampling technique %s is not defined." %p1.sampler)






if __name__=="__main__":

    profile=False
    if not profile:
        main()
    else:
        import cProfile
        cProfile.run('main()')
