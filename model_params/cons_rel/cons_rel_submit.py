#!/bin/python

import subprocess
import sys
import numpy as np

p_list = [2.0, 1.0, 1.5]
prior_list = ["unif", "log"]
first_nf_list = [20, 60, 100]
second_nf_list = [100, 250, 500, 750, 1000]
#approx_list = ["deltaN", "HCA"]
approx_list = ["deltaN"]

range_unif = ["1.0e-14", "1.0e-13"]
range_log = ["-14", "-12"]

njobs = 2

def run_cmd(cmd):
    proc = subprocess.Popen(cmd,
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    proc.wait()

def modify_pfile(nf, prior, approx, p):
    #Modify the parameter file
    pfile = "parameters_multimodecode.txt"

    #Number of fields
    sed_cmd = "sed -i \"s/num_inflaton =.*/num_inflaton = %s/g\" " %nf +pfile
    #print sed_cmd
    run_cmd(sed_cmd)

    #Approximation scheme
    if approx is "HCA":
        sed_cmd = "sed -i \"s/use_horizon_cross_approx =.*/" \
            "use_horizon_cross_approx = %s/g\" " %(".true.") +pfile
    else:
        sed_cmd = "sed -i \"s/use_horizon_cross_approx =.*/" \
            "use_horizon_cross_approx = %s/g\" " %(".false.") +pfile
    #print sed_cmd
    run_cmd(sed_cmd)

    if prior is "unif":
        #Unif sampling

        sed_cmd = "sed -i \"s/param_sampling =.*/"\
                "param_sampling = %s/g\" " %(2) +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        #p exponent
        sed_cmd = "sed -i \"s/vp_prior_min(2,:) =.*/"\
                "vp_prior_min(2,:) = %s/g\" " %p +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        sed_cmd = "sed -i \"s/vp_prior_max(2,:) =.*/"\
                "vp_prior_max(2,:) = %s/g\" " %p +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        #Param priors
        sed_cmd = "sed -i \"s/vp_prior_min(1,:) =.*/"\
                "vp_prior_min(1,:) = %s/g\" " %range_unif[0] +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        sed_cmd = "sed -i \"s/vp_prior_max(1,:) =.*/"\
                "vp_prior_max(1,:) = %s/g\" " %range_unif[1] +pfile

        #print sed_cmd
        run_cmd(sed_cmd)


    elif prior is "log":
        #Log sampling

        sed_cmd = "sed -i \"s/param_sampling =.*/"\
                "param_sampling = %s/g\" " %(3) +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        #p exponent
        sed_cmd = "sed -i \"s/vp_prior_min(2,:) =.*/"\
                "vp_prior_min(2,:) = %s/g\" " %(np.log10(p)) +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        sed_cmd = "sed -i \"s/vp_prior_max(2,:) =.*/"\
                "vp_prior_max(2,:) = %s/g\" " %(np.log10(p)) +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        #Param priors
        sed_cmd = "sed -i \"s/vp_prior_min(1,:) =.*/"\
                "vp_prior_min(1,:) = %s/g\" " %range_log[0] +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

        sed_cmd = "sed -i \"s/vp_prior_max(1,:) =.*/"\
                "vp_prior_max(1,:) = %s/g\" " %range_log[1] +pfile

        #print sed_cmd
        run_cmd(sed_cmd)

for nf in first_nf_list:
    for prior in prior_list:
        for approx in approx_list:

	    newdir = "newrun_p%s_Nf%s_%s_%s" %(2.0,nf,prior,approx)
	    sedcmd = "sed -i \"s/multifield.*/multifield\/NEWRUN_CONSREL_SMALL\/%s/g\" "%newdir + "slurm_job.sh"
            run_cmd(sedcmd)

            modify_pfile(nf, prior, approx, p=2.0)

            #Submit the jobs
            submitcmd = ". ./submit.sh " + str(njobs)

            #print submitcmd
            run_cmd(submitcmd)

#for nf in second_nf_list:
#    for prior in prior_list:
#        for approx in approx_list:
#            for p in p_list:
#
#                newdir = "newrun_p%s_Nf%s_%s_%s" %(p,nf,prior,approx)
#                sedcmd = "sed -i \"s/multifield.*/multifield\/NEWRUN_CONSREL_BIG\/%s/g\" "%newdir + "slurm_job.sh"
#                run_cmd(sedcmd)
#
#                modify_pfile(nf, prior, approx, p)
#
#                #Submit the jobs
#                submitcmd = ". ./submit.sh " + str(njobs)
#
#                #print submitcmd
#                run_cmd(submitcmd)
