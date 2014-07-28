#!/usr/bin/python

"""
Module that contains the potential function and its derivatives.
"""

import numpy as np
import sys

def check_choice(choice):
    TypeError("Model choice ", choice, "not implemented")

def potential(phi, choice, **params):

    if choice == "Nquad":
        return 0.5*np.sum(params["m2"]*phi**2)
    elif choice == "Nmono":
        p = params["p"]
        return (1.0/p)*np.sum(params["lambd"]*phi**p)
    else:
        check_choice(choice)

def dVdphi(phi,choice,**params):

    if choice == "Nquad":
        return params["m2"]*phi
    elif choice == "Nmono":
        p = params["p"]
        return params["lambd"]*phi**(p-1)
    else:
        check_choice(choice)

def d2Vdphi2(phi,choice,**params):

    if choice == "Nquad":
        return np.diag(params["m2"])
    elif choice == "Nmono":
        p = params["p"]
        d2V = np.zeros((phi.size,phi.size))
        d2V[np.diag_indices_from(d2V)] = (p-1.0)*params["lambd"]*phi**(p-2)
        return d2V
    else:
        check_choice(choice)

def d3Vdphi3(phi,choice,**params):

    if choice == "Nquad":
        return np.zeros((phi.size,phi.size,phi.size))
    elif choice == "Nmono":
        p = params["p"]
        d3V = np.zeros((phi.size,phi.size,phi.size))
        d3V[np.diag_indices_from(d2V)] = (p-1.0)*(p-2.0)*params["lambd"]*phi**(p-3)

        return d3V
    else:
        check_choice(choice)


def V_i(phi, choice,**params):
    """For sum-separable potentials, returns V_i(phi_i)."""

    if choice == "Nquad":
        return 0.5*params["m2"]*phi**2
    elif choice == "Nmono":
        p = params["p"]
        return (1.0/p)*params["lambd"]*phi**p
    else:
        check_choice(choice)
