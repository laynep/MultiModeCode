#!/usr/bin/python

"""
Module that contains the potential function and its derivatives.
"""

import numpy as np
import sys

def check_choice(choice):
    print "Choice=", choice, "not implemented"
    sys.exit()

def potential(phi, choice="Nquad", **params):

    if choice == "Nquad":
        return 0.5*np.sum(params["m2"]*phi**2)
    else:
        check_choice(choice)


def dVdphi(phi,choice="Nquad",**params):

    if choice == "Nquad":
        return params["m2"]*phi
    else:
        check_choice(choice)

def d2Vdphi2(phi,choice="Nquad",**params):

    if choice == "Nquad":
        return params["m2"]*np.identity(phi.size)
    else:
        check_choice(choice)

def d3Vdphi3(phi,choice="Nquad",**params):

    if choice == "Nquad":
        return np.zeros((phi.size,phi.size,phi.size))
    else:
        check_choice(choice)


def V_i(phi, choice="Nquad",**params):
    """For sum-separable potentials, returns V_i(phi_i)."""

    if choice == "Nquad":
        return 0.5*params["m2"]*phi**2
    else:
        check_choice(choice)


def H_horiz_cross(scale, choice="Nquad",**params):
    """Hubble parameter as scale leaves horizon."""
    print "test H_horiz_cross"
    return
