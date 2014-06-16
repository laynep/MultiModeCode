#!/usr/bin/python

"""
Module that defines all the classes, etc for inflationary calculations.
"""

import numpy as np
import potential as pot
import hc_sample as sample

import sys


class inflation_model:
    """
Class for inflationary models, contains a potential with derivatives, methods to set parameters (priors), and methods to calculate observables.
    """

    model = None
    nfields = None

    params = {
            "Nquad":{ "m2":None }
            }

    def __init__(self, model, nfields):
        self.model = model
        self.nfields = nfields

    def load_params(self, **params):
        """Load the model parameters."""
        self.params[self.model] = params

    #Wrappers to potential.py
    def V(self,phi):
        """Potential energy density."""
        return pot.potential(phi,self.model, **self.params[self.model])

    def dV(self,phi):
        """Derivative of potential wrt fields."""
        return pot.dVdphi(phi,self.model, **self.params[self.model])

    def d2V(self,phi):
        """Second derivative of potential wrt fields."""
        return pot.d2Vdphi2(phi,self.model, **self.params[self.model])

    def d3V(self,phi):
        """Third derivative of potential wrt fields."""
        return pot.d3Vdphi3(phi,self.model, **self.params[self.model])

    def V_i(self, phi):
        """For sum-separable potentials, returns V_i(phi_i)."""
        return pot.V_i(phi,choice=self.model, **self.params[self.model])

    #Slow-roll parameters.
    def eps_i(self, phi):
        """For sum-separable potentials, returns eps_i = (1/2) (dV_i/V)**2."""

        dV = self.dV(phi)
        V = self.V(phi)
        return 0.5e0*(dV**2)/V**2

    def eta_i(self, phi):
        """For sum-separable potentials, returns eta_i = d2V_ii**2/V"""
        d2V = self.d2V(phi)
        V = self.V(phi)
        return np.diagonal(d2V)**2/V





class deltaN_model(inflation_model):
    """
    The \delta N formulation for observables from Vernizzi-Wands (astro-ph/0603799) and Battefeld-Easther (astro-ph/0610296).  Requires sum-separable potential and assumes massless, Gaussian random fields at horizon crossing.  If using the horizon crossing approximation (HCA), then ignores the final surface.
    """

    HCA = False
    phi_zero = None

    def __init__(self, HC_approx, **infl_args):
        inflation_model.__init__(self,**infl_args) #init the parent class
        self.HCA = HC_approx
        self.phi_zero = np.zeros(self.nfields)


    def Z_i(self, phi_end=phi_zero):
        """The function Z from Battefeld-Easther (astro-ph/0610296) that encodes the dependence of the observables on the end of inflation surface.  This contribution is ignored in the horizon crossing approximation (HC_approx)."""
        if self.HCA:
            return np.zeros(self.nfields)
        else:
            eps_i = self.eps_i(phi_end)
            eps = np.sum(eps_i)
            V = self.V(phi_end)
            return V*eps_i/eps - self.V_i(phi_end)

    def dZ_ij(self, phi_hc, phi_end=phi_zero):
        """Derivative of Z_i with respect to fields at horizon crossing."""
        if self.HCA:
            return np.zeros((self.nfields,self.nfields))
        else:

            eps_end = self.eps_i(phi_end)
            eps_hc = self.eps_i(phi_hc)
            eps_t_end = sum(eps_end)
            V_end = self.V(phi_end)
            V_hc = self.V(phi_hc)
            eta_end = self.eta_i(phi_end)

            delta = np.identity(phi_hc.size)

            mat1 = np.ones((phi_hc.size,phi_hc.size))*(eps_end/eps_t_end) - delta
            mat1=mat1.T
            mat2 = eps_end*(1.0-eta_end/eps_t_end)

            nfields = xrange(self.nfields)

            dZ_jkl = np.array([[[
                np.sqrt(2.0/eps_end[k])*
                mat1[l,j]*mat1[k,j]*mat2[j]
                for k in nfields]
                for l in nfields]
                for j in nfields])

            dZ_jkl = (-V_end**2/V_hc)*np.sum(dZ_jkl,0)

            return dZ_jkl

    def dNdphi(self, phi_hc, phi_end=phi_zero):
        """Derivative of N_total with respect to the field position at horizon crossing."""
        eps_i = self.eps_i(phi_hc)
        V_i = self.V_i(phi_hc)
        Z_i = self.Z_i(phi_end)
        V = self.V(phi_hc)
        return ((1.0/np.sqrt(2.0*eps_i))/V)*(V_i + Z_i)

    def d2Ndphi2(self, phi_hc, phi_end=phi_zero):
        """Second derivative of N_total with respect to the field position at horizon crossing."""

        eta_hc = self.eta_i(phi_hc)
        eps_hc = self.eps_i(phi_hc)
        V_hc = self.V(phi_hc)
        V_i_hc = self.V_i(phi_hc)
        Z_i = self.Z_i(phi_end)
        dZ_ij = self.dZ_ij(phi_hc, phi_end)

        delta = np.identity(self.nfields)

        nfields = xrange(self.nfields)
        d2N = np.array(
                [[delta[k,l]*
                    (1.0 -(eta_hc[l]/2.0/eps_hc[l])*
                        (V_i_hc[l]+Z_i[l])/V_hc) +
                (1.0/np.sqrt(2.0*eps_hc[l])/V_hc)*dZ_ij[l,k]
                for k in nfields]
                for l in nfields])

        return d2N


    def PR(self, phi_hc, phi_end=phi_zero):
        """Calculates the amplitude of the scalar power spectrum as a function of the horizon crossing and ending phi values for a given choice of potential."""

        V_hc = self.V(phi_hc)
        H_hc = np.sqrt(V_hc/3.0)
        dN = self.dNdphi(phi_hc,phi_end)

        massless_P = (H_hc/2.0/np.pi)**2

        return np.sum(dN*dN)*massless_P

    def r(self, phi_hc, phi_end=phi_zero):
        """Tensor-to-scalar ratio."""

        dN = self.dNdphi(phi_hc,phi_end)
        return 8.0/np.sum(dN*dN)

    def n_t(self, phi_hc):
        """Tensor spectral tilt."""

        eps_i = self.eps_i(phi_hc)
        return -2.0*np.sum(eps_i)/(1.0 - np.sum(eps_i))

    def f_NL(self, phi_hc, phi_end=phi_zero):
        """Local non-Gaussianity parameter f_NL for the bispectrum with phi ~ f_NL[ phi^2 - <phi^2>]."""

        dN = self.dNdphi(phi_hc,phi_end)
        d2N = self.d2Ndphi2(phi_hc,phi_end)

        fnl = np.einsum('i,j,ij',dN,dN,d2N)
        fnl *= (-5.0/6.0)/(sum(dN*dN))**2

        return fnl

    def tau_NL(self, phi_hc, phi_end=phi_zero):
        """Local non-Gaussianity parameter tau_NL for the trispectrum."""

        dN = self.dNdphi(phi_hc,phi_end)
        d2N = self.d2Ndphi2(phi_hc,phi_end)
        taunl = np.einsum('ab,ac,b,c', d2N, d2N, dN, dN)
        taunl /= np.sum(dN*dN)**3

        return taunl

    def n_s(self, phi_hc, phi_end=phi_zero):
        """Scalar spectral index."""

        dV = self.dV(phi_hc)
        d2V = self.d2V(phi_hc)
        eps = np.sum(self.eps_i(phi_hc))
        dN = self.dNdphi(phi_hc,phi_end)
        V = self.V(phi_hc)

        ns = np.einsum('ij,i,j', d2V, dV, dV)
        ns *= (1.0/V/np.sum(dN*dN))

        ns += 1.0 - 2.0*eps - (2.0/sum(dN*dN))

        return ns

    def alpha_s(self, phi_hc, phi_end=phi_zero):
        """Running of the scalar spectral index."""

        V = self.V(phi_hc)
        dV = self.dV(phi_hc)
        d2V = self.d2V(phi_hc)
        d3V = self.d3V(phi_hc)
        dN = self.dNdphi(phi_hc,phi_end)

        term = []
        term.append(np.einsum('a,b,ab', dV, dV, d2V)*(-2.0/V**3))

        term.append(np.einsum('a,a', dV, dV)**2*(2.0/V**4))

        term.append((4.0/V/np.sum(dN*dN)**2)*(V - np.einsum('a,b,ab',dN,dN,d2V))**2)

        term.append((2.0/V/np.sum(dN*dN))*np.einsum('a,b,c,abc',dN,dN,dV,d3V))

        term.append((4.0/V/np.sum(dN*dN))* (np.einsum('c,b,bc',dV,dN,d2V) -
                np.einsum('a,ac,b,bc',dN,d2V,dN,d2V)))

        return sum(term)


    def calc_observs(self,phi_hc,phi_end=phi_zero):
        self.observ={}

        self.observ["PR"]=self.PR(phi_hc, phi_end)
        self.observ["r"]=self.r(phi_hc,phi_end)
        self.observ["n_s"]=self.n_s(phi_hc,phi_end)
        self.observ["alpha_s"]=self.alpha_s(phi_hc,phi_end)
        self.observ["n_t"]=self.n_t(phi_hc)
        self.observ["f_NL"]=self.f_NL(phi_hc,phi_end)
        self.observ["tau_NL"]=self.tau_NL(phi_hc,phi_end)


class universe(inflation_model):
    """
Realization of a universe given an inflationary model.
Extends inflation_model to include a set of parameters and resulting observables.
    """

