#!/usr/bin/python

"""
Module that defines all the classes, etc for inflationary calculations.
"""

import numpy as np
import potential as pot

import sys


class inflation_model:
    """
Class for inflationary models, contains a potential with derivatives,  and methods to calculate observables.
    """

    def __init__(self, model, nfields):
        self.model = model
        self.nfields = nfields

    def load_params(self, **params):
        """Load the model parameters."""
        self.params = {}
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

    def xi_i(self, phi):
        """For sum-separable potentials, returns xi_i = dV_i*d3V_iii/V**2"""
        d3V = self.d3V(phi)
        dV = self.dV(phi)
        V = self.V(phi)

        d3V = [d3V[i,i,i] for i in xrange(self.nfields)]

        return dV*d3V/V**2




class deltaN_model(inflation_model):
    """
    The \delta N formulation for observables from Vernizzi-Wands (astro-ph/0603799) and Battefeld-Easther (astro-ph/0610296).  Requires sum-separable potential and assumes massless, Gaussian random fields at horizon crossing.  If using the horizon crossing approximation (HCA), then ignores the final surface.
    """

    phi_zero = None

    def __init__(self, HC_approx, **infl_args):
        inflation_model.__init__(self,**infl_args) #init the parent class
        self.HCA = HC_approx
        #self.phi_zero = np.zeros(self.nfields)
        #self.phi_zero = np.ones(self.nfields)*1e-2


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

            #Terrible...
            dZ_jkl = np.array([[[
                np.sqrt(2.0/eps_end[k])* \
                mat1[l,j]*mat1[k,j]*mat2[j] \
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

        coeff = 0.5/V_hc
        coeff2 = 1.0/np.sqrt(2.0)/V_hc
        d2N = np.array( delta - coeff*np.diag(eta_hc*V_i_hc/eps_hc) \
                - coeff*np.diag(eta_hc*Z_i/eps_hc) \
                + coeff2*np.dot( np.diag(1.0/np.sqrt(eps_hc)),dZ_ij))

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

    def n_t(self, phi_hc, phi_end=phi_zero):
        """Tensor spectral tilt."""

        eps = np.sum(self.eps_i(phi_hc))
        return -2.0*eps/(1.0 - eps)
        #return -2.0*eps

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

        ns = np.einsum('ij,i,j', d2V, dN, dN)
        ns *= (2.0/V/np.sum(dN*dN))

        ns += 1.0 - 2.0*eps - (2.0/sum(dN*dN))

        return ns

    def alpha_s(self, phi_hc, phi_end=phi_zero):
        """Running of the scalar spectral index."""

        eps_i = self.eps_i(phi_hc)
        eps = np.sum(eps_i)
        eta_i = self.eta_i(phi_hc)
        xi_i = self.xi_i(phi_hc)
        u_i = (self.V_i(phi_hc) + self.Z_i(phi_end))/self.V(phi_hc)

        sum1 = np.sum(u_i**2/eps_i)

        print eps
        print eps_i
        print eta_i
        print xi_i
        print u_i
        print sum1

        alpha_s = -8.0*eps**2
        alpha_s += 4.0*np.sum(eps_i*eta_i)
        alpha_s += (-16.0/sum1**2)*(1.0 - np.sum(eta_i*u_i**2/2.0/eps_i))**2
        alpha_s += (-8.0/sum1)*np.sum(eta_i*u_i*(1.0-eta_i*u_i/2.0/eps_i))
        alpha_s += (4.0*eps/sum1)*np.sum(eta_i*u_i**2/eps_i)
        alpha_s += (-2.0/sum1)*np.sum(xi_i*u_i**2/eps_i)

        return alpha_s


    def calc_observs(self, phi_hc,phi_end=phi_zero, obs_to_calc=None):
        """Calculate the slow-roll observables for a given value of the fields at horizon crossing and the end of inflation (ignored if using horizon crossing approximation)."""

        obs_dict={"PR":self.PR,
                "r":self.r,
                "n_s":self.n_s,
                "alpha_s":self.alpha_s,
                "n_t":self.n_t,
                "f_NL":self.f_NL,
                "tau_NL":self.tau_NL}

        if obs_to_calc == None:
            obs_to_calc = obs_dict.keys()

        self.observ={}

        for key in obs_to_calc:
            self.observ[key]=obs_dict[key](phi_hc,phi_end)


class universe(deltaN_model):
    """
Realization of a universe given an inflationary model.  Extends inflation_model to include methods to set parameters (priors) and resulting observables.

Samples the horizon crossing surface for N-quadratic inflation and uses the horizon crossing approximation (HCA) to calculate observables in the \delta N formalism.  Only works for N-quadratic at the moment.

Can put arbitrary prior on the initial conditions or masses.  Builds PDFs for the observables using a histogram estimator and can calculate the derivatives with respect to the "prior parameters" of the expected value of the log-likelihood.
    """


    def __init__(self, sampler=None, N_pivot=55.0, HC_approx=True, **infl_args):
        deltaN_model.__init__(self,HC_approx,**infl_args)

        self.N_pivot = N_pivot

        self.load_sampler(sampler)

    #Sampling routines
    def load_sampler(self, sampler):
        """Load the sampling routines."""
        if sampler==None:
            self.sampler=None
            return

        if self.model != "Nquad":
            raise TypeError("Model not implemented in sampler.")

        #Choices for sampling techniques
        sampling_techn = {
                "constant":self.constant,
                "MP_and_uniformsphere":self.MP_and_uniformsphere
                }

        self.sampler = sampling_techn[sampler]

    def get_new_params(self, **samp_params):
        """Get a new set of parameters and horizon crossing field values."""

        if self.model != "Nquad":
            raise TypeError("Trying to set horizon crossing field value \
                    in model that isn't N-quadratic.  Not implemented.")

        if self.sampler != None:
            params, self.phi_hc = self.sampler(**samp_params)

        self.load_params(m2=params)


    def constant(self):
        """Returns the initial conditions and masses set to a constant value."""

        phi = np.sqrt(4.0*self.N_pivot/self.nfields)
        phi = np.array([phi for i in xrange(self.nfields)])

        phi_end = 1e-2*phi

        m2 = 10.0**np.array([-10.0+2*i for i in xrange(self.nfields)])

        return m2, phi

    def MP_and_uniformsphere(self,  nmoduli, radius, m_avg=1.5e-5):
        """ Samples the Marcenko-Pasteur distribution with parameter beta=#fields/(#fields+#moduli) by building an (N+P)xN random matrix, for N=#axions (fields) and P=#moduli, with entries drawn from a Gaussian with zero mean and variance sigma^2, which is fixed by requiring COBE normalization and noting sigma^2=<m^2>.

Does a uniform sampling of a sphere with given radius for initial conditions."""

        #Masses:
        #GR Matrix
        #sigma^2 = <m^2> ~ (1.5e-5 Mpl)^2
        sigma = m_avg

        mat=np.random.normal(0.0, sigma,
                (nmoduli, self.nfields))
        mat=np.dot(mat.T,mat)
        m2, eigvect = np.linalg.eigh(mat) #m. faster than eigvals

        #ICs on sphere
        mat = np.random.normal(0.0, 1.0, self.nfields)
        ICs = (radius/np.sqrt(np.sum(mat*mat)))*mat

        return m2, ICs
