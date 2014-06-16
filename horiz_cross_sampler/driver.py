#!/usr/bin/python

"""
Simple driver function.
"""

import numpy as np
from classes import *

#TEST
print "testing classes.py"
N_piv=55.0
nfields = 2
phi = np.sqrt(4.0*N_piv/nfields)
phi = np.array([phi for i in xrange(nfields)])
phi_end = np.array([1.0e0, 2.0e-1])
m2 = 10.0**np.array([-10.0, -8.0])


x=deltaN_model(HC_approx=False, model="Nquad",nfields=1)
x.load_params(m2=m2)

x.calc_observs(phi, phi_end)
print "observables:", x.observ
