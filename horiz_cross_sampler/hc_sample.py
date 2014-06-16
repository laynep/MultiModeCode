#!/usr/bin/python

"""
Samples the horizon crossing surface for N-quadratic inflation and uses the horizon crossing approximation (HCA) to calculate observables in the \delta N formalism.  Only works for N-quadratic at the moment.

Can put arbitrary prior on the initial conditions or masses.  Builds PDFs for the observables using a histogram estimator and can calculate the derivatives with respect to the "prior parameters" of the expected value of the log-likelihood.
"""

import numpy as np
import sys


