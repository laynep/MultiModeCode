import numpy as np
import sys
import scipy.stats as stats


import matplotlib.pyplot as plt

#Cosmology data
dataEE  =np.loadtxt('100field_nsralpha_isoE.txt')
dataisoN=np.loadtxt('100field_nsralpha_isoN_300.txt')
dataisoN60  =np.loadtxt('100field_nsralpha_isoN_60.txt')

nsEE=dataEE[:,0]
nsisoN=dataisoN[:,0]
nsisoN60=dataisoN60[:,0]

print '100 field KS EE-isoN:', stats.ks_2samp(nsEE,nsisoN)
print '100 field KS EE-isoN60:', stats.ks_2samp(nsEE,nsisoN60)
print '100 field KS isoN-isoN60:', stats.ks_2samp(nsisoN60,nsisoN)

dataEE  =np.loadtxt('3field_nsralpha.txt')
dataisoN=np.loadtxt('3field_nsralpha_isoN.txt')
dataSR  =np.loadtxt('3field_nsralpha_SR.txt')

nsEE=dataEE[:,0]
nsisoN=dataisoN[:,0]
nsSR=dataSR[:,0]

print '3 field KS EE-isoN:', stats.ks_2samp(nsEE,nsisoN)
print '3 field KS EE-SR:', stats.ks_2samp(nsEE,nsSR)
print '3 field KS isoN-SR:', stats.ks_2samp(nsSR,nsisoN)


