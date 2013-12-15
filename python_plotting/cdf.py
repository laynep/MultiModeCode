import matplotlib.pyplot as plt
import numpy as np
import pyplotsetup
import sys


#!!!!!!!!!!!
#Figure options

#figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 2.0*pyplotsetup.fig_height))
figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 1.0*pyplotsetup.fig_height))

#print 'Aspect Ratio: ', pyplotsetup.fig_width / (1.0*pyplotsetup.fig_height)
#sys.exit()

adjustprops = dict(left=0.13, bottom=0.15, right=0.95, top=1.0,
    wspace=0.1, hspace=0.30)

fig = plt.figure(**figprops)
fig.subplots_adjust(**adjustprops)

ax1 = fig.add_subplot(211)


#CDF data
cdf_100_isoE   =np.loadtxt('100field_cdf_ns_isoE.txt')
cdf_100_isoN   =np.loadtxt('100field_cdf_ns_isoN_300.txt')
cdf_100_isoN_60=np.loadtxt('100field_cdf_ns_isoN_60.txt')

cdf_3_isoE   =np.loadtxt('3field_cdf_ns_isoE.txt')
cdf_3_isoN   =np.loadtxt('3field_cdf_ns_isoN.txt')
cdf_3_SR=np.loadtxt('3field_cdf_ns_SR.txt')

plt.ylim([0,1.1])

p1, = plt.plot(cdf_100_isoE[:,0], cdf_100_isoE[:,1],'r')
p2, = plt.plot(cdf_100_isoN[:,0], cdf_100_isoN[:,1],'b')
p3, = plt.plot(cdf_100_isoN_60[:,0], cdf_100_isoN_60[:,1],'k')


plt.rc('legend',**{'fontsize':8})
plt.legend([p1,p2,p3],
        [r'iso-$E$; $E_0=0.1\, \mathrm{M_{Pl}}$',
            r'iso-$N_e$; $N_e = 300$',
            r'iso-$N_e$; $N_e = 60$'],
        loc = 'upper left')


ax1.set_ylabel(r'$P$')

#plt.ylabel(r'$P$')
#plt.xlabel(r'$n_s$')


ax2 = fig.add_subplot(212)

plt.ylim([0,1.1])

p1, = plt.plot(cdf_3_isoE[:,0], cdf_3_isoE[:,1],'r')
p2, = plt.plot(cdf_3_isoN[:,0], cdf_3_isoN[:,1],'b')
p3, = plt.plot(cdf_3_SR[:,0],   cdf_3_SR[:,1],'k')


plt.rc('legend',**{'fontsize':8})
plt.legend([p1,p2,p3],
        [r'iso-$E$; $E_0=0.1\, \mathrm{M_{Pl}}$',
            r'iso-$N_e$; $N_e = 100$',
            r'slow-roll'],
        loc = 'upper left')


ax2.set_ylabel(r'$P$')
ax2.set_xlabel(r'$n_s$')


save=False
if save:
    #name = 'cdf_100field'
    name = 'cdf_100_3field'
    direct='/home/lpri691/LaTex/multifield_modecode/ics_and_preds/plots/'
    #direct='./'
    plt.savefig(direct+name+'.png', dpi=250)
    plt.savefig(direct+name+'.pdf')


plt.show()
