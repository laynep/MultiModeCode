import numpy as np
import matplotlib.pyplot as plt

import pyplotsetup

import sys

##################
#Figure options
figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 0.66*pyplotsetup.fig_height))

adjustprops = dict(left=0.1, bottom=0.145, right=0.97, top=0.96,
    wspace=0.1, hspace=0.02)

fig = plt.figure(**figprops)
fig.subplots_adjust(**adjustprops)
legendparams = {'legend.fontsize': 8,
          'legend.linewidth': 3}
plt.rcParams.update(legendparams)

##################
#Data
data = np.loadtxt('/home/lpri691/Cosmology_Research/F_projects/ModeCode/ModeCode_multifield/output.txt')

a_pivot =  1.4976849E-55

#Make plots
ax1 = fig.add_subplot(111)

plots=dict([])
leg={6:r'$\mathcal P_\mathcal{S}$',
     7:r'$\mathcal P_\mathrm{ent}$',
     8:r'$\mathcal P_\mathcal{R}$'}

color={6:'k',
       7:'b',
       8:'g'}

line= {6:2.0,
       7:1.5,
       8:1.5}

for i in leg:
    plots[i],=plt.plot(np.log(data[:,0]/a_pivot),abs(data[:,i]),
            color=color[i],
            linewidth=line[i])

ax1.legend([v for k, v in plots.iteritems()],
        [values for keys, values in leg.iteritems()],
        loc='lower left')

plt.xlim(0.0,55.0)
plt.ylim(10.0**-52,10**-7)


ax1.set_yscale('log')

ax1.set_xlabel(r'$N_e$')
ax1.set_ylabel(r'$\mathcal P$')

save=True
name = 'iso_vs_ent'
if save:
    direct='/home/lpri691/LaTex/multifield_modecode/multimodecode/plots/'
    plt.savefig(direct+name+'.pdf')

plt.show()
