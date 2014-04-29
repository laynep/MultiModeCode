import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats

import pyplotsetup

import sys

##################
#Figure options
figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 1.0*pyplotsetup.fig_height))

adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.96,
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


##################
#Make plot
gs = gridspec.GridSpec(2,1,height_ratios=[4,1])

#ax1 = fig.add_subplot(111)
ax1 = plt.subplot(gs[0])

plots=dict([])
leg={1:r'$|\delta P|^2$',
     2:r'$|\delta P_\mathrm{nad}|^2$',
     3:r'$|\delta P_\mathrm{ad}|^2$',
     4:r'$|\delta P_\mathrm{nad}\delta P_\mathrm{ad}^*|$'}
     #5:'$|\deltaP_\mathrm{ad}P_\mathrm{nad}^*|$'}

color={1:'k',
       2:'b',
       3:'r',
       4:'c'}

line= {1:1.5,
       2:2.0,
       3:1.5,
       4:1.5}

linestyle= {1:'-',
            2:'-',
            3:'-',
            4:'-'}

#plt.plot(np.log(data[:,0]/a_pivot),abs(data[:,1]*10**-12.0),
#        'k--', linewidth=1.,alpha=0.25)
plt.fill_between(np.log(data[:,0]/a_pivot),abs(data[:,1]*10**-15.0),1e-60,
        color='k',
        alpha=0.1)


for i in leg:
    plots[i],=plt.plot(np.log(data[:,0]/a_pivot),abs(data[:,i]),
            color=color[i],
            linewidth=line[i],
            linestyle=linestyle[i])


ax1.legend([v for k, v in plots.iteritems()],
        [values for keys, values in leg.iteritems()],
        loc='lower left')

ax1.set_yscale('log')
#ax1.set_xscale('log')


plt.xlim(0.0,55.0)
plt.ylim(10.0**-48,10**-27)

ax1.set_ylabel(r'$\mathcal P$')

plt.setp(ax1.get_xticklabels(), visible=False)


#Residuals
ax2 = plt.subplot(gs[1])

plt.fill_between(np.log(data[:,0]/a_pivot),abs(data[:,1]*10**-15.0),1e-70,
        color='k',
        alpha=0.1)

plt.plot(np.log(data[:,0]/a_pivot),abs(data[:,1]-data[:,3]-data[:,4]-data[:,5]),
        'b--',
        alpha=0.8)
ax2.set_yscale('log')

ax2.set_xlabel(r'$N_e$')
ax2.set_ylabel(r'$\mathcal P$')

plt.xlim(0.0,55.0)
plt.ylim(10**-50,10**-28)
#ax2.set_yticks([ 10**-36, 10**-48, 10**-60])
ax2.set_yticks([ 10**-36, 10**-48])

save=True
name = 'isocurv_compare'
if save:
    direct='/home/lpri691/LaTex/multifield_modecode/multimodecode/plots/'
    #plt.savefig(direct+name+'.png', dpi=250)
    plt.savefig(direct+name+'.pdf')


plt.show()
