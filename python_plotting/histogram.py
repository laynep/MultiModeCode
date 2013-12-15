import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from scipy import stats

import pyplotsetup

import sys


iso_E = True

if iso_E:

    #Cosmology data
    #data=np.loadtxt('nsralpha.txt')
    data=np.loadtxt('100field_nsralpha_isoE.txt')

    ns=data[:,0]
    r=data[:,1]
    alpha=1000*data[:,2]

    ns_min, ns_max = min(ns),max(ns)
    #ns_min, ns_max = 0.945, 0.955
    #r_min, r_max =min(r)-0.001,max(r)+0.001

    r_min, r_max =min(r),max(r)

    rzoom_min, rzoom_max = 0.1440, 0.146
    alpha_min, alpha_max =min(alpha),max(alpha)

    r_tick_min, r_tick_max, r_units = -0.8405, r_max, 0.0005

else:

    #Cosmology data
    data=np.loadtxt('nsralpha_isoN_100fields.txt')

    ns=data[:,0]
    r=data[:,1]
    #r=data[:,1]
    alpha=1000*data[:,2]

    #ns_min, ns_max = min(ns),max(ns)
    ns_min, ns_max = 0.962,1.0001*max(ns)
    #ns_min, ns_max = 0.9625, 0.9645

    #r_min, r_max =min(r),max(r)
    r_min, r_max =np.log10([0.144030, 0.1441])


    rzoom_min, rzoom_max = 0.1440, 0.146
    alpha_min, alpha_max =min(alpha),max(alpha)

    #r_tick_min, r_tick_max, r_units = -0.8405, r_max, 0.0005


#Histogram
doing_histogram = True
if doing_histogram:
    bins=120
    binsnsr=[90,84] #Scott method from Mathematica
    binsnsalpha=[90,164] #Scott method from Mathematica

    norm=False

    nsr, xedges1, yedges1 = np.histogram2d(ns, r,
            range=[[ns_min, ns_max],[r_min, r_max]],
            bins=binsnsr, normed=norm)

    nsalpha, xedges2, yedges2 = np.histogram2d(ns, alpha,
            range=[[ns_min, ns_max],[alpha_min, alpha_max]],
            bins=binsnsalpha, normed=norm)

    nsr_zoom, xedges_zoom, yedges_zoom = np.histogram2d(ns, r,
            range=[[ns_min, ns_max],[rzoom_min, rzoom_max]],
            bins=bins, normed=norm)

    extent_nsr = [xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]]

    extent_nsr_zoom = [xedges_zoom[0], xedges_zoom[-1],
            yedges_zoom[0], yedges_zoom[-1]]
    extent_nsalpha = [xedges1[0], xedges1[-1], yedges2[0], yedges2[-1]]


    #Normalize to "probability"
    nsr = nsr/float(len(ns))
    nsalpha = nsalpha/float(len(ns))


# Perform a kernel density estimate (KDE) on the data
doing_kde=False
if doing_kde:
    values_nsr = np.vstack([ns, r])
    values_nsalpha = np.vstack([ns, alpha])

    kernel_nsr = stats.gaussian_kde(values_nsr)
    kernel_nsalpha = stats.gaussian_kde(values_nsalpha)

    x_nsr, y_nsr = np.mgrid[ns_min:ns_max:100j, r_min:r_max:100j]
    x_nsalpha, y_nsalpha = np.mgrid[ns_min:ns_max:100j, alpha_min:alpha_max:100j]

    positions_nsr = np.vstack([x_nsr.ravel(), y_nsr.ravel()])
    positions_nsalpha = np.vstack([x_nsalpha.ravel(), y_nsalpha.ravel()])

    f_nsr = np.reshape(kernel_nsr(positions_nsr).T, x_nsr.shape)
    f_nsalpha = np.reshape(kernel_nsalpha(positions_nsalpha).T, x_nsalpha.shape)



#!!!!!!!!!!!
#Figure options

#figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 2.0*pyplotsetup.fig_height))
figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 1.85*pyplotsetup.fig_height))

adjustprops = dict(left=0.19, bottom=0.09, right=0.97, top=0.96,
    wspace=0.1, hspace=0.02)

imshowasp = 'auto'

#Got these from http://matplotlib.org/examples/color/colormaps_reference.html
color_map = plt.get_cmap('jet')  #DEFAULT

#color_map = plt.get_cmap('PuBu')
#color_map = plt.get_cmap('seismic')
#color_map = plt.get_cmap('hot_r')
#color_map = plt.get_cmap('gist_heat')
#color_map = plt.get_cmap('binary')

#Interpolation for imshow
interp = 'nearest'


#-----------

fig = plt.figure(**figprops)
fig.subplots_adjust(**adjustprops)

ax1 = fig.add_subplot(211)

ax1.set_ylabel(r'$r$')

if doing_kde:
    plt.imshow(f_nsr.T, origin='lower', extent=extent_nsr, aspect=imshowasp,
            cmap=color_map, interpolation=interp)
else:
    #plt.yticks(np.arange(r_tick_min, r_tick_max, r_units))

    plt.imshow(nsr.T, origin='lower', extent=extent_nsr, aspect=imshowasp,
            cmap=color_map, interpolation=interp)


#-----------

ax2=plt.subplot(212)
ax2.set_ylabel(r'$10^3 \alpha$')
ax2.set_xlabel(r'$n_s$')

if doing_kde:
    plt.imshow(f_nsalpha.T, origin='lower', extent=extent_nsalpha, aspect=imshowasp,
            cmap=color_map, interpolation=interp)
else:
    plt.imshow(nsalpha.T, origin='lower', extent=extent_nsalpha, aspect=imshowasp,
            cmap=color_map, interpolation=interp)
    #inset_plot = fig.add_axes([0.275,0.675,0.4,0.25])
    #plt.imshow(nsr_zoom.T, origin='lower', extent=extent_nsr_zoom,
    #        aspect=imshowasp, cmap=color_map, interpolation=interp)

#Remove x-axis label top fig (ax1)
plt.setp(ax1.get_xticklabels(), visible=False)

save=True
if iso_E:
    if doing_kde:
        name = 'nsralpha_100field_kde'
    else:
        name = 'nsralpha_100field'
else:
    if doing_kde:
        name = 'nsralpha_100field_kde_isoN'
    else:
        name = 'nsralpha_100field_isoN'
if save:
    direct='/home/lpri691/LaTex/multifield_modecode/ics_and_preds/plots/'
    #direct='./'
    plt.savefig(direct+name+'.png', dpi=250)
    plt.savefig(direct+name+'.pdf')


plt.show()
