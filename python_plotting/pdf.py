import matplotlib.pyplot as plt
import numpy as np
import pyplotsetup


#!!!!!!!!!!!
#Figure options

figprops = dict(figsize=(1.0*pyplotsetup.fig_width, 0.8*pyplotsetup.fig_height))

adjustprops = dict(left=0.15, bottom=0.2, right=0.97, top=1.0,
    wspace=0.1, hspace=0.02)

fig = plt.figure(**figprops)
fig.subplots_adjust(**adjustprops)

ax1 = fig.add_subplot(111)


#CDF data
pdf_100_isoE=np.loadtxt('risoN_300_PDF.txt')
pdf_100_isoN=np.loadtxt('risoE_PDF.txt')
pdf_100_isoN_60=np.loadtxt('risoN_60_PDF.txt')

plt.ylim([0,1200.])


p1, = plt.plot(pdf_100_isoN[:,0], pdf_100_isoN[:,1]      ,color='r', linestyle='-')
p2, = plt.plot(pdf_100_isoE[:,0], pdf_100_isoE[:,1]      ,color='b', linestyle='--')
p3, = plt.plot(pdf_100_isoN_60[:,0], pdf_100_isoN_60[:,1],color='k', linestyle=':')


plt.rc('legend',**{'fontsize':8})
plt.legend([p1,p2,p3],
        [r'iso-$E$; $E_0=0.1\, \mathrm{M_{Pl}}$',
            r'iso-$N_e$; $N_e = 300$',
            r'iso-$N_e$; $N_e = 60$'])


plt.ylabel(r'$P$')
plt.xlabel(r'$r_\mathrm{iso}$')

plt.xticks(np.arange(0.00, 0.01, 0.004))
plt.yticks(np.arange(0.00, 1200, 550))

#ax2 = fig.add_subplot(212)
#plt.plot([0.94,0.95,0.96])
#plt.ylabel('test')


save=True
if save:
    name = 'ratio_iso_100field'
    #direct='/home/lpri691/LaTex/multifield_modecode/ics_and_preds/plots/'
    direct='/home/lpri691/LaTex/multifield_modecode/ics_and_preds/prl/prl_submit/'
    #direct='./'
    plt.savefig(direct+name+'.png', dpi=250)
    plt.savefig(direct+name+'.pdf')


plt.show()
