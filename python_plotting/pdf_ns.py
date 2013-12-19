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

adjustprops = dict(left=0.14, bottom=0.15, right=0.95, top=0.98,
    wspace=0.1, hspace=0.30)

fig = plt.figure(**figprops)
fig.subplots_adjust(**adjustprops)

ax1 = fig.add_subplot(212)


#PDF data
pdf_100_isoE   =np.loadtxt('100field_pdf_ns_isoE.txt')
pdf_100_isoN   =np.loadtxt('100field_pdf_ns_isoN_300.txt')
pdf_100_isoN_60=np.loadtxt('100field_pdf_ns_isoN_60.txt')

pdf_3_isoE   =np.loadtxt('3field_pdf_ns_isoE.txt')
pdf_3_isoN   =np.loadtxt('3field_pdf_ns_isoN_300.txt')
pdf_3_isoN_60=np.loadtxt('3field_pdf_ns_isoN_60.txt')

#plt.xlim([0.93,0.96001])
plt.xlim([0.94,0.96001])

p1, = plt.plot(pdf_100_isoE[:,0], pdf_100_isoE[:,1],color='r',linestyle='-')
p2, = plt.plot(pdf_100_isoN[:,0], pdf_100_isoN[:,1],color='b',linestyle='--')
p3, = plt.plot(pdf_100_isoN_60[:,0], pdf_100_isoN_60[:,1],color='k',linestyle=':')

#plt.rc('legend',**{'fontsize':8})
#plt.legend([p1,p2,p3],
#        [r'iso-$E$; $E_0=0.1\, \mathrm{M_{Pl}}$',
#            r'iso-$N_e$; $N_e = 300$',
#            r'iso-$N_e$; $N_e = 60$'],
#        loc = 'upper left')


ax1.set_ylabel(r'$P$')
ax1.set_xlabel(r'$n_s$')

#plt.ylabel(r'$P$')
#plt.xlabel(r'$n_s$')


ax2 = fig.add_subplot(211)

#plt.ylim([0,1.1])
plt.xlim([0.94,0.970])

p1, = plt.plot(pdf_3_isoE[:,0], pdf_3_isoE[:,1],color='r',linestyle='-')
p2, = plt.plot(pdf_3_isoN[:,0], pdf_3_isoN[:,1],color='b',linestyle='--')
p3, = plt.plot(pdf_3_isoN_60[:,0],   pdf_3_isoN_60[:,1],color='k',linestyle=':')


plt.rc('legend',**{'fontsize':8})
plt.legend([p1,p2,p3],
        [r'iso-$E$; $E_0=0.1\, \mathrm{M_{Pl}}$',
            r'iso-$N_e$; $N_e = 300$',
            #r'slow-roll'],
            r'iso-$N_e$; $N_e=60$'],
        loc = 'upper left')


ax2.set_ylabel(r'$P$')
#ax2.set_xlabel(r'$n_s$')


save=True
if save:
    #name = 'pdf_100field'
    name = 'pdf_100_3field'
    #direct='/home/lpri691/LaTex/multifield_modecode/ics_and_preds/plots/'
    direct='/home/lpri691/LaTex/multifield_modecode/ics_and_preds/prl/prl_submit/'
    #direct='./'
    plt.savefig(direct+name+'.png', dpi=250)
    plt.savefig(direct+name+'.pdf')


plt.show()
