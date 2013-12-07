import numpy as np
import matplotlib.pyplot as plt
import pyplotsetup

figprops = dict(figsize=(1.0*pyplotsetup.fig_width,
    1.0*pyplotsetup.fig_height))
fig = plt.figure(**figprops)
adjustprops = dict(left=0.16, bottom=0.18, right=0.97, top=0.93, 
    wspace=0.1, hspace=0.1)
fig.subplots_adjust(**adjustprops)
ax = fig.add_subplot(111)

t = np.arange(0,2,0.01)

ax.plot(t,np.sin(t))
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$\sin(t)$')

plt.savefig('test.eps')
plt.savefig('test.pdf')

