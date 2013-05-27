
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('phiarr.txt', unpack=True)

size = np.shape(data)
nfields = size[0]
nsteps = size[1]

x = data[0,:]
colors = ['b', 'r', 'b', 'g', 'w']
print x

plt.figure()
for icol in range(1, nfields):
    y = data[icol,:]
    print y
    plt.plot(x, y, color = colors[icol])
    #plt.xlim((0, 10))
    #plt.ylim((-2, 2))
    plt.xlabel('N')
    plt.ylabel('Field')
    plt.title('Field evolution')
plt.show()

plt.figure()
for icol in range(1, nfields-1):
    plt.subplot(2, 2, icol)
    y1 = data[icol,:]
    y2 = data[icol+1,:]
    plt.plot(y1, y2)
    plt.xlabel('Field %s'%icol)
    plt.ylabel('Field %s'%(icol+1))
plt.show()

