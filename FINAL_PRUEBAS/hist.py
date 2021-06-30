import matplotlib.pyplot as plt
import numpy as np
import qexpy as q
import sys

def maxwell_boltzmann(v,a,b):
    return a * v * q.exp(-b*v**2)

vel = np.loadtxt("velocities_data.txt", unpack=True)
Nbins = int(sys.argv[1])

n,bins,patches = plt.hist(vel,density=True,edgecolor='black',bins=Nbins)
plt.xlabel('Velocidad')
plt.ylabel('Frecuencia')
plt.savefig('Histogramas.pdf')
plt.show()

xmid = np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])

fit = q.fit(xdata=xmid, ydata=n, model=maxwell_boltzmann, parguess=[1,1])

print(fit)

print(n)
print(bins)
