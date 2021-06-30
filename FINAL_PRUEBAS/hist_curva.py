import matplotlib.pyplot as plt
import numpy as np
import qexpy as q
import sys
import qexpy.plotting as qplt

def maxwell_boltzmann(v,a,b):
    return a * v * q.exp(-b*v**2)

vel = np.loadtxt("velocities_data.txt", unpack=True)
#Nbins = 10

n,bins,figure = qplt.hist(vel,density=True,edgecolor='black')
xmid = np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])

#fit = q.fit(xdata=xmid, ydata=n, model=maxwell_boltzmann, parguess=[1,1])
fit = figure.fit(model=maxwell_boltzmann, parguess=[1,1])
figure.show()
print(fit)

print(n)
print(bins)
