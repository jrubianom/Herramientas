import matplotlib.pyplot as plt
import numpy as np
import sys


for i in range(6):
    name = "velocities_data"+str(i)+"-100.txt"
    vel = np.loadtxt(name, unpack=True)
    plt.hist(vel,density=True,edgecolor='black')
    name_hist = "Histograma"+str(i)+"-100.pdf"
    plt.savefig(name_hist)
