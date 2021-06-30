import matplotlib.pyplot as plt
import numpy as np
import sys


vel = np.loadtxt(sys.argv[1], unpack=True)
plt.hist(vel,edgecolor='black')
plt.xlabel('Velocidad')
plt.ylabel('Frecuencia')
plt.savefig('ULTIMO.pdf')
