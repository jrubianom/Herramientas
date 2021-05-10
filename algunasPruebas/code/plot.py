import matplotlib.pyplot as plt
import numpy as np
import sys

N, Dir0, Ddir0, Dir3, Ddir3, Bloc0, Dbloc0, Bloc3, Dbloc3, Eig0, Deig0, Eig3, Deig3, Arm0, Darm0, Arm3, Darm3 = np.loadtxt(sys.argv[1], unpack=True)

#N, Dir0, Ddir0 = np.loadtxt(sys.argv[1], unpack=True)
#N1=[]


plt.errorbar(N, Dir0, yerr=Ddir0, label='Directo', fmt='-o' ,uplims=False)
#plt.errorbar(N, Dir3, yerr=Ddir3, label='Directo', fmt='-o' ,uplims=False)
plt.errorbar(N, Bloc0, yerr=Dbloc0, label='Blocking',fmt='-o', uplims=False)
#plt.errorbar(N, Bloc3, yerr=Dbloc3, label='Blocking', fmt='-o', uplims=False)
plt.errorbar(N, Eig0, yerr=Deig0, label='Eigen',fmt='-o', uplims=False)
#plt.errorbar(N, Eig3, yerr=Deig3, label='Eigen', fmt='-o',  uplims=False)
plt.errorbar(N, Arm0, yerr=Darm0, label='Armadillo', fmt='-o', uplims=False)
#plt.errorbar(N, Arm3, yerr=Darm3, label='Armadillo', fmt='-o', uplims=False)

plt.ylabel('Time')
plt.xlabel('N')
plt.yscale('log')
plt.grid()
plt.legend()
plt.savefig('Mult_Time_O0.pdf')
#plt.show()
