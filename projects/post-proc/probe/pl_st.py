import numpy as np
import fileinput
import matplotlib.pyplot as pl
from matplotlib import rc

p = np.loadtxt('./frecuency.txt', skiprows=0);

U_mean = (0.5 + 0.02*p[:,0])
Re = U_mean/0.001*0.1
St = p[:,1]*0.1/U_mean

print(U_mean)
print(Re)
print(St)

pl.plot(Re, St, '*')
pl.plot(Re, St)
pl.legend(loc='best')
pl.xlabel("Re")
pl.ylabel("St")
pl.show()
