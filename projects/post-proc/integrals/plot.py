import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc

def wildcard(idIter):
  if (idIter < 10):
   iteration = "0000" + str(idIter);
  elif (idIter < 100):
    iteration = "000" + str(idIter);
  elif (idIter < 1000):
    iteration = "00" + str(idIter);
  elif (idIter < 10000):
    iteration = "0" + str(idIter);
  elif (idIter < 100000):
    iteration = str(idIter);
  return iteration

fontSize=14
pl.figure(num=None, figsize=(4.5, 3))

p = np.loadtxt('./brain/p_int.txt', skiprows=0);
#q = np.loadtxt('./brain/p_int_rec.txt', skiprows=0);

pl.plot(p[1,:])
pl.grid()
#pl.legend(loc='best')
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$\\beta(V_n,W)$", fontsize=fontSize)
pl.tight_layout()
pl.show()
