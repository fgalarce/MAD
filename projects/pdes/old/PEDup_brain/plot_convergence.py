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

p = np.loadtxt('./results_inclusion_3kv/norm_solution.txt', skiprows=0);
q = np.loadtxt('./results_inclusion_13kv/norm_solution.txt', skiprows=0);
r = np.loadtxt('./results_inclusion_47kv/norm_solution.txt', skiprows=0);

fontSize=16
pl.figure(num=None, figsize=(6, 4))
pl.rcParams.update({'font.size': fontSize})
pl.plot(p[:,0], label='3kv', linewidth=2.0)
pl.plot(q[:,0], label='13kv', linewidth=2.0)
pl.plot(r[:,0], label='47kv', linewidth=2.0)
pl.legend(loc='best')
pl.xlabel("Time", fontsize=fontSize)
pl.ylabel("Norm solution", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()
