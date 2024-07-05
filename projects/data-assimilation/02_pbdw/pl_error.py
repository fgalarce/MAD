import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import scienceplots
import os
dir_results = os.environ['MAD_RESULTS']


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


pod = np.loadtxt("/Users/fgalarce/home/01_research/simulations/pod_windows/pod_00923_00933/eigenValues.txt")
beta = np.loadtxt('./carotid/G_sv.txt', skiprows=0);
error = np.zeros(100)

for i in range(100):
  error[i] = np.sqrt(np.sum(pod[i:]) / np.sum(pod))

print(error/beta)

#pl.style.use(['science', 'nature'])
pl.style.use('dark_background')
pl.plot(np.arange(1,101), error/beta, label='$\\beta(V_n, W_m)^{-1} dist(u,V_n)$')
pl.plot(np.arange(1,101), beta, label='$\\beta(V_n, W_m)$')
pl.plot(np.arange(1,101), error, label='$dist(u,V_n)$')
pl.xlabel("$n$")
pl.yscale("log")
pl.legend(loc='best')
#pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.tight_layout()
#pl.grid("on")
pl.show()
