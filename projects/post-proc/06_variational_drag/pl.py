
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

p = np.loadtxt('energy_k.txt', skiprows=500);

pl.style.use(['science', 'nature'])
pl.plot(p, label='12123')
pl.xlabel("")
pl.ylabel("")
#pl.yscale("log")
pl.legend(loc='best')
#pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#pl.tight_layout()
#pl.grid("on")
pl.show()
