import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import seaborn as sns
sns.set_theme()
fontSize=18
pl.rcParams.update({'font.family': 'serif'})
pl.rcParams.update({'font.size': fontSize})

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


m = np.arange(6,27,2)
error = np.array([0.20676335064793608, 0.17885166273317799, 0.08040996111337173, 0.06354313299229775, 0.04264982603620565, 0.03142607684911316, 0.01809169572807876, 0.015043809056288828, 0.013972933765276963, 0.010976618953602413, 0.010047394855737095])

ancho_fig = 6
alto_fig = 4
pl.figure(num=None, figsize=(ancho_fig, alto_fig))
pl.plot(m, error, linewidth=2.0)
pl.legend(loc='best')
pl.xlabel("$m$")
pl.ylabel("max. error")
pl.yscale("log")
pl.grid("on")
pl.tight_layout()
pl.show()
