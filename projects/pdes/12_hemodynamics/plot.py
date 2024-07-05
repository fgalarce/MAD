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

fontSize=16
pl.figure(num=None, figsize=(6, 4))
pl.rcParams.update({'font.size': fontSize})
p = np.loadtxt('./carotid2/windkessel.txt', skiprows=1);

fig, axs = pl.subplots(1,len(p[0,:]))
#fig.suptitle('Vertically stacked subplots')
for i in range(len(p[0,:])):
  axs[i].plot(p[120:len(p[:,i]),i])

pl.show()
pl.savefig("windk.pdf")
