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


mu = 0
sigma = 0.05
t = np.arange(0,np.pi,0.01)
signal = np.sin(t)

noise = np.zeros(len(t))
alpha = 0.2
# noise generator
for i in range(len(t)):
  noise[i] =  np.random.normal(-signal[i]/np.max(signal)*alpha, sigma, 1)

fontSize=16
pl.figure(num=None, figsize=(6, 4))
pl.rcParams.update({'font.size': fontSize})
#pl.plot(x, label='<++>', linewidth=2.0)
pl.plot(t, signal, label='signal', linewidth=2.0)
pl.plot(t, signal + noise, label='noisy signal', linewidth=1.0)
pl.legend(loc='best')
#pl.xlabel("<++>", fontsize=fontSize)
#pl.ylabel("<++>", fontsize=fontSize)
#pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()
