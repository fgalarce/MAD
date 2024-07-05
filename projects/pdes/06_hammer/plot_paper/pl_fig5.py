import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
from numpy import linalg as la
import scienceplots
import pandas as pd

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

firstTime = 0
timeSteps = 50000
dt=1e-3

pl.style.use(['science', 'nature'])
#Bi = np.arange(0.0, 1.1, 0.2)
Bi = np.array([0.0, 0.2, 0.6, 1.0, 1.4])
pl.ylabel("$u^*$")
pl.xlabel("$t^*$")
for i in range(len(Bi)):
  p = np.loadtxt('./resultsBi_' + str(int(Bi[i]*100)) + '/ctrlU.txt');
  pl.plot(dt*np.arange(firstTime, timeSteps), p[firstTime:timeSteps,0], label="$Bi=$" + " %.1f" % (Bi[i]))

pl.ylim([-0.8,1.0])
#pl.ylim([-30,40.0])
pl.tight_layout()
pl.legend(loc="upper right",fontsize=9)
pl.savefig("u_vs_Bi_2.pdf", format="pdf", bbox_inches="tight")
pl.show()
