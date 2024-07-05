import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
from numpy import linalg as la
import scienceplots
import pandas as pd
import os

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
timeSteps = 10000
dt=1e-3

dir_results = os.environ['MAD_RESULTS'] + "/hammer/"
pl.style.use(['science', 'nature'])
T = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
#T = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
pl.ylabel("$p^*$")
pl.xlabel("$t^*$")
p = np.loadtxt(dir_results + '/results/ctrlP.txt')[:,2];
q = np.loadtxt(dir_results + '/results2/ctrlP.txt')[:,2];
pl.plot(dt*np.arange(firstTime, timeSteps),p)
pl.plot(dt*np.arange(firstTime, timeSteps),q)

#pl.ylim([-0.8,1.0])
#pl.xlim([0.0,10.0])
#pl.xlim([0.0,40.0])
pl.tight_layout()
pl.legend(loc="upper right",fontsize=9)
#pl.savefig("u_vs_Bi_2.pdf", format="pdf", bbox_inches="tight")
pl.show()
