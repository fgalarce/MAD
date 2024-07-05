import numpy as np
import scienceplots
import matplotlib.pyplot as pl
from matplotlib import rc
import os

T = np.arange(0.0, 3, 0.25)
Bi = np.arange(0.0, 0.9, 0.1)
position = 1

pmaxArr = np.zeros((len(Bi),len(T)))
pl.style.use(['science', 'nature'])
dir_results = os.environ['MAD_RESULTS'] + "/hammer/plot_paper/contour_Bi_T/"

for i in range(len(Bi)):
  for j in range(len(T)):
    dirResults = dir_results + "/resultsTclosure_" + str(int(100*T[j])) + "/Bi_" + str(int(Bi[i]*10)) + "/"
    p = np.loadtxt(dirResults + 'ctrlP.txt', skiprows=0);
    pmax = max(abs(p[:,position]))
    pmaxArr[i,j] = pmax

corner_masks = [False, True]
x,y = np.meshgrid(T, Bi)
cs = pl.contourf(x, y, pmaxArr)
pl.colorbar()
pl.xlabel("$T_c^*$")
pl.ylabel("Bi")
pl.tight_layout()
pl.show()
