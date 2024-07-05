import numpy as np
import scienceplots
import matplotlib.pyplot as pl
from matplotlib import rc
import os

T = np.arange(0.0, 3, 0.25)
Bi = np.arange(0.0, 0.9, 0.1)
position = 2

pmaxArr = np.zeros((len(Bi),len(T)))
pl.style.use(['science', 'nature'])
dir_results = os.environ['MAD_RESULTS'] + "/hammer/picard/"

Lenght = np.arange(50,220, 20)
timeSteps=10000
time = np.arange(timeSteps)*1e-3
l = 190
c = 1000
B = 0.0
position = 2
dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
p = np.loadtxt(dirResults + 'CtrlP.txt', skiprows=0)[:,position];
dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
q = np.loadtxt(dirResults + 'CtrlP.txt', skiprows=0)[:,position];

pl.plot(time, p, "r--", label="convective", linewidth=2.0)
pl.plot(time, q, "k--", label="$\\lambda_2 = 0$", linewidth=2.0)
pl.ylabel("$p^*$")
pl.xlabel("$t^*$")
pl.legend(loc='best')
pl.tight_layout()
pl.savefig("pmax_conv_L" + str(l) + ".pdf", format="pdf", bbox_inches="tight")
pl.show()

