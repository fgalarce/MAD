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
pmax_withConv = np.zeros(len(Lenght))
pmax_withoutConv = np.zeros(len(Lenght))
l = 110
c = 400
B = 0.6
position = 2
for i in range(len(Lenght)):
  dirResults=dir_results + "/with_convection/Lenght_" + str(int(Lenght[i])) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
  p = np.loadtxt(dirResults + 'CtrlP.txt', skiprows=0)[:,position];
  pmax_withConv[i] = max(abs(p))
  dirResults=dir_results + "/without_convection/Lenght_" + str(int(Lenght[i])) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
  q = np.loadtxt(dirResults + 'CtrlP.txt', skiprows=0)[:,position];
  pmax_withoutConv[i] = max(abs(q))

pl.plot(Lenght, pmax_withConv, "r--", label="convective", linewidth=2.0)
pl.plot(Lenght, pmax_withoutConv, "k--", label="$\\lambda_2 = 0$", linewidth=2.0)
pl.ylabel("Adimensional Pressure Peak")
pl.xlabel("Adimensional Pipe Lenght")
pl.legend(loc='best')
pl.tight_layout()
#pl.savefig("pmax_conv_Bi" + str(B*10) + ".pdf", format="pdf", bbox_inches="tight")
pl.show()

