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

Bi = np.arange(0.0, 0.9, 0.1)
Lenght = np.arange(50,220, 20)
wave_velocity = np.arange(400,1300,100)

nbSnapshots = len(Bi) * len(Lenght) * len(wave_velocity)

#c = 400
#B = 0.5
#l=210
timeSteps=10000
time = np.arange(timeSteps)*1e-3
average_withoutConvection = np.zeros(timeSteps)
average_withConvection = np.zeros(timeSteps)
for l in Lenght:
  for c in wave_velocity:
    for B in Bi:
      dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
      p = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
      average_withConvection = average_withConvection + p
      dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
      q = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
      average_withoutConvection = average_withoutConvection + q
#      pl.plot(time, p, "r", alpha=0.5)
#      pl.plot(time, q, "k", alpha=0.5)

#average_withConvection = average_withConvection / float(len(Lenght))
#average_withoutConvection = average_withoutConvection / float(len(Lenght))
pl.plot(time, average_withConvection*average_withConvection, "r--", label="convective", linewidth=2.0)
pl.plot(time, average_withoutConvection*average_withoutConvection, "k--", label="$\\lambda_2 = 0$", linewidth=2.0)
pl.ylabel("Average velocity norm")
pl.xlabel("$t^*$")
pl.yscale("log")
pl.show()
