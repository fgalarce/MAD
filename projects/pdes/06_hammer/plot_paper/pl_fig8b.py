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

c = 400
#B = 0.5
l = 110
timeSteps=10000
time = np.arange(timeSteps)*1e-3
average_withoutConvection = np.zeros(timeSteps)
average_withConvection = np.zeros(timeSteps)
envolvente_max_conv = np.zeros(timeSteps)
envolvente_min_conv = np.zeros(timeSteps)
envolvente_max_without_conv = np.zeros(timeSteps)
envolvente_min_without_conv = np.zeros(timeSteps)
for B in Bi:
  dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
  p = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
  average_withConvection = average_withConvection + p
  dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
  q = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
  average_withoutConvection = average_withoutConvection + q
#  pl.plot(time, p*p, "r", alpha=0.2)
#  pl.plot(time, q*q, "k", alpha=0.2)

B=0.0
dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
p = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
average_withConvection = average_withConvection + p
dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
q = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
#pl.plot(time, p*p, "r")
#pl.plot(time, q*q, "k")
pl.plot(time, abs(p*p - q*q), label="$\\text{Bi}=0.0$")

#B=0.4
#dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
#p = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
#average_withConvection = average_withConvection + p
#dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
#q = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
##pl.plot(time, p*p, "r")
##pl.plot(time, q*q, "k")
#pl.plot(time, abs(p*p - q*q), label="$\\text{Bi}=0.4$")
#
#B=0.6
#dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
#p = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
#average_withConvection = average_withConvection + p
#dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
#q = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
##pl.plot(time, p*p, "r")
##pl.plot(time, q*q, "k")
#pl.plot(time, abs(p*p - q*q), label="$\\text{Bi}=0.6$")

B=0.8
dirResults=dir_results + "/without_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
p = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
average_withConvection = average_withConvection + p
dirResults=dir_results + "/with_convection/Lenght_" + str(int(l)) + "/Bi_" + str(int(B*10)) + "/c_" + str(int(c)) + "/"
q = np.loadtxt(dirResults + 'normSolution.txt', skiprows=0);
#pl.plot(time, p*p, "r*")
#pl.plot(time, q*q, "k*")
pl.plot(time, abs(p*p - q*q), label="$\\text{Bi}=0.8$")

average_withConvection = average_withConvection / float(len(Bi))
average_withoutConvection = average_withoutConvection / float(len(Bi))
p = average_withConvection
q = average_withoutConvection
#pl.plot(time, average_withConvection * average_withConvection, "r--", label="convective", linewidth=2.0)
#pl.plot(time, average_withoutConvection * average_withoutConvection, "k--", label="$\\lambda_2 = 0$")
pl.plot(time, abs(p*p - q*q), label="Average")
pl.ylabel("Pseudo Kinetic Energy Discrepancy")
pl.xlabel("$t^*$")
pl.legend()
pl.tight_layout()
#pl.yscale("log")
pl.show()
