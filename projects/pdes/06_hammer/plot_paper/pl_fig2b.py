import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import scienceplots
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

L = 200
a = 1000
rho = 2000
g = 9.81
h = 100
pr = rho * g * h

dir_results = os.environ['MAD_RESULTS'] + "/hammer/plot_paper/"
fft2 = np.loadtxt(dir_results + 'fft/fft_2.txt', skiprows=0)[0:100]*pr;
fft2_f = np.loadtxt(dir_results + 'fft/fft_2_f.txt', skiprows=0)[0:100];
fft3 = np.loadtxt(dir_results + 'fft/fft_3.txt', skiprows=0)[0:100]*pr;
fft3_f = np.loadtxt(dir_results + 'fft/fft_3_f.txt', skiprows=0)[0:100];
fft4 = np.loadtxt(dir_results + 'fft/fft_4.txt', skiprows=0)[0:100]*pr;
fft4_f = np.loadtxt(dir_results + 'fft/fft_4_f.txt', skiprows=0)[0:100];
fft5 = np.loadtxt(dir_results + 'fft/fft_4.txt', skiprows=0)[0:100]*pr;
fft5_f = np.loadtxt(dir_results + 'fft/fft_4_f.txt', skiprows=0)[0:100];

pl.style.use(['science', 'nature'])
pl.plot(fft2_f, fft2, alpha=0.3, label='$\Delta t=2 \\times 10^{-3}$ s')
pl.plot(fft3_f, fft3, alpha=0.5, label='$\Delta t=2 \\times 10^{-4}$ s')
pl.plot(fft4_f, fft4, label='$\Delta t=2 \\times 10^{-5}$ s')
pl.plot(fft5_f, fft5, "--k", alpha= 0.4, label='$\Delta t=2 \\times 10^{-6}$ s')
pl.xlabel("Frecuency [Hz]")
pl.ylabel("Amplitude [Pa]")
pl.yscale("log")
pl.legend()
pl.tight_layout()
pl.savefig("convergence_frecuency_zoom.pdf", format="pdf", bbox_inches="tight")
pl.show()
