import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc

p = np.loadtxt('./results/sim00030/pressure_at_probe.txt', skiprows=0);

dt = 0.01

pl.plot(np.arange(0,100)*dt, p)
pl.xlabel("Time [sec]")
pl.ylabel("Pressure [Pa]")
pl.show()
