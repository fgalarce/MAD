import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc

fontSize=16
pl.figure(num=None, figsize=(6, 4))
pl.rcParams.update({'font.size': fontSize})
p = np.loadtxt('./carotid2/windkessel.txt', skiprows=1);

first=20

pl.plot(p[first:len(p[:,0]),0])
pl.plot(p[first:len(p[:,0]),1])
pl.plot(p[first:len(p[:,0]),2])
pl.plot(p[first:len(p[:,0]),3])
pl.grid()

pl.show()
pl.savefig("windk.pdf")
