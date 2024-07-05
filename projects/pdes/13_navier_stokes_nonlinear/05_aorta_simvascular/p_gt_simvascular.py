from numpy import *
c = 1.0/4/pi
rd = 1.0
rp = 0.1
q0 = 10

time = arange(0,1.0,0.001)
P = zeros(len(time))
for i in range(len(time)):
  t_star = time[i]/rd/c
  P[i] = (rp/rd + 0.5) * ((sin(t_star/2))**2) + 0.25 * (1.0 - exp(-t_star) - sin(t_star))
  P[i] = P[i] * rd * q0

import matplotlib.pyplot as pl
pl.plot(P)
pl.show()
