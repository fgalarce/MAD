import numpy as np
import matplotlib.pyplot as plt

# Outlets
outlets = [3, 5, 6, 7]
legends = ['Outlet {:d}'.format(i) for i in outlets]

# Import windkessels
wk = np.loadtxt("wkssl_pressure.txt")

# Convert to mmHg
mmHg2cgs = 1333.22
cgs2mmHg = 1.0/mmHg2cgs

# Generate time array
dt = 0.000937
t = np.linspace(dt, dt*(wk.shape[0] + 1), wk.shape[0])

# Plot pressures
plt.figure(1)
for i in range(wk.shape[1]):
  plt.plot(t, cgs2mmHg*wk[:,i])
plt.xlabel('time (s)')
plt.ylabel('proximal pressure (mmHg)')
plt.legend(legends)
plt.show()
