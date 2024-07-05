import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc

p = np.loadtxt("./coordinates_brain.txt")

x_prime = p[:,2]
y_prime = p[:,1]
z_prime = -p[:,0]

q = np.zeros(np.shape(p))

q[:,0] = x_prime
q[:,1] = y_prime
q[:,2] = z_prime
q[:,3] = p[:,3]

print(q)

np.savetxt("./coordinates_brain_rotated.txt", q)
