import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
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

plt.style.use(['science', 'nature'])

fig, ax = plt.subplots()

timeSteps = 50000
nbVertices=501
dx=0.002
length=200
space = np.arange(nbVertices)*length
dirResults1 = os.environ['MAD_RESULTS'] + "/hammer/results_withConv/"
dirResults2 = os.environ['MAD_RESULTS'] + "/hammer/results_withoutConv/"

x = np.arange(0, 1+dx, dx)
#line, = ax.plot(x, np.loadtxt(dirResults + 'sol.' + wildcard(0) + '.txt')[nbVertices:2*nbVertices+1])
#line, = ax.plot(x, np.loadtxt(dirResults1 + 'sol.' + wildcard(0) + '.txt')[0:nbVertices])
#ax.set_ylim([-1.2,1.2])
#ax.set_ylim([-0.5,0.5])
#ax.set_ylim([-20.0,20.0])

p = np.loadtxt(dirResults1 + 'sol.' + wildcard(100) + '.txt');
q = np.loadtxt(dirResults2 + 'sol.' + wildcard(100) + '.txt');
plt.plot(p, 'k', label='p')
plt.plot(q, 'k--', label='q')

p = np.loadtxt(dirResults1 + 'sol.' + wildcard(400) + '.txt');
q = np.loadtxt(dirResults2 + 'sol.' + wildcard(400) + '.txt');
plt.plot(p, 'r', label='p')
plt.plot(q, 'r--', label='q')

p = np.loadtxt(dirResults1 + 'sol.' + wildcard(1000) + '.txt');
q = np.loadtxt(dirResults2 + 'sol.' + wildcard(1000) + '.txt');
plt.plot(p, 'g', label='p')
plt.plot(q, 'g--', label='q')

p = np.loadtxt(dirResults1 + 'sol.' + wildcard(5000) + '.txt');
q = np.loadtxt(dirResults2 + 'sol.' + wildcard(5000) + '.txt');
plt.plot(p, 'b', label='p')
plt.plot(q, 'b--', label='q')

plt.show()
