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
dirResults = os.environ['MAD_RESULTS'] + "/hammer/results/"

x = np.arange(0, 1+dx, dx)
#line, = ax.plot(x, np.loadtxt(dirResults + 'sol.' + wildcard(0) + '.txt')[nbVertices:2*nbVertices+1])
line, = ax.plot(x, np.loadtxt(dirResults + 'sol.' + wildcard(0) + '.txt')[0:nbVertices])
#ax.set_ylim([-1.2,1.2])
#ax.set_ylim([-0.5,0.5])
ax.set_ylim([-20.0,20.0])

def animate(i):
    p = np.loadtxt(dirResults + 'sol.' + wildcard(i) + '.txt');
    line.set_ydata(p[nbVertices:nbVertices*2]) 
#    line.set_ydata(p[0:nbVertices]) 
    return line,

ani = animation.FuncAnimation(fig, animate, interval=2, blit=True, save_count=50, frames=np.arange(timeSteps))
plt.show()
