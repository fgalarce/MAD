from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *
from matplotlib import rc

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

rc('text', usetex=True)
style.use('ggplot')

fontSize=25
tickSize=20
figure(num=None, figsize=(9, 6))
rcParams.update({'font.size': fontSize})

p = loadtxt('./Y.m');

for i in range(16):
  if (i < 7):
    plot(p[0, i], p[1,i], 'o', label=str(i), linewidth=2.0)
  elif (i < 14 and i >= 7):
    plot(p[0, i], p[1,i], 'x', label=str(i), linewidth=2.0)
  else:
    plot(p[0, i], p[1,i], '.', label=str(i), linewidth=2.0)

legend(loc='best',fontsize=12)
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
show()
#ax.set_xlabel("<++>", fontsize=fontSize)
#ax.set_ylabel("<++>", fontsize=fontSize)
#yscale("log")
