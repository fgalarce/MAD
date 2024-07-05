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

figure(num=None, figsize=(10, 6))
rcParams.update({'font.size': 16})

p = loadtxt('./pdrop_convergence.txt', skiprows=0);
fileName = 'pdrop.pdf'

v = arange(0,4)
v[0] = 2
v[1] = 3
v[2] = 4
v[3] = 5
theo_pdrop = arange(0,4)
theo_pdrop[0] = 3
theo_pdrop[1] = 3
theo_pdrop[2] = 3
theo_pdrop[3] = 3
print(v)
plot(p[:,0]/5.0, label='Monolithic', linewidth=2.0)
plot(p[:,1]/5.0, label='CT', linewidth=2.0)
plot(theo_pdrop, label='Theoretical', linewidth=2.0, linestyle='--')
legend(loc='best')
grid()
xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax = axes()
ax.set_xlabel("Thousands of vertices")
ax.set_ylabel("Pressure jump")
savefig(fileName)
clf()

call('open ' + fileName, shell=True)
