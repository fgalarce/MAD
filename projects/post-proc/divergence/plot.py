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

p = loadtxt('./results/norm_div_u.txt', skiprows=0)[0:30];
q = loadtxt('./results_02/norm_div_u.txt', skiprows=0)[0:30];
r = loadtxt('./results_no_piola/norm_div_u.txt', skiprows=0)[0:30];
fileName = 'div.pdf'

time = arange(p.size)
time = time * 0.01
plot(time, p, linewidth=2.0, label='u $\Omega_i$')
plot(time, q, linewidth=2.0, label='u $\Omega_j$')
plot(time, q, linewidth=2.0, label='P u $\Omega_j$')
legend(loc='best')
grid()
xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax = axes()
ax.set_xlabel("Time [sec]")
ax.set_ylabel("L2 norm of velocity divergence")
savefig(fileName)
clf()

call('open ' + fileName, shell=True)
