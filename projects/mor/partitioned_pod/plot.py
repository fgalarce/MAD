from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *

from matplotlib import rc
rc('text', usetex=True)

style.use('ggplot')

figure(num=None, figsize=(10, 6))
rcParams.update({'font.size': 16})

p = loadtxt('./results/eigenValues.txt', skiprows=0);
q = loadtxt('./results_mdma/eigenValues.txt', skiprows=0);

nbModes = p.size

error_felisce = zeros((nbModes,))
error_mdma = zeros((nbModes,))
for i in range(nbModes):
    error_felisce[i] = sqrt(sum(p[i:nbModes]) / sum(p))
    error_mdma[i] = sqrt(sum(q[i:nbModes]) / sum(q))

fileName = 'ev.pdf'

plot(error_felisce, label='felisce', linewidth=2.0)
plot(error_mdma, label='mdma', linewidth=3.0, linestyle='-.')
legend(loc='best')
grid()
yscale("log")
xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax = axes()
ax.set_xlabel("$n$")
ax.set_ylabel("$dist(u, V_n)$")
savefig(fileName)
clf()

call('open ' + fileName, shell=True)
