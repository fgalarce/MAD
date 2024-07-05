from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *

from matplotlib import rc
rc('text', usetex=True)

figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
rcParams.update({'font.size': 20})
latexSize=20
tickSize=20
style.use('ggplot')

p = loadtxt('./results_noise_n_star_cls/error_pbdw.txt', skiprows=1);
q = loadtxt('./results_noise_n_star/error_pbdw.txt', skiprows=1);
fileName = 'cls_nstar.pdf'

plot(p[:,0]/p[:,1], label='cls', linewidth=3.0)
plot(q[:,0]/q[:,1], label='free', linewidth=3.0)
#yscale('log')
legend(loc='best')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
ax.set_xlabel("Time", fontsize=latexSize)
ax.set_ylabel("Error", fontsize=latexSize)
savefig(fileName)
clf()

call('evince ' + fileName, shell=True)
