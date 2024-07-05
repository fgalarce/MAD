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

n = 20
p = zeros(n-1)
for j in range(25):
  for i in range(1,n):
    p[i-1] = loadtxt('./results/mu' + str(i) + '.txt', skiprows=0)[j];
  plot(range(1,n), p, linewidth=1.5)

for i in range(5):
  for j in range(5):
    ev = loadtxt('/local/fgalarce/4d-flow-ultrasound/data/model/pod_mixed_linear_H1/window_time' + str(i) + '_HR' + str(j) + '/window_time' + str(i) + '_HR' + str(j) + "/eigenValues.txt")
    plot()

fileName = 'bound.pdf'
#yscale('log')
legend(loc='best')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
ax.set_xlabel("$n$", fontsize=latexSize)
ax.set_ylabel("$2 \epsilon_n \mu $", fontsize=latexSize)
savefig(fileName)
clf()

call('evince ' + fileName, shell=True)
