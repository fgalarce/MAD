from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *
import fileinput

from matplotlib import rc
rc('text', usetex=True)
figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
rcParams.update({'font.size': 30})
latexSize=35
tickSize=30
style.use('ggplot')

curves= zeros((150, 10))
it = 0
for i in range(5):
  for j in range(2):
    p = loadtxt('./window_time' + str(i) + "_HR" + str(j) + "/eigenValues.txt")
    for n in range(150): 
      curves[n, it] = sqrt( sum(p[n:p.size])/p[0] )
    it = it + 1
    
plot(curves, linewidth=2.0)
grid()
yscale('log')
xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax = axes()
ax.set_xlabel("n")
ax.set_ylabel("Error")
savefig("ev_u.pdf")
