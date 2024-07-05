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

nbMeshes = 60

i = 61

e_learned = zeros((nbMeshes))
p = loadtxt('../../data/venturi/mesh/v' + wildcard(i) + "/g.txt", skiprows=0);
print('Geometry ' + str(i))
for j in range(nbMeshes):
  q = loadtxt('./results/weights'+ wildcard(j) + '.txt');
  e_learned[j] = dot(p,q)

figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
rcParams.update({'font.size': 20})
latexSize=20
tickSize=20
style.use('ggplot')

plot(e_learned, linewidth=3.0)
#yscale('log')
legend(loc='best')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
#ax.set_xlabel("<+xlabel+>", fontsize=latexSize)
#ax.set_ylabel("<+ylabel+>", fontsize=latexSize)
savefig('./e_learned.pdf')
clf()

# box plot
#boxplot(<+data+>, 0, '') # disable outliers
#yscale('log')
#xlabel(<+xlabel+>, fontsize=latexSize)
#ylabel(<+ylabel+>, fontsize=latexSize)
#axes = gca()
#axes.set_ylim([0.0001, 0.1])
#grid()
#xticks(arange(1,8,1), arange(10,80,10), fontsize=tickSize)
#yticks(fontsize=tickSize, rotation = 45)
#savefig(fileName)
#clf()

#call('evince ' + fileName, shell=True)
