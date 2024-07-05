from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import *
from subprocess import *
from matplotlib import cm
from numpy import *

style.use('ggplot')

fig = figure(num=None, figsize=(10, 6))
ax = fig.gca(projection='3d')
rcParams.update({'font.size': 16})

p = loadtxt('./Y.txt', skiprows=0);
fileName = 'Y.pdf'

ax.scatter(p[:,0], p[:,1], p[:,2], marker='o')


xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax.set_xlabel("$y_1$")
ax.set_ylabel("$y_2$")
ax.set_zlabel("$y_3$")
#plt.savefig(fileName)
show()
#clf()

#call('open ' + fileName, shell=True)
