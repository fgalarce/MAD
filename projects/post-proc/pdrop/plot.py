from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *

from matplotlib import rc
rc('text', usetex=True)

style.use('ggplot')

figure(num=None, figsize=(10, 6))
rcParams.update({'font.size': 16})

p = loadtxt('./results/p_drop.txt', skiprows=0);
q = loadtxt('./results/p_drop_GT.txt', skiprows=0);
fileName = 'imrp.pdf'

plot(p[:,0] * 760.0 / ( 10.0 * 101325.0), label='$\Gamma_o^1$', linewidth=3.0)
plot(p[:,1] * 760.0 / ( 10.0 * 101325.0), label='$\Gamma_o^2$', linewidth=3.0)
#plot(q * 760.0 / ( 10.0 * 101325.0), label='gt', linewidth=3.0)
legend(loc='best')
grid()
xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax = axes()
ax.set_xlabel("Time it")
ax.set_ylabel("Pressure drop [mmHg]")
savefig(fileName)
clf()

call('open ' + fileName, shell=True)

import fileinput

def getFromDataFile(dataFilePath, variableName):
  datafile = open(dataFilePath, 'r')
  line = datafile.readline()
  while line:
    if (line.split(" = ")[0] == variableName):
      return line.split("=")[1]
    line = datafile.readline()

# search for line top and replace it by bottom
def fileReplace(top, bottom, fileName):
  tempFile =  open(fileName, 'r+' )
  for line in fileinput.input( fileName ):
      tempFile.write( line.replace(top, bottom)  )
  tempFile.close()

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
