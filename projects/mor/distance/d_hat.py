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

p = loadtxt('./dd.txt');
q = zeros((16,16));

for i in range(16):
  for j in range(16):
    q[i,j] = loadtxt('matrix_n5_haus/haus' + str(i) + "_" + str(j) + ".txt");

print(p)
print(q)

savetxt("dhatd.txt", abs(p - q))
    
