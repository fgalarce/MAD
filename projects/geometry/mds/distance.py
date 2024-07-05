from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput

def setParameter(parameter, value, fileName):
  datafile = open(fileName, 'r+')
  succes = False
  for line in fileinput.input( fileName ):
    if (line.split("=")[0] == parameter):
      datafile.write(line.replace(line[0:len(line)-1], line.split("=")[0] + "=" +  str(value)))
      succes = True
    else:
      datafile.write(line)
  if not succes:
    datafile.write(parameter +"=" +  str(value))
  datafile.close()

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

nbGeometries=16
mpi="/Home/flow/galarce/MAD/petsc/build/bin/mpirun"

for i in range(nbGeometries):
  for j in range(nbGeometries):
    if (i != j):
      dirResults="./results/" + str(i) + "_" + str(j) + "/"
      call("mkdir -p " + dirResults, shell=True)
      
      dataFile = dirResults + "par_distances"
      space_1="../../../../data/fpc/geo" + wildcard(i) + "/pod/velocity_mode"
      space_2="../../../../data/fpc/geo" + wildcard(i) + "/pod/piola/" + wildcard(j) + "/I_Tv"
      
      call("cp par_distances " + dataFile, shell=True)
      setParameter("dirResults", dirResults, dataFile) 
      setParameter("space_1", space_1 ,dataFile)
      setParameter("space_2", space_2 ,dataFile)
      cmd = mpi + " -np 28 ./distance.exe " + dataFile

      print(cmd)
      call(cmd, shell=True)
      call("rm " + dirResults + "/distances.geo", shell=True)
      call("rm " + dirResults + "/distances.case", shell=True)
