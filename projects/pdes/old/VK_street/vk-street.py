from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput

def getParameter(dataFilePath, variableName):
  datafile = open(dataFilePath, 'r')
  line = datafile.readline()
  while line:
    if (line.split(" = ")[0] == variableName):
      return line.split("=")[1]
    line = datafile.readline()

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

# Define the parser
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs = int(args.n_procs)

nbMeshes = 8;
nbSimulations = 32;
Re = arange(300, 299, -3.125)

diameters_file = loadtxt('../../../../big_data/vk-street/mesh/ellipse_info.txt')
diameters = arange(nbMeshes)
for i in range(nbMeshes):
  if (diameters_file[i][1] > diameters_file[i][2]):
    diameters[i] = 2.0 * diameters_file[i][1]
  else:
    diameters[i] = 2.0 * diameters_file[i][2]

def loop(j):

  for i in range(nbSimulations):

    dirResults='../../../../big_data/vk-street/data_base_/mesh' + wildcard(j) + '/sim' + wildcard(i) + '/'
    call("mkdir -p " + dirResults, shell=True)
    
    dataFile = dirResults + "parameters"
    
    call("cp parameters " + dataFile, shell=True)

    u0 = 5.0
    viscosity = diameters[j] * u0 / Re[i]
    
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("geometryData", '../../../../big_data/vk-street/mesh/ellipse' + wildcard(j) + ".mesh", dataFile) 
    setParameter("viscosity", viscosity, dataFile) 
   
     
    if n_procs == 1:
        cmd = "/local/fgalarce/MAD/mad/external_software/petsc/build/bin/mpirun -np 30 ./ns.exe " + dataFile
    else:
        cmd = "/local/fgalarce/MAD/mad/external_software/petsc/build/bin/mpirun -np 10 ./ns.exe " + dataFile + " > " + dirResults + "/vk-street.log"
    
    print(cmd)
    call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(4,nbMeshes)) 
