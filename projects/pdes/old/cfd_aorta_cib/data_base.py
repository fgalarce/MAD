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

#alpha = arange(0.9, 1.15, 0.05)
alpha = random.uniform(0.9, 1.15, 128)
i = 0

for u0 in alpha:
  for Rd in alpha:

    dirResults='./phantom_leo4/sim' + wildcard(i) + '/'
    call("mkdir -p " + dirResults, shell=True)
    
    dataFile = dirResults + "parameters"
    
    call("cp parameters " + dataFile, shell=True)
    
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("inlet_u0", u0, dataFile) 
    setParameter("distalResistances", str(8000.0*Rd) + " " + str(8000.0*Rd) + " " + str(8000.0*Rd) + " " + str(1000.0*Rd), dataFile)
     
    cmd = "/Home/flow/galarce/MAD/mad/external_software/petsc/build/bin/mpirun -np 36 ./cib.exe " + dataFile + " > " + dirResults + "/cib.log"
    
    print(cmd)
    call(cmd, shell=True) 
    i = i + 1
