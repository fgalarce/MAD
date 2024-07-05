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

# Define the parser
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs =int(args.n_procs)

def loop(j):

  dataFile = "./results/data" + str(j)
  
  call("cp data " + dataFile, shell=True)
  
  setParameter("nbModes", j, dataFile) 
  
  if n_procs == 1:
      cmd = "./bound_pdrop.exe " + dataFile
  else:
      cmd = "./bound_pdrop.exe " + dataFile + " > ./results/" + str(j)  + ".log"
  
  print("=======================================")
  print(cmd)
  print("=======================================")
  call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(1,50):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(1,50)) 
