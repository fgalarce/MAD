from subprocess import *
from numpy import *
from sys import *
from multiprocessing import Pool
import fileinput
import argparse

parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs = int(args.n_procs)

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
  
def loop_time(i):
  for j in range(5):

    dirResults="../../data/model/pod_mixed_linear_H1/window_time" + str(i) + "_HR" + str(j) + "/"

    call("mkdir -p " + dirResults, shell=True)
    dataFile = dirResults + "data"
    call("cp data_G " + dataFile, shell=True)

    setParameter("dirResults", dirResults, dataFile)

    if n_procs == 1:
      cmd = "./computeG.exe " + dataFile
    else:
      cmd = "./computeG.exe " + dataFile + " > " + dirResults + "or.log"

    print(cmd)
    call(cmd, shell=True)
    
if n_procs == 1:
  for j in range(5):
    loop_time(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop_time, range(5))
