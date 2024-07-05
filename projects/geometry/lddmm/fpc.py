from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput, sys
sys.path.insert(1, '../../../io/')
from tools import *

# Define the parser
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()
n_procs=int(args.n_procs)

first=0
nbMeshes=72
radius=arange(0.05, 0.1, 0.0007)

def loop(idMesh):
  for j in range(nbMeshes):
    ensName="lddmm_" + wildcard(j)
    dirResults="../../../../data/fpc/geo" + wildcard(idMesh) + "/mappings/"
    if (j != idMesh): 
      call("mkdir -p " + dirResults, shell=True)
    source = "../../../../data/mesh/fpc/fpc." + wildcard(idMesh) + ".mesh"
    target = "../../../../data/mesh/fpc/fpc." + wildcard(j) + ".mesh"

    if radius[j] < radius[idMesh]:
      run = "python3 LDDMM.py -scale=1 -save_folder=" + dirResults + " -name=" + ensName + " -source=" + source + " -target=" + target + " -dimension=2 -correct_fpc=1 -transversal_size=" + str(radius[idMesh])
      print(run)
    else:
#      run = "python3 LDDMM.py -scale=1 -save_folder=" + dirResults + " -name=" + ensName + " -source=" + source + " -target=" + target + " -dimension=2"
      run = "  "

    if n_procs == 1:
      cmd = run
    else:
      cmd = run + " > " + dirResults + "/lddmm" + wildcard(j) + ".log"
   
    if (j != idMesh): 
      print(cmd)
      call(cmd, shell=True) 

if n_procs == 1:
  for j in range(first, nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(first, nbMeshes)) 
