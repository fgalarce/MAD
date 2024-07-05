from subprocess import *
import numpy as np
from multiprocessing import Pool
import fileinput
import sys
import mesh_io as io

dirName = '/Users/galarce/research/mad/source/projects/other/field_element/'
filenameMesh = './mesh/brain_alfonso/brain_3d_fine.mesh'
#boundaries = np.array(['MRE_band.csv', 'neck.csv'])
#idLabel=np.array([1,2])
boundaries = np.array(['el_pulso_bien.csv', 'el_cuello_bien.csv'])
idLabel=np.array([1,2])

vertices, tria, tetrahedra = io.importInriaMesh(filenameMesh)

#tria = np.zeros((len(tria[:,0]),4))
tria[:,3] = 0.0

for i in range(len(idLabel)):
  trianglesBd = np.loadtxt(dirName + boundaries[i], skiprows=1)
  print(trianglesBd)
  tria[trianglesBd.astype(int), 3] = idLabel[i]

np.savetxt('tria_labeled', tria, fmt='%i')
