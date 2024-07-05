from subprocess import *
import matplotlib.pyplot as pl
import numpy as np
from multiprocessing import Pool
import fileinput
import scipy
from scipy import linalg
import seaborn as sns
sns.set_theme()
fontSize=18
pl.rcParams.update({'font.family': 'serif'})
pl.rcParams.update({'font.size': fontSize})

EPSILON = 1e-6

def getCenter(meshF2, fieldF2):
  x_start = 0.0
  i_start = 0
  for i in range(len(meshF2)):
    if (fieldF2[i] > 0):
      x_start = meshF2[i]
      i_start = i
      break

  x_end = 0.0
  for i in range(i_start, len(meshF2)):
    if (fieldF2[i] == 0):
      x_end = meshF2[i]
      i_end = i
      break;
  dirac = np.zeros(len(meshF2))
  dirac[i_start] = 1.0
  dirac[i_end] = 1.0
  center = x_start + (x_end - x_start)/2.0
  i_center = i_end - i_start
#  pl.plot(meshF2, dirac)
  return center

def getIcenter(meshF2, fieldF2):
  x_start = 0.0
  i_start = 0
  for i in range(len(meshF2)):
    if (fieldF2[i] > 0):
      x_start = meshF2[i]
      i_start = i
      break

  x_end = 0.0
  for i in range(i_start, len(meshF2)):
    if (fieldF2[i] == 0):
      x_end = meshF2[i]
      i_end = i
      break;
  dirac = np.zeros(len(meshF2))
  dirac[i_start] = 1.0
  dirac[i_end] = 1.0
  center = x_start + (x_end - x_start)/2.0
  i_center = int((i_end - i_start)/2.0)
#  pl.plot(meshF2, dirac)
  return i_start
    
def is_ortho(A):
  isortho = 1
  for i in range(len(A[0,:])):
    for j in range(len(A[0,:])):
      dp = np.dot(A[:,i], A[:,j])
      if (dp < - EPSILON and dp > EPSILON and i != j):
        isortho = 0
  return isortho

def is_normal(A):
  isnormal = 1
  for i in range(len(A[0,:])):
    dp = np.dot(A[:,i], A[:,i])
    if (dp < 1 - EPSILON and dp > 1 + EPSILON):
      isnormal = 0
  return isnormal

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

def measure_empirical():

  return

