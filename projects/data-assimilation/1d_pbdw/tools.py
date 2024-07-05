from subprocess import *
import matplotlib.pyplot as pl
import numpy as np
import fileinput
import seaborn as sns
sns.set_theme()
fontSize=18
pl.rcParams.update({'font.family': 'serif'})
pl.rcParams.update({'font.size': fontSize})

EPSILON = 1e-6

#def pde(parameters):

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

def orthonormalize(A):
  for i in range(np.shape(A)[1]):
    for j in range(i):
      A[:,j] = A[:,j] - proj(A[:,i], A[:,j])
    A[:,i] = A[:,i] / np.sqrt(np.dot(A[:,i],A[:,i]))
  return A
 
def proj(u,v):
  return np.dot(u,v)*u / np.sqrt(np.dot(u,u))


def fastplot(x, y, label_x=" ", label_y=" ", label=" "):
  fontSize=16
  pl.figure(num=None, figsize=(12, 8))
  pl.plot(x, y, label=label)
  if (label_x != " "):
    pl.xlabel("x", fontsize=fontSize)
  if (label_y != " "):
    pl.ylabel("eta", fontsize=fontSize)
  if (label != " "):
    pl.legend(loc='best')
  pl.tight_layout()
  pl.grid()
  pl.show()

def show(legend=0):
  if (legend == 1):
    pl.legend(loc='best')
  pl.tight_layout()
  pl.grid()
  pl.show()

def isinside(number, array):
  for i in range(len(array)):
    if (number == array[i]):
      return 1
  return 0

