import random
import numpy as np
import time
from f_toolbox import sls as sls

"""
print("matA = [",end="")
for i in range(0,20):
  print("[",end="")
  for j in range(0,20):
    print(random.randint(0,1000)/100,end="")
    if(j<19):
      print(",",end="")
  print("]",end="")
  if(i<19):
    print(",",end="")
print("]")
"""

mat_size = 100

A = np.random.rand(mat_size, mat_size)
x = np.random.rand(mat_size)



start = time.time()
for i in range(10):
  A_inv = np.linalg.inv(A)
  y = np.dot(A_inv, x)
print(time.time() - start)
print(y[0:10])

start = time.time()
for i in range(10):
  y = np.linalg.solve(A, x)
print(time.time() - start)
print(y[0:10])

start = time.time()
for i in range(10):
  matC = sls.solve(A, x)
print(time.time() - start)
print(matC[0:10])


#print(matC)
