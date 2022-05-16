import random
import numpy
import time
from f_toolbox import interp as interp
import matplotlib.pyplot as plt


def f(x):
  return 4 + 1.5 * x - 0.5 * x**2 + 0.04 * x**3

def df(x):
  return 1.5 - x + 0.12 * x**2




x = numpy.linspace(0.0, 10.0, 20)
y = f(x)
dy = df(x)

yi = interp.interpolate(x[10], x[:], y[:], 4, 0) 
print(yi, y[10])
yi = interp.interpolate(x[10]+0.001, x[:], y[:], 4, 0) 
print(y[10], yi, y[10])

expanded = numpy.zeros((1001, 4), dtype=numpy.float64, order='F')
interp.fill(x, y, 5, expanded)



plt.plot(x, y, ls='None', marker='x')
plt.plot(x, dy, ls='None', marker='+')
plt.plot(expanded[:,0], expanded[:,1], ls='solid')
plt.plot(expanded[:,0], expanded[:,2], ls='dashed')
plt.show()






