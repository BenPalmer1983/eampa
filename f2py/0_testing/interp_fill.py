import random
import numpy
import time
from f_toolbox import interp as interp
import matplotlib.pyplot as plt



d = numpy.genfromtxt('data.csv', delimiter=',')
dnew = numpy.zeros((len(d), 3), dtype=numpy.float64, order='F')

interp.fill(d[:,0],d[:,1],5, dnew)

print(d)
print(dnew)



plt.plot(d[:,0], d[:,1], ls='dashed', marker="+")
plt.plot(d[:,0], d[:,2], ls='dashed', marker="+")
plt.plot(dnew[:,0], dnew[:,1], ls='solid')
plt.plot(dnew[:,0], dnew[:,2], ls='solid')
plt.show()
