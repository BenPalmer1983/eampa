import numpy



for n in range(10):
  x = numpy.random.normal()
  y = numpy.random.normal()
  z = numpy.random.normal()
  r = 0.025 * numpy.random.uniform()
  x, y, z = r*x, r*y, r*z

  print(x, y, z)
