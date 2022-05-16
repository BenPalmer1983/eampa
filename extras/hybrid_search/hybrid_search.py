#!/bin/python3

import numpy
from scipy.optimize import minimize


class hybrid_search:


  test_counter = 0
  pool = None
  samples = None


  @staticmethod
  def run(f, p, sample_size=100, minima_size=10):
    """###########################################
    Make a pool of parameters to use for new random
    ###########################################"""
    hybrid_search.make_pool(f, p)

    """###########################################
    Search over a wide parameter range to sample
    points within set threshold for minimization
    ###########################################"""

    threshold = 10.0 * hybrid_search.pool[-1,0]
    hybrid_search.samples = numpy.zeros((sample_size, len(p)+1,), dtype=numpy.float64)

    sn = 0
    while(sn<sample_size):
      p_new = hybrid_search.rand_p(p)
      rss = f(p_new)
      if(rss < threshold):
        hybrid_search.samples[sn,0] = rss
        hybrid_search.samples[sn,1:] = p_new
        sn = sn + 1

    # Order
    hybrid_search.samples = hybrid_search.samples[hybrid_search.samples[:, 0].argsort()]
    
    keys = numpy.arange(len(hybrid_search.samples))
    numpy.random.shuffle(keys) 

    rss_best  = hybrid_search.samples[0,0]  
    mn = 0
    kn = 0
    while(mn < minima_size):
      if(numpy.random.uniform() < (rss_best / hybrid_search.samples[keys[kn],0])):
        res = minimize(f, hybrid_search.samples[keys[kn], 1:], method='BFGS', options={'gtol':  1.0e-8, 'maxiter': 3, })
        hybrid_search.samples[keys[kn], 0] = res['fun']
        hybrid_search.samples[keys[kn], 1:] = res['x'][:]
        mn = mn + 1
      kn = (kn + 1) % len(keys)

    hybrid_search.samples = hybrid_search.samples[hybrid_search.samples[:, 0].argsort()]
    hybrid_search.samples = hybrid_search.genetic(f, hybrid_search.samples, gens=2, inner=1)

    mn = 0
    while(mn < minima_size):
      res = minimize(f, hybrid_search.samples[mn, 1:], method='BFGS', options={'gtol':  1.0e-8, 'maxiter': 3, })
      hybrid_search.samples[mn, 0] = res['fun']
      hybrid_search.samples[mn, 1:] = res['x'][:]
      mn = mn + 1
    
    hybrid_search.samples = hybrid_search.samples[hybrid_search.samples[:, 0].argsort()]
    res = minimize(f, hybrid_search.samples[0, 1:], method='BFGS', options={'gtol':  1.0e-8, 'maxiter': 3, })
    
    return res['x']



  """ MAKE A POOL OF PARAMETERS """

  @staticmethod
  def make_pool(f, p):
    pool_size = 20
    scratch_size = 100
    threshold = 1.0e10
    wide_narrow = 0.4

    scratch = numpy.zeros((scratch_size, len(p)+1,), dtype=numpy.float64) 

    # Input parameters
    sn = 0
    scratch[sn, 0] = f(p[:])
    scratch[sn, 1:] = numpy.copy(p[:])
    hybrid_search.rand_wide(p)
    sn = sn + 1

    while(sn < scratch_size):
      if(numpy.random.uniform() < wide_narrow):
        p_new = hybrid_search.rand_wide(p)
        rss = f(p_new[:])
        if(rss<threshold):
          scratch[sn, 0] = rss
          scratch[sn, 1:] = numpy.copy(p_new[:])
          sn = sn + 1
      else:
        p_new = hybrid_search.rand_narrow(p)
        rss = f(p_new[:])
        if(rss<threshold):
          scratch[sn, 0] = rss
          scratch[sn, 1:] = numpy.copy(p_new[:])
          sn = sn + 1
    
    # Sort scratch
    scratch = scratch[scratch[:, 0].argsort()]

    # Breed results
    hybrid_search.genetic(f, scratch, gens=1, inner=2)
  
    # Save pool
    hybrid_search.pool = numpy.copy(scratch[:pool_size,:])

 
    

  """ GENETIC """

  @staticmethod
  def genetic(f, pool, gens=5, inner=5):
    pl = len(pool)
    pw = len(pool[0])
    temp = numpy.zeros(((inner + 1) * pl, pw), dtype=numpy.float64,)
 
    temp[:pl, :] = numpy.copy(pool[:, :])

    parents = numpy.arange(pl)
    for gn in range(gens):
      n = pl
      for nn in range(inner):
        numpy.random.shuffle(parents)
        for pn in range(pl // 2):
          ca, cb = hybrid_search.breed(temp[parents[2*pn],1:], temp[parents[2*pn+1],1:])
          ca = hybrid_search.no_clones(ca, temp[:n, 1:])
          temp[n, 1:] = numpy.copy(ca)
          temp[n, 0] = f(ca)
          n = n + 1
          cb = hybrid_search.no_clones(cb, temp[:n, 1:])
          temp[n, 1:] = numpy.copy(cb)
          temp[n, 0] = f(cb)
          n = n + 1
      temp = temp[temp[:, 0].argsort()]
    pool = temp[:pl, :]
    return pool
    
  @staticmethod
  def no_clones(p, pop):
    for n in range(len(pop)):
      if(p[:].all() == pop[n,:].all()):
        p[:] = numpy.random.uniform() * p[:]
        return p[:]
    return p[:]

  @staticmethod
  def breed(pa, pb):
    ca = numpy.zeros((len(pa),), dtype=numpy.float64,)  
    cb = numpy.zeros((len(pa),), dtype=numpy.float64,)  
    state = 0  
    for n in range(len(pa)):
      if(state == 0):
        ca[n] = pa[n]
        cb[n] = pb[n]
      else:
        cb[n] = pa[n]
        ca[n] = pb[n]
      if(numpy.random.uniform() > 0.5):
        if(state == 0): 
          state = 1
        else:
          state = 0
    return ca, cb






  """ RANDOM PARAMETER FUNCTIONS """


  @staticmethod
  def rand_wide(p): 
    return numpy.random.normal(size=len(p[:])) * 10**(5*numpy.random.normal(size=len(p[:])))
  

  @staticmethod
  def rand_narrow(p, m=0.01):
    r = numpy.random.uniform() 
    if(r < 0.50):
      return p[:] * (1.0 + m * numpy.random.normal(size=len(p[:])))
    else:
      return p[:] + m * numpy.random.normal(size=len(p[:]))


  @staticmethod
  def rand_p(p, m=0.1): 
    p_new = hybrid_search.p_from_pool(p)
    r = numpy.random.uniform() 
    if(r < 0.50):
      return p_new * (1.0 + m * numpy.random.normal(size=len(p_new)))
    else:
      return p_new + m * numpy.random.normal(size=len(p_new))


  @staticmethod
  def p_from_pool(p): 
    rss_best = hybrid_search.pool[0,0]
    keys = numpy.arange(len(hybrid_search.pool))
    numpy.random.shuffle(keys)    
    for kn in range(len(keys)):
      if(numpy.random.uniform() < (rss_best / hybrid_search.pool[keys[kn],0])):
        return hybrid_search.pool[keys[kn],1:]










      
  """ TESTING """

  @staticmethod
  def test_objective(p):
     x = p[0]
     y = p[1]
     hybrid_search.test_counter = hybrid_search.test_counter + 1
     return (x**2 + y - 11)**2 + (x + y**2 -7)**2

  @staticmethod
  def test():

    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import Axes3D
 
    from numpy import arange
    from numpy import meshgrid
 
 
    # define range for input
    r_min, r_max = -5.0, 5.0
    # sample input range uniformly at 0.1 increments
    xaxis = arange(r_min, r_max, 0.1)
    yaxis = arange(r_min, r_max, 0.1)
    # create a mesh from the axis
    x, y = meshgrid(xaxis, yaxis)
    # compute targets
    p = [x, y]
    results = hybrid_search.test_objective(p)
    # create a surface plot with the jet color scheme
    figure = pyplot.figure()
    axis = figure.gca(projection='3d')
    axis.plot_surface(x, y, results, cmap='jet')
    # show the plot
    #pyplot.show()

    p = [-3.5, 7.8]
    #for n in range(20):
    #  print(hybrid_search.rand_p(p))
    hybrid_search.test_counter = 0
    p = hybrid_search.run(hybrid_search.test_objective, p, sample_size=100)
    #print(hybrid_search.test_rss(hybrid_search.test_objective, p))
    print(p)
    print(hybrid_search.test_counter)


if __name__ == "__main__":
  hybrid_search.test()    
