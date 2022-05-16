#!/bin/python3

import numpy
import time
import copy

from scipy.optimize import minimize

from g import g
from ds import ds

class hybrid_search:

  samples = None
  samples_rss = None
  minima = None
  minima_rss = None
  scartch = None
  scartch_rss = None
  test_counter = 0
  pool = None

  @staticmethod
  def run(f, p, threshold=1000.0, sample_size=100, minima_size=10):

    # Create arrays
    hybrid_search.samples = numpy.zeros((sample_size, len(p),), dtype=numpy.float64)
    hybrid_search.samples_rss = numpy.zeros((sample_size), dtype=numpy.float64)
    hybrid_search.minima = numpy.zeros((minima_size, len(p),), dtype=numpy.float64)
    hybrid_search.minima_rss = numpy.zeros((minima_size), dtype=numpy.float64)
    hybrid_search.pool = numpy.zeros((1, len(p),), dtype=numpy.float64)
    hybrid_search.pool[0,:] = p

    rss = f(p)
    threshold = rss * threshold

    for i in range(3):
      if(i == 0):
        search = 'w'
      else:
        search = 'n'
      p = hybrid_search.inner(f, p, threshold, search)
      threshold = threshold * 0.9
    return p 

  @staticmethod
  def inner(f, p, threshold, search):
    """###########################################
    Search over a wide parameter range to sample
    points within set threshold for minimization
    ###########################################"""
    rss = f(p)
    hybrid_search.samples[0,:] = numpy.copy(p)
    hybrid_search.samples_rss[0] = rss
    
    sn = 1
    while(sn<len(hybrid_search.samples)):
      p_new = hybrid_search.rand_p(search)      
      rss = f(p_new)
      if(rss < threshold):
        hybrid_search.samples[sn,:] = numpy.copy(p_new)
        hybrid_search.samples_rss[sn] = rss
        sn = sn + 1

    
    """###########################################
    Randomly choose points to minimise
    ###########################################"""
    sn = 0
    mn = 0
    while(mn<len(hybrid_search.minima)):
      if(numpy.random.uniform() < 0.05):
        res = minimize(f, hybrid_search.samples[sn, :], method='BFGS', options={'gtol':  1.0e-8, 'maxiter': 3, })
        hybrid_search.minima[mn, :] = res['x'][:]
        hybrid_search.minima_rss[mn] = res['fun']
        mn = mn + 1
      sn = (sn + 1) % len(hybrid_search.samples)

    
    """###########################################
    Merge samples so far
    ###########################################"""

    temp, temp_rss = hybrid_search.merge_and_sort(hybrid_search.minima, hybrid_search.minima_rss, hybrid_search.samples, hybrid_search.samples_rss)

    
    """###########################################
    Breed best
    ###########################################"""
    
    minima_size = len(hybrid_search.minima_rss)
    pop, pop_rss = hybrid_search.genetic(f, temp[:2*minima_size,:], temp_rss[:2*minima_size], gens=5)
    hybrid_search.pool = numpy.copy(pop)


    """###########################################
    Breed best
    ###########################################"""
    res = minimize(f, pop[0,:], method='BFGS', options={'gtol':  1.0e-8, 'maxiter': 6, })
    p = numpy.copy(res['x'])

    return p


  @staticmethod
  def sort_rss(a, a_rss):
    ah = len(a)
    aw = len(a[0])
    t = numpy.zeros((ah, aw + 1,), dtype=numpy.float64)
    t[:,0] = numpy.copy(a_rss[:])
    t[:,1:] = numpy.copy(a[:,:])
    t = t[t[:, 0].argsort()]
    return t[:,1:], t[:,0]


  @staticmethod
  def merge_and_sort(a, a_rss, b, b_rss):
    ax = a.shape[0]
    ay = a.shape[1]
    bx = b.shape[0]
    by = b.shape[1]
    
    t = numpy.zeros((ax + bx, by + 1,), dtype=numpy.float64,) 
    t[:ax,0] = a_rss[:]
    t[ax:,0] = b_rss[:]
    t[:ax,1:] = a[:,:]
    t[ax:,1:] = b[:,:]
    
    t = t[t[:, 0].argsort()]

    out_rss = numpy.copy(t[:,0])
    out = numpy.copy(t[:,1:])
    return out, out_rss





  @staticmethod
  def genetic(f, pop, pop_rss, gens=5, inner=5):
    pl = len(pop)
    pw = len(pop[0])
    temp = numpy.zeros(((inner + 1) * pl, pw), dtype=numpy.float64,)
    temp_rss = numpy.zeros(((inner + 1) * pl), dtype=numpy.float64,)
 
    temp[:pl, :] = numpy.copy(pop[:, :])
    temp_rss[:pl] = numpy.copy(pop_rss[:])

    parents = numpy.arange(pl)
    for gn in range(gens):
      n = pl
      for nn in range(inner):
        numpy.random.shuffle(parents)
        for pn in range(pl // 2):
          ca, cb = hybrid_search.breed(temp[parents[2*pn],:], temp[parents[2*pn+1],:])
          ca = hybrid_search.no_clones(ca, temp[:n, :])
          temp[n, :] = numpy.copy(ca)
          temp_rss[n] = f(ca)
          n = n + 1
          cb = hybrid_search.no_clones(cb, temp[:n, :])
          temp[n, :] = numpy.copy(cb)
          temp_rss[n] = f(cb)
          n = n + 1
      temp, temp_rss = hybrid_search.sort_rss(temp, temp_rss)
    pop = temp[:pl, :]
    pop_rss = temp_rss[:pl]
    return pop, pop_rss


    
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
  def rand_p(search): 
    p = hybrid_search.pool[numpy.random.randint(0, len(hybrid_search.pool))]

    if(search == 'w'):
      p_seed = numpy.zeros((len(p),),dtype=numpy.float64)
      for n in range(len(p)):
        if(p[n] == 0.0):
          p_seed[n] = 0.1
        else:
          p_seed[n] = p[n]
      r = numpy.random.uniform()
      # 20% chance to randomly scale and change sign
      if(r < 0.20):
        p_seed[:] = 2.0 * (0.5 - numpy.random.rand(len(p[:]))) * p_seed[:]
      elif(r < 0.30):
        p_seed[:] = 0.1 * (0.5 - numpy.random.rand(len(p[:])))
      elif(r < 0.60):
        p_seed[:] = 10.0 * (0.5 - numpy.random.rand(len(p[:])))
      else:
        p_seed[:] = (1.0 + 10**(-2 + 3 * numpy.random.uniform()) * (0.5 - numpy.random.rand(len(p[:])))) * p_seed[:]
      return p_seed

    if(search == 'n'):
      p_seed = numpy.copy(p)
      r = numpy.random.uniform()
      # 5% chance to randomly change sign of parameters
      if(r < 0.05):
        for n in range(len(p_seed)):
          if(numpy.random.uniform() < 0.5):
            p_seed[n] = -1.0 * p_seed[n]
      r = numpy.random.uniform()
      if(r < 0.30):
        p_seed[:] = p_seed[:] + 0.0001 * (0.5 - numpy.random.rand(len(p[:])))
      elif(r < 0.60):
        p_seed[:] = p_seed[:] + 0.01 * (0.5 - numpy.random.rand(len(p[:])))
      elif(r < 0.80):
        p_seed[:] = p_seed[:] + (0.5 - numpy.random.rand(len(p[:]))) * p_seed[:]
      else:
        p_seed[:] = p_seed[:] + 0.01 * (0.5 - numpy.random.rand(len(p[:]))) * p_seed[:]
      return p_seed
      
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
    p = hybrid_search.run(hybrid_search.test_objective, p, threshold=100.0, sample_size=100)
    #print(hybrid_search.test_rss(hybrid_search.test_objective, p))
    print(p)
    print(hybrid_search.test_counter)


if __name__ == "__main__":
  hybrid_search.test()    
