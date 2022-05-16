#!/bin/python3

import numpy
import time
import copy
from scipy.optimize import minimize

from g import g

class ga:

  f = None
  p_in = None
  p_best = None
  rss_best = None
  psize = 0
  popsize = 0
  c_popsize = 0
  pop = None
  pop_rss = None
  f_pop = None
  f_pop_rss = None
  f_popsize = 0
  c_pop = None
  c_pop_rss = None
  threshold = 1.0e20
  parents = None
  populate_counter = 0
  evolve_period = 10
  evolve_amount = 3
  search = None

  @staticmethod
  def run(f, p, niter=1000, popsize=1000, fresh=0.1, threshold=1.0e20, min_final=True, search=None):
    ga.f = f
    ga.p_in = numpy.copy(p)
    ga.psize = len(p)
    ga.popsize = popsize
    ga.c_popsize = popsize
    ga.f_popsize = int(popsize * fresh)
    if(search is None):
      ga.search = 'w'
    else:
      ga.search = search

    ga.setup()
    ga.pop, ga.pop_rss = ga.populate(ga.pop, ga.pop_rss)

    for gn in range(niter):
      ga.generation(gn)
    p = numpy.copy(ga.pop[0, :])
    if(min_final):
      res = minimize(ga.rss, p, method='BFGS', options={'gtol':  1.0e-8, 'maxiter': 10, })
      p = res['x']

    # Clear note
    g.displaynote = ''
    return p
    
  @staticmethod 
  def rss(p):
    rss = ga.f(p)
    if(ga.rss_best is None):
      ga.p_best = numpy.copy(p)
      ga.rss_best = rss
    if(rss < ga.rss_best):
      ga.p_best = numpy.copy(p)
      ga.rss_best = rss
    return rss

  @staticmethod
  def generation(gn):  
    # Shuffle parent list  
    ga.parents = numpy.arange(ga.popsize)
    numpy.random.shuffle(ga.parents)

    # Breed
    
    g.displaynote = 'gen ' + str(gn) + ' breed parents'
    for n in range(ga.popsize // 2):
      pa = numpy.copy(ga.pop[ga.parents[2 * n], :])
      pb = numpy.copy(ga.pop[ga.parents[2 * n + 1], :])
      ca, cb = ga.breed(pa, pb)
      ca = ga.no_clones(ca)
      cb = ga.no_clones(cb)
      ga.c_pop[2 * n, :] = numpy.copy(ca[:])
      ga.c_pop[2 * n+1, :] = numpy.copy(cb[:])
      ga.c_pop_rss[2 * n] = ga.rss(ca)      
      ga.c_pop_rss[2 * n+1] = ga.rss(cb)
            
    # Merge
    ga.merge(ga.pop, ga.pop_rss, ga.c_pop, ga.c_pop_rss)
   
    
    # Shuffle parent list
    ga.parents = numpy.arange(ga.popsize)
    numpy.random.shuffle(ga.parents)
    
    # Fresh
    g.displaynote = 'gen ' + str(gn) + ' populate fresh'
    ga.f_pop, ga.f_pop_rss = ga.populate(ga.f_pop, ga.f_pop_rss)
    
    g.displaynote = 'gen ' + str(gn) + ' breed fresh'
    # Breed
    for n in range(ga.f_popsize):
      loop = True
      while(loop):
        try:
          pa = numpy.copy(ga.pop[ga.parents[n], :])
          pb = numpy.copy(ga.f_pop[n, :])
          ca, cb = ga.breed(pa, pb)
          ca = ga.no_clones(ca)
          cb = ga.no_clones(cb)
          ga.c_pop[2 * n, :] = numpy.copy(ca[:])
          ga.c_pop[2 * n+1, :] = numpy.copy(cb[:])
          ga.c_pop_rss[2 * n] = ga.rss(ca)      
          ga.c_pop_rss[2 * n+1] = ga.rss(cb)
          loop = False
        except:
          pass
      
    # Merge  
    ga.merge(ga.pop, ga.pop_rss, ga.c_pop[:2*ga.f_popsize, :], ga.c_pop_rss[:2 * ga.f_popsize])
    
    if(ga.evolve_period is not None):
      if(gn > 0 and (gn+1)%ga.evolve_period == 0):

        # Shuffle parent list
        ga.parents = numpy.arange(ga.popsize)
        numpy.random.shuffle(ga.parents[2:])
        
        for en in range(ga.evolve_amount):
          g.displaynote = 'gen ' + str(gn) + ' evolve ' + str(en)
          pa = numpy.copy(ga.pop[ga.parents[en], :])
          res = minimize(ga.rss, pa, method='CG', options={'gtol':  1.0e-4, 'maxiter': 1, })
          pe = res['x']
          ga.e_pop[en, :] = pe[:]
          ga.e_pop_rss[en] = ga.rss(pe[:])
          
        # Merge  
        ga.merge(ga.pop, ga.pop_rss, ga.e_pop[:, :], ga.e_pop_rss[:])
          
  @staticmethod
  def breed(pa, pb):
    ca = numpy.zeros((ga.psize,), dtype=numpy.float64,)  
    cb = numpy.zeros((ga.psize,), dtype=numpy.float64,)  
    state = 0  
    for n in range(ga.psize):
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


  @staticmethod
  def no_clones(p):
    for n in range(ga.popsize):
      if(p[:].all() == ga.pop[n,:].all()):
        p[:] = numpy.random.uniform() * p[:]
        return p[:]
    return p[:]

  @staticmethod
  def setup():
    ga.pop = numpy.zeros((ga.popsize, ga.psize,), dtype=numpy.float64,)  
    ga.pop_rss = numpy.zeros((ga.popsize,), dtype=numpy.float64,)  
    ga.c_pop = numpy.zeros((ga.popsize, ga.psize,), dtype=numpy.float64,)  
    ga.c_pop_rss = numpy.zeros((ga.popsize,), dtype=numpy.float64,)  
    ga.f_pop = numpy.zeros((ga.f_popsize, ga.psize,), dtype=numpy.float64,)  
    ga.f_pop_rss = numpy.zeros((ga.f_popsize,), dtype=numpy.float64,)  
    ga.e_pop = numpy.zeros((ga.evolve_amount, ga.psize,), dtype=numpy.float64,)  
    ga.e_pop_rss = numpy.zeros((ga.evolve_amount,), dtype=numpy.float64,)  
  
 
  @staticmethod
  def populate(pop, pop_rss):
    g.displaynote = 'populating (' + str(len(pop)) + ')'
    ga.populate_counter = ga.populate_counter + 1
    ns = 0
    if(ga.populate_counter == 1):
      ns = 1
      pop[0,:] = ga.p_in[:]
      pop_rss[0] = ga.rss(pop[0,:])
    for n in range(ns, len(pop_rss)):
      loop = True
      while(loop):
        r = numpy.random.uniform()
        # n narrow search
        if(ga.search.lower() == 'n'):
          if(r < 0.05):
            p_new = ga.p_in + 10**(-4*numpy.random.uniform()) * ga.p_in * (0.5 - numpy.random.rand(ga.psize))
          elif(r < 0.25):
            p_new = ga.p_best + 0.01 * ga.p_best * (0.5 - numpy.random.rand(ga.psize))
          elif(r < 0.75):
            p_new = ga.p_best + 10**(-1-4*numpy.random.uniform()) * ga.p_best * (0.5 - numpy.random.rand(ga.psize))
          elif(r < 0.90):
            p_new = ga.p_best + 10**(-1-6*numpy.random.uniform()) * ga.p_best * (0.5 - numpy.random.rand(ga.psize))
          # 10% from random
          else:
            p_new = (0.5 - numpy.random.rand(ga.psize)) * 10**(7.0 * (0.5 - numpy.random.rand(ga.psize)))
        # w  wider search
        else:
          if(r < 0.05):
            p_new = ga.p_in + 10**(-4*numpy.random.uniform()) * ga.p_in * (0.5 - numpy.random.rand(ga.psize))
          elif(r < 0.12):
            p_new = ga.p_best + 0.01 * ga.p_best * (0.5 - numpy.random.rand(ga.psize))
          elif(r < 0.24):
            p_new = ga.p_best + 10**(-1-4*numpy.random.uniform()) * ga.p_best * (0.5 - numpy.random.rand(ga.psize))
          elif(r < 0.3):
            p_new = ga.p_best + 10**(-1-6*numpy.random.uniform()) * ga.p_best * (0.5 - numpy.random.rand(ga.psize))
          # 10% small parameters
          elif(r < 0.4):
            p_new = 1.0e-5 * (0.5 - numpy.random.rand(ga.psize))
          # 10% large parameters
          elif(r < 0.5):
            p_new = 1.0e3 * (0.5 - numpy.random.rand(ga.psize))
          # 25% random
          elif(r < 0.75):
            p_new = (0.5 - numpy.random.rand(ga.psize)) * 10**(5.0 * (0.5 - numpy.random.uniform()))
          # 25% random
          elif(r <= 1.0):
            p_new = (0.5 - numpy.random.rand(ga.psize)) * 10**(7.0 * (0.5 - numpy.random.rand(ga.psize)))
        rss = ga.rss(p_new)

        if(rss < ga.threshold or numpy.isnan(rss)):
          loop = False
      pop[n,:] = p_new[:] 
      pop_rss[n] = rss  

    return pop, pop_rss

 
  @staticmethod
  def merge(pop_a, rss_a, pop_b, rss_b):
    ax = pop_a.shape[0]
    ay = pop_a.shape[1]
    bx = pop_b.shape[0]
    by = pop_b.shape[1]
    
    t = numpy.zeros((ax + bx, by + 1,), dtype=numpy.float64,) 
    t[:ax,0] = rss_a[:]
    t[ax:,0] = rss_b[:]
    t[:ax,1:] = pop_a[:,:]
    t[ax:,1:] = pop_b[:,:]
    
    t = t[t[:, 0].argsort()]
    rss_a[:] = t[:ax,0]
    pop_a[:,:] = t[:ax, 1:]
    
    
    rss_b[:] = 0.0
    pop_b[:,:] = 0.0
    
    