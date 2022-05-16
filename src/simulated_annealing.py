#!/bin/python3

import numpy
import time
import copy

from g import g
from ds import ds

class simulated_annealing:


  @staticmethod
  def run(f, p, niter=1000, titer=10, tstart=100, tend=0.01, pfact=0.7, pvar=None, vartype=None, gaussian=False):

    if(pvar == None or vartype == None):
      g.displaynote = 'Initializing pvar'
      p, pvar = simulated_annealing.pvar(f, p)
    elif(vartype == 'add'):
      pvar = pvar + numpy.zeros((len(p),), dtype=numpy.float64)
    elif(vartype == 'mult'):
      pvar = pvar * p
    p = simulated_annealing.optimize(f, p, pvar, niter, titer, tstart, tend, pfact, gaussian) 
    g.displaynote = ''
    return p
    
  @staticmethod
  def optimize(f, p, pvar, niter=1000, titer=10, tstart=100, tend=0.01, pfact=0.7, gaussian=False):
    
    temp = numpy.linspace(tstart, tend, titer)

    rss = f(p)
    rss_best = rss
    p_best = numpy.copy(p)

    for tn in range(titer):
      t = temp[tn]
      dp = pfact**(tn) * 10.0 * pvar[:]

      for n in range(niter):
        g.displaynote = 'T Loop ' + str(tn) + "    N Loop " + str(n) + "  t="+str("{0:7.3f}".format(t)) + "  pfact="+str("{0:7.3f}".format(pfact**(tn)))
        p_new = numpy.copy(p)

        if(gaussian):
          p_new[:] = p_new[:] + dp * numpy.random.normal(0.0, 0.1, len(p_new[:]))
        else:
          p_new[:] = p_new[:] + dp * (0.5 - numpy.random.rand(len(p_new[:]))) 


        rss_new = f(p_new)
        if(rss_new < rss or numpy.random.uniform() < numpy.exp((rss_new - rss) / t)):
          rss = rss_new
          p = numpy.copy(p_new)
        if(rss < rss_best):
          rss_best = rss_new
          p_best = numpy.copy(p_new)
      p = numpy.copy(p_best)
      rss = rss_best
    return p_best


  @staticmethod
  def pvar(f, p):
    p_upper = numpy.zeros((len(p),), dtype=numpy.float64)
    p_lower = numpy.zeros((len(p),), dtype=numpy.float64)

    rss_initial = f(p)
    var = numpy.zeros((14,), dtype=numpy.float64)
    for n in range(14):
      var[n] = 10**(-7 + n)


    for pn in range(len(p)):
      for n in range(14):
        p_new = numpy.copy(p)
        p_new[pn] = p_new[pn] + var[n]
        rss = f(p_new)
        if(rss > 10.0 * rss_initial):
          p_upper[pn] = p_new[pn]
          break
    for pn in range(len(p)):
      for n in range(14):
        p_new = numpy.copy(p)
        p_new[pn] = p_new[pn] - var[n]
        rss = f(p_new)
        if(rss > 10.0 * rss_initial):
          p_lower[pn] = p_new[pn]
          break
    
 
    pvar = p_upper - p_lower
    p = p_lower + 0.5 * pvar
    return p, pvar






  """




  @staticmethod
  def opt(p, p_var, frss, loops=1000, pfact=0.99, t_start=100, t_end=1.0, t_steps=10):

    simanneal.pnew = p
    simanneal.pbest = simanneal.pnew
    simanneal.rssbest = frss(p)
    simanneal.rssnew = simanneal.rssbest
    
    t = t_start
    tn = 0
    while(t >= t_end):      
      for n in range(loops):
        p = simanneal.pnew + pfact**(tn) * p_var * (0.5 - numpy.random.rand(len(simanneal.pnew))) 
        rss = frss(p)

        # If better, accept
        if(rss < simanneal.rssnew):
          simanneal.pnew = p
          simanneal.rssnew = rss
        else:
          # If worse, accept only with a small probability
          if(numpy.random.uniform() < numpy.exp((simanneal.rssnew - rss) / t)):
            simanneal.pnew = p
            simanneal.rssnew = rss
        if(rss < simanneal.rssbest):
          simanneal.pbest = p
          simanneal.rssbest = rss
 
      # Update new with best
      simanneal.pnew = simanneal.pbest
      simanneal.rssnew = simanneal.rssbest

      # Cool
      t = t - (t_start - t_end) / (t_steps-1)
      tn = tn + 1

    return simanneal.pbest
  """