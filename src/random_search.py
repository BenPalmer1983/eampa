#!/bin/python3

import numpy
import time
import copy

from g import g
from ds import ds

class random_search:


  @staticmethod
  def run(f, p, niter=1000):

    p[:] = 0.0
    rss = f(p)    
    rss_best = rss
    p_best = numpy.copy(p) 


    for n in range(niter):
      p = random_search.rand(p, n, niter)
      rss = f(p)    
      if(not numpy.isnan(rss) and rss < rss_best):
        rss_best = rss
        p_best = numpy.copy(p) 

    return p_best


  @staticmethod
  def rand(p, n, niter):
    r = numpy.random.uniform()
    if(r < 0.2):
      return (0.5 - numpy.random.rand(len(p)))
    elif(r < 0.3):
      return 0.01 * (0.5 - numpy.random.rand(len(p)))
    elif(r < 0.4):
      return 10.0 * (0.5 - numpy.random.rand(len(p)))
    elif(r < 0.6):
      return 2.0 * (numpy.random.rand(len(p)) - 0.5) * 10**(((niter - n) /niter) * (-7 + 15 * numpy.random.rand(len(p))))
    elif(r < 0.9):
      return 2.0 * (numpy.random.rand(len(p)) - 0.5) * 10**((-4 + 6 * numpy.random.rand(len(p))))
    else:
      return 2.0 * (numpy.random.rand(len(p)) - 0.5) * 10**((-7 + 15 * numpy.random.rand(len(p))))

