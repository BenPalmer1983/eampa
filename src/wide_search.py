#!/bin/python3

import numpy
import time
import copy

from scipy.optimize import minimize

from g import g
from ds import ds

class wide_search:


  @staticmethod
  def run(f, p):
    threshold = 1.0e10

    p_upper = numpy.zeros((len(p),), dtype=numpy.float64)
    p_lower = numpy.zeros((len(p),), dtype=numpy.float64)

    
    var = numpy.zeros((14,), dtype=numpy.float64)
    for n in range(14):
      var[n] = 10**(-7 + n)


    for pn in range(len(p)):
      for n in range(14):
        p_new = numpy.copy(p)
        p_new[pn] = p_new[pn] + var[n]
        rss = f(p_new)
        if(rss > threshold):
          p_upper[pn] = p_new[pn]
          break


    for pn in range(len(p)):
      for n in range(14):
        p_new = numpy.zeros((len(p),), dtype=numpy.float64)
        p_new[pn] = p_new[pn] - var[n]
        rss = f(p_new)
        if(rss > threshold):
          p_lower[pn] = p_new[pn]
          break

    



