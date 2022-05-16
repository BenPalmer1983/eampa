#!/bin/python3

import numpy
import time
import os
import matplotlib.pyplot as plt

from g import g
from ds import ds


class search_step:

  @staticmethod
  def run(f, p):
    rss_initial = f(p)
    pvar_min = numpy.zeros((len(p),), dtype=numpy.float64)
    pvar_max = numpy.zeros((len(p),), dtype=numpy.float64)

    for i in range(-7, 7):
      p_new = numpy.copy(p)
      pvar = numpy.zeros((len(p),), dtype=numpy.float64)
      pvar[:] = 10**i
      p_new = p_new + pvar
      rss = f(p_new)
      if(rss < 10 * rss_initial and sum(pvar) > sum(pvar_max)):
        pvar_max = numpy.copy(pvar)

    for i in range(-7, 7):
      p_new = numpy.copy(p)
      pvar = numpy.zeros((len(p),), dtype=numpy.float64)
      pvar[:] = 10**i
      p_new = p_new + pvar
      rss = f(p_new)
      if(rss < 10 * rss_initial and sum(pvar) > sum(pvar_min)):
        pvar_min = numpy.copy(pvar)


    pvar = pvar_max
    if(sum(pvar_min) > sum(pvar_max)):
      pvar = pvar_min

    return pvar