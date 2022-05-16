#!/bin/python3

import numpy

"""
Rose equation of state (see G. Bonny code manual)


######################################################"""

class rose: 

  @staticmethod
  def calc(a, a0, b0, ecoh, vol):
    if(ecoh < 0.0):
      ecoh = -1.0 * ecoh

    alpha = numpy.sqrt(9 * vol * b0 / ecoh)
    abar = alpha * (a / a0 - 1)
    e = -1.0 * ecoh * (1+abar) * numpy.exp(-abar)
    return e

