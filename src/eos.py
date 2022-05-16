#!/bin/python3

import numpy
import time

from scipy.optimize import minimize

from f_toolbox import polyfit

from g import g
from ds import ds
from units import units



class eos: 

  v = None
  e = None

  @staticmethod
  def fit(v, e):
    eos.v = v
    eos.e = e
    eos_p = {"e0": 0.0, "v0": 0.0, "b0": 0.0, "b0p": 0.0, "b0_gpa": 0.0,}

    poly = polyfit.fit(v, e, 2)
    if(poly[2] == 0):
      poly[2] = 1.0e-20
    eos_p['v0'] = (-1 * poly[1]) / (2 * poly[2])
    eos_p['e0'] = (poly[2] * eos_p['v0'] * eos_p['v0']) + (poly[1] * eos_p['v0']) + poly[0]
    eos_p['b0'] = 2.0 * poly[2] * eos_p['v0']
    eos_p['b0p'] = 2.0

    eos_p['b0_gpa'] = units.convert('EV/ANG3', 'GPA', eos_p['b0'])
    return eos_p

  
  @staticmethod
  def bm_rss(p):
    return sum((eos.e[:] - eos.bm_calc(p, eos.v[:]))**2)


  @staticmethod
  def bm_calc(p, V):
    V0 = p[0]
    E0 = p[1]
    B0 = p[2]
    B0P = p[3]
    if(V0 == 0.0):
      V0 = 1.0e-20
    try:
      eta = (V/V0)**(1.0/3.0)
    except:
      eta = 1.0e10
    return E0 + (9.0/16.0) * (B0 * V0) * ((eta*eta - 1)*(eta*eta - 1)) * (6.0 + B0P * (eta * eta - 1) - 4.0 * eta * eta ) 



