#!/bin/python3

import numpy

"""
! Mehl Singh Papaconstantopoulos strains 1993
! Properties of ordered intermetallic alloys: first-principles
! and approximate methods
!
! Mehl Klein Papaconstantopoulos 1993
! First principles calculations of elastic properties of metals
######################################################"""

class msp: 

  @staticmethod
  def calc_c11_c12(s, e0, vol, c11, c12):
    if(e0 > 0.0):
      e0 = -1.0 * e0
    e = e0 + vol * (c11-c12) * s**2
    return e

  @staticmethod
  def calc_c44(s, e0, vol, c44):
    if(e0 > 0.0):
      e0 = -1.0 * e0
    e = e0 + (vol/2) * c44 * s**2
    return e
