#!/bin/python3
"""
relax
"""

import numpy
import os
import time

from scipy.optimize import minimize
from scipy.optimize import basinhopping

from f_toolbox import math

from g import g
from ds import ds
from std import std
from output import output
from label import label
from configs import configs
from nl import nl
from md import md


class relax:

  d = {'config_id': None, 'minimizer': 'bfgs'}

  @staticmethod
  def run():
    output.log("Relax", verbose=2)
    

  @staticmethod
  def set_minimizer(inp):
    if(inp.lower() == 'bfgs'):
      relax.d['minimizer'] = 'bfgs'
    elif(inp.lower() == 'cg'):
      relax.d['minimizer'] = 'cg'


  @staticmethod
  def run_relax(config_id, pmin='e', pvar='a0', options=''):
    ta = time.time()
    output.log("", verbose=2)
    output.log("Relax", verbose=2)


    # Optimise a0 by minimizing energy    
    #####################################################################

    if(pvar == 'a0' and pmin == 'e'):
      output.log("Optimise a0 by minimizing energy", verbose=2)
      ct_original = g.configs[config_id]['calc_type']
      g.configs[config_id]['calc_type'] = 'e'
      relax.d['config_id'] = config_id
      a0 = g.configs[config_id]['a0']
      ca = g.stats['config_calc_count']

      p = relax.minimize(relax.minimise_energy_a0, [a0])
      g.configs[config_id]['a0'] = p[0]

      output.log("Final calc", verbose=2)
      nl.change_cell(config_id, a0=g.configs[config_id]['a0'])
      e = relax.calc_energy(config_id)
      output.log("Optimise a0 by minimizing energy", verbose=2)
      output.log("minimizer: " + relax.d['minimizer'], verbose=2)
      output.log("a0:        " + str(g.configs[config_id]['a0']), verbose=2)
      output.log("c_energy:  " + str(g.configs[config_id]['c_energy']), verbose=2)
      output.log("steps:     " + str(g.stats['config_calc_count'] - ca), verbose=2)
      output.log("time:      " + str(time.time() - ta), verbose=2)
      g.configs[config_id]['calc_type'] = ct_original


    # Optimise a0 by minimizing forces    
    #####################################################################

    if(pvar == 'a0' and pmin == 'f'):
      output.log("Optimise a0 by minimizing forces", verbose=2)
      ct_original = g.configs[config_id]['calc_type']
      g.configs[config_id]['calc_type'] = 'ef'
      relax.d['config_id'] = config_id
      a0 = g.configs[config_id]['a0']
      ca = g.stats['config_calc_count']
      p = relax.minimize(relax.minimise_forces_a0, [a0])
      g.configs[config_id]['a0'] = p[0]
     
      output.log("Final calc", verbose=2)
      nl.change_cell(config_id, a0=g.configs[config_id]['a0'])
      e = relax.calc_energy(config_id)
      output.log("Optimise a0 by minimizing energy", verbose=2)
      output.log("minimizer: " + relax.d['minimizer'], verbose=2)
      output.log("a0:        " + str(g.configs[config_id]['a0']), verbose=2)
      output.log("c_energy:  " + str(g.configs[config_id]['c_energy']), verbose=2)
      output.log("steps:     " + str(g.stats['config_calc_count'] - ca), verbose=2)
      output.log("time:      " + str(time.time() - ta), verbose=2)
      g.configs[config_id]['calc_type'] = ct_original

    # Optimise uv by minimizing energy    
    ##################################################################### 
 
    elif(pvar == 'uv' and pmin == 'e'):
      output.log("Optimise UV by minimizing energy")
      ct_original = g.configs[config_id]['calc_type']
      g.configs[config_id]['calc_type'] = 'e'
      relax.d['config_id'] = config_id
      uv = g.configs[config_id]['uv']
      uvp = relax.uvto9(uv)      

      ca = g.stats['config_calc_count']
      p = relax.minimize(relax.minimise_energy_uv, uvp)
      uv = relax.uvto33(p)
      g.configs[config_id]['uv'] = uv
      
      output.log("Final calc", verbose=2)
      nl.change_cell(config_id, uv=g.configs[config_id]['uv'], scale_uv=False)
      e = relax.calc_energy(config_id)
      output.log("Optimise uv by minimizing energy", verbose=2)
      output.log("minimizer: " + relax.d['minimizer'], verbose=2)
      output.log("uv:        " + str(g.configs[config_id]['uv'][0,0]) + " " + str(g.configs[config_id]['uv'][0,1]) + " " + str(g.configs[config_id]['uv'][0,2]), verbose=2)
      output.log("           " + str(g.configs[config_id]['uv'][1,0]) + " " + str(g.configs[config_id]['uv'][1,1]) + " " + str(g.configs[config_id]['uv'][1,2]), verbose=2)
      output.log("           " + str(g.configs[config_id]['uv'][2,0]) + " " + str(g.configs[config_id]['uv'][2,1]) + " " + str(g.configs[config_id]['uv'][2,2]), verbose=2)
      output.log("c_energy:  " + str(g.configs[config_id]['c_energy']), verbose=2)
      output.log("steps:     " + str(g.stats['config_calc_count'] - ca), verbose=2)
      output.log("time:      " + str(time.time() - ta), verbose=2)
      g.configs[config_id]['calc_type'] = ct_original


    # Optimise uv and a0 by minimizing energy    
    ##################################################################### 

    elif(pvar == 'uva0' and pmin == 'e'):
      output.log("Optimise UV by minimizing energy", verbose=2)
      ct_original = g.configs[config_id]['calc_type']
      g.configs[config_id]['calc_type'] = 'e'
      relax.d['config_id'] = config_id
      uv = g.configs[config_id]['uv']
      uvp = relax.uvto9(uv)
      x = uv[0, 0]

      ca = g.stats['config_calc_count']
      p = relax.minimize(relax.minimise_energy_uv, uvp)
      uv = relax.uvto33(p)
      a0 = g.configs[config_id]['a0']


      y =  uv[0,0]
      if(x != 0.0 and y != 0.0):
        uv = (x / y) * (uv[:,:])
        a0 = a0 * (y / x)

      g.configs[config_id]['uv'] = uv
      g.configs[config_id]['a0'] = a0
 
      output.log("Final calc", verbose=2)
      nl.change_cell(config_id, uv=g.configs[config_id]['uv'], scale_uv=False)
      e = relax.calc_energy(config_id)
      output.log("Optimise uv by minimizing energy", verbose=2)
      output.log("minimizer: " + relax.d['minimizer'], verbose=2)
      output.log("a0:        " + str(g.configs[config_id]['a0']), verbose=2)
      output.log("uv:        " + str(g.configs[config_id]['uv'][0,0]) + " " + str(g.configs[config_id]['uv'][0,1]) + " " + str(g.configs[config_id]['uv'][0,2]), verbose=2)
      output.log("           " + str(g.configs[config_id]['uv'][1,0]) + " " + str(g.configs[config_id]['uv'][1,1]) + " " + str(g.configs[config_id]['uv'][1,2]), verbose=2)
      output.log("           " + str(g.configs[config_id]['uv'][2,0]) + " " + str(g.configs[config_id]['uv'][2,1]) + " " + str(g.configs[config_id]['uv'][2,2]), verbose=2)
      output.log("c_energy:  " + str(g.configs[config_id]['c_energy']), verbose=2)
      output.log("steps:     " + str(g.stats['config_calc_count'] - ca), verbose=2)
      output.log("time:      " + str(time.time() - ta), verbose=2)
      g.configs[config_id]['calc_type'] = ct_original

  
    elif(pvar == 'coords' and pmin == 'e'):
      output.log("Optimise coords by minimizing energy")
      ct_original = g.configs[config_id]['calc_type']
      g.configs[config_id]['calc_type'] = 'e'
      relax.d['config_id'] = config_id
      ca = g.stats['config_calc_count']

      c_original = numpy.copy(g.configs[config_id]['coords'])
      
      md.relax(config_id, rcut=7.0, steps=100, niter=3)

      output.log("Final calc", verbose=2)
      nl.rebuild_nl(config_id)
      e = relax.calc_energy(config_id)

      output.log("Optimise coords by minimizing energy", verbose=2)
      output.log("minimizer: " + relax.d['minimizer'], verbose=2)
      output.log("a0:        " + str(g.configs[config_id]['a0']), verbose=2)
      output.log("uv:        " + str(g.configs[config_id]['uv'][0,0]) + " " + str(g.configs[config_id]['uv'][0,1]) + " " + str(g.configs[config_id]['uv'][0,2]), verbose=2)
      output.log("           " + str(g.configs[config_id]['uv'][1,0]) + " " + str(g.configs[config_id]['uv'][1,1]) + " " + str(g.configs[config_id]['uv'][1,2]), verbose=2)
      output.log("           " + str(g.configs[config_id]['uv'][2,0]) + " " + str(g.configs[config_id]['uv'][2,1]) + " " + str(g.configs[config_id]['uv'][2,2]), verbose=2)
      output.log("c_energy:  " + str(g.configs[config_id]['c_energy']), verbose=2)
      output.log("steps:     " + str(g.stats['config_calc_count'] - ca), verbose=2)
      output.log("time:      " + str(time.time() - ta), verbose=2)
      g.configs[config_id]['calc_type'] = ct_original
     
       
    elif(pvar == 'coords' and pmin == 'ef'):
      # Not Working
      output.log("Optimise coords by minimizing energy", verbose=2)
      ct_original = g.configs[config_id]['calc_type']
      g.configs[config_id]['calc_type'] = 'ef'
      relax.d['config_id'] = config_id
      ca = g.stats['config_calc_count']

      configs.calc(config_id)
      p = relax.c_to_p(g.configs[config_id]['coords'][1:,:])

      res = minimize(relax.minimise_energy_coords, p, method='BFGS', jac=relax.f_jac, 
        options={'gtol':  1.0e-6, 'maxiter': 3, })




      g.configs[config_id]['calc_type'] = ct_original


    output.log("", verbose=2)

  @staticmethod
  def minimize(f, p):
    if(relax.d['minimizer'] == 'cg'):
      res = minimize(f, p, method='CG',)
    elif(relax.d['minimizer'] == 'bfgs'):
      res = minimize(f, p, method='BFGS',
        options={'gtol':  1.0e-4, 'maxiter': 7, })
    elif(relax.d['minimizer'] == 'basinhopping'):
      res = basinhopping(f, p, niter=4, T=1.0, stepsize=0.5)
    return res['x']

  @staticmethod
  def minimise_energy_coords(p):
    p[:] = p[:] % 1.0
    config_id = relax.d['config_id']
    g.configs[config_id]['coords'][:,:] = relax.p_to_c(p)    
    nl.build_nl(config_id, rebuild=True)
    return relax.calc_energy(config_id)
    

  @staticmethod
  def minimise_energy_a0(a0):
    config_id = relax.d['config_id']
    nl.change_cell(config_id, a0=a0[0])
    return relax.calc_energy(config_id)    


  @staticmethod
  def minimise_energy_a0(a0):
    config_id = relax.d['config_id']
    nl.change_cell(config_id, a0=a0[0])
    return relax.calc_energy(config_id)


  @staticmethod
  def minimise_forces_a0(a0):
    config_id = relax.d['config_id']
    nl.change_cell(config_id, a0=a0[0])
    return relax.calc_total_force(config_id)


  @staticmethod
  def minimise_energy_uv(uvp):
    config_id = relax.d['config_id']
    uv = relax.uvto33(uvp)
    nl.change_cell(config_id, uv=uv, scale_uv=False)
    return relax.calc_energy(config_id)


  @staticmethod
  def calc_energy(config_id):
    configs.calc(config_id)
    return g.configs[config_id]['c_energy']    


  @staticmethod
  def calc_total_force(config_id):
    configs.calc(config_id)
    return g.configs[config_id]['c_forces_total']





  @staticmethod
  def f_jac(p):
    config_id = relax.d['config_id']
    p = p % 1.0
    ps = len(p)
    fs = ps // 3
    jac = numpy.zeros((ps, ps, ), dtype=numpy.float64)
    relax.calc_energy(config_id)
    r = numpy.copy(g.configs[config_id]['c_forces'][:,:])
    h = 1.0e-5
    for n in range(ps):
      p_new = numpy.copy(p)
      p_new[n] = (p_new[n] + h) % 1.0
      g.configs[config_id]['coords'][1:,:] = relax.p_to_c(p_new) 
      nl.build_nl(config_id, rebuild=True)  
      relax.calc_energy(config_id)
      jac[:fs, n] = (g.configs[config_id]['c_forces'][1:,0] - r[1:,0]) / h
      jac[fs:2*fs, n] = (g.configs[config_id]['c_forces'][1:,1] - r[1:,1]) / h
      jac[2*fs:,n] = (g.configs[config_id]['c_forces'][1:,2] - r[1:,2]) / h
    return jac




  # Other functions

  @staticmethod
  def uvto33(uvin):
    uv = ds.config_uv(None)
    uv[0, :] = uvin[:3]
    uv[1, :] = uvin[3:6]
    uv[2, :] = uvin[6:]
    return uv

  @staticmethod
  def uvto9(uvin):
    uv = numpy.zeros((9,),dtype=numpy.float64)
    uv[:3] = uvin[0, :]
    uv[3:6] = uvin[1, :]
    uv[6:] = uvin[2, :]
    return uv

  @staticmethod
  def mat3to9(minp):
    mout = numpy.zeros((3,3,),dtype=numpy.float64)
    mout[0,0] = minp[0]
    mout[1,1] = minp[1]
    mout[2,2] = minp[2]
    return mout


  @staticmethod
  def c_to_p(c):
    cs = len(c)
    p = numpy.zeros((3 * cs,), dtype=numpy.float64)
    p[:cs] = c[:,0]
    p[cs:2*cs] = c[:,1]
    p[2*cs:] = c[:,2]
    return p


  @staticmethod
  def p_to_c(p):
    cs = len(p) // 3
    c = numpy.zeros((cs,3,), dtype=numpy.float64)
    c[:,0] = p[:cs] 
    c[:,1] = p[cs:2*cs]
    c[:,2] = p[2*cs:]
    return c


  """
  @staticmethod
  def minimise_energy_uv_jac(uvp):
    config_id = relax.d['config_id']
    jac = numpy.zeros((9,),dtype=numpy.float64)

    for n in range(9):
      uvp_m = numpy.copy(uvp)
      uvp_m[n] = uvp_m[0] + 1.1
      uv = relax.uvto33(uvp_m)
      print(uv)
      nl.change_cell(config_id, uv=uv, scale_uv=False)
      jac[n] = relax.calc_energy(config_id)
      print(jac[n])
    print(jac)
  """