#!/bin/python3
"""
md
"""

import numpy
import os
import time

from f_toolbox import math
from f_toolbox import atom

from g import g
from ds import ds
from std import std
from output import output
from label import label
from configs import configs
from nl import nl


class md:

  log = []

  @staticmethod
  def run(config_id):
    # Log
    output.log("MD", verbose=0)

    # Make dirs
    md.make_dirs()



  @staticmethod
  def go(config_id):

    nl.build(config_id)
    configs.calc(config_id)
    print(g.configs[config_id]['c_energy'])

    clabels = g.configs[config_id]['labels']
    cin = g.configs[config_id]['coords']
    a0 = g.configs[config_id]['a0']
    uv = g.configs[config_id]['uv']
    nl_id = g.configs[config_id]['nl_id']

    nlabels = g.nl[nl_id]['labels']
    r = g.nl[nl_id]['r']
    rvec = g.nl[nl_id]['rvec']
    inhalo = g.nl[nl_id]['inhalo']
    vinit = numpy.zeros((len(clabels), 3), dtype=numpy.float64, order='F')
    vout = numpy.zeros((len(clabels), 3), dtype=numpy.float64, order='F')
    cout = numpy.zeros((len(clabels), 3), dtype=numpy.float64, order='F')
    atom.verlet_step(clabels, cin, a0, uv, 7.0, vinit, 0.95, 500, 10, 0.01, True, cout, vout)
    md.make_xyz(atom.chistory, clabels, a0, uv)
    md.make_e(atom.ehistory)


    g.configs[config_id]['coords'] = cout
    nl.rebuild_nl(config_id)
    configs.calc(config_id)
    output.log("End", verbose=0)



  @staticmethod
  def relax(config_id, rcut=7.0, steps=100, niter=1):
    # Relax coords
    clabels = g.configs[config_id]['labels']
    a0 = g.configs[config_id]['a0']
    uv = g.configs[config_id]['uv']
    nl_id = g.configs[config_id]['nl_id']

    nlabels = g.nl[nl_id]['labels']
    r = g.nl[nl_id]['r']
    rvec = g.nl[nl_id]['rvec']
    inhalo = g.nl[nl_id]['inhalo']
    vinit = numpy.zeros((len(clabels), 3), dtype=numpy.float64, order='F')
    vout = numpy.zeros((len(clabels), 3), dtype=numpy.float64, order='F')
    cout = numpy.zeros((len(clabels), 3), dtype=numpy.float64, order='F')

    for n in range(niter):
      cin = g.configs[config_id]['coords']
      atom.relax(clabels, cin, a0, uv, rcut, steps, 10, 0.01, cout)
      g.configs[config_id]['coords'] = cout
      nl.rebuild_nl(config_id)



  # Make additional directories for efs
  @staticmethod
  def make_dirs():   
    g.dir['md'] = os.path.join(g.dir['out'], "md")
    for k in g.dir.keys():
      os.makedirs(g.dir[k], exist_ok=True)


  @staticmethod
  def make_xyz(cin, clabels, a0, uv):
    fh = open(os.path.join(g.dir['md'], 'history.xyz'), 'w')
    for n in range(len(cin)):
      fh.write(str(len(cin[n, :, 0])) + '\n')
      fh.write('\n')
      for i in range(len(cin[n, :, 0])):
        c = a0 * numpy.matmul(uv, cin[n, i, :])
        fh.write(label.get(clabels[i]) + '  ' + str(c[0]) + '  ' + str(c[1]) + '  ' + str(c[2]) + '\n')
    fh.close()


  @staticmethod
  def make_e(e):
    fh = open(os.path.join(g.dir['md'], 'energy.dat'), 'w')
    for n in range(len(e)):
      n_str = str(n)
      while(len(n_str)<8):
        n_str = n_str + " "
      fh.write(n_str + str(e[n]) + '\n')
    fh.close()
  
