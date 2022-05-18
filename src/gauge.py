#!/bin/python3

import numpy
import time
import os

from f_toolbox import atom
from f_toolbox import math
from f_toolbox import transforms
from f_toolbox import interp

from g import g
from ds import ds
from configs import configs
from potential import potential
from nl import nl
from output import output
from plot import plot
from relax import relax
from eos import eos
from units import units
from rose import rose
from msp import msp

class gauge: 

  log = []


  @staticmethod
  def run():
    output.log("Run Effective Gauge Transformation", verbose=0)

    gauge.make_dirs()

    gauge.process()


    exit() 


  # Make additional directories for bp
  @staticmethod
  def make_dirs():    

    g.dir['gauge'] = os.path.join(g.dir['out'], "gauge")
    g.dir['gauge_original'] = os.path.join(g.dir['gauge'], "original")
    g.dir['gauge_transformed'] = os.path.join(g.dir['gauge'], "transformed")
    g.dir['gauge_original_plot'] = os.path.join(g.dir['gauge'], "original_plot")
    g.dir['gauge_transformed_plot'] = os.path.join(g.dir['gauge'], "transformed_plot")

    for k in g.dir.keys():
      os.makedirs(g.dir[k], exist_ok=True)



  @staticmethod
  def process(): 

    potential.tab_for_output(10.0)
    output.potentials_tabulated(g.dir['gauge_original'])
    plot.potentials(g.dir['gauge_original_plot'], source='tab_for_output')

    output.log("Test structures.", verbose=0)
    structures = g.gauge['structures']
    label = g.gauge['label']
    a0 = g.gauge['a0']
    uv = g.gauge['uv']
    size = g.gauge['size']
    rcut = g.gauge['rcut']

    relaxing = []
    for structure in structures:
      config_id = configs.add_common(structure, label, a0, uv, size, rcut, 'gauge', 'e')
      relaxing.append({
                       'config_id': config_id,
                       'structure': structure,
                       'a0': None,
                       'e0': None,
                      })

    m_rn = None
    m_val = 0.0
    for rn in range(len(relaxing)):
      config_id = relaxing[rn]['config_id']
      relax.run_relax(config_id, pmin='e', pvar='a0')
      relaxing[rn]['a0'] = g.configs[config_id]['a0']
      relaxing[rn]['e0'] = g.configs[config_id]['c_energy'] / g.configs[config_id]['count']
      if(m_rn is None):
        m_rn = rn
        m_val = relaxing[rn]['e0']
      elif(relaxing[rn]['e0'] < m_val):
        m_rn = rn

    for rn in range(len(relaxing)):
      flag = ""
      if(rn == m_rn):
        flag = "***"
      line = "{:4d} {:8s} {:8.2f} {:8.2f} {:3s}".format(relaxing[rn]['config_id'], relaxing[rn]['structure'], relaxing[rn]['a0'], relaxing[rn]['e0'], flag)
      output.log(line, verbose=0)

    # Get config id
    config_id = relaxing[m_rn]['config_id']

    # Equilibrium density
    output.log("Get equilibirum density.", verbose=0)
    density = configs.calc_density(config_id)    
    rho_eq = sum(density[:,0]) / g.configs[config_id]['count']
    output.log("Rho eq " + str(rho_eq), verbose=0)
    S = 1.0 / rho_eq
    
    # F gradient at 
    potential.tab_for_output(rcut_none_set=(10.0))

    fgrad = None
    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      if(pfn['ftype'] == 3):
        fgrad = interp.interpolate(rho_eq, g.potential['tab_for_output'][fn][:,0] , g.potential['tab_for_output'][fn][:,2], 5, 0)
        output.log("F grad at Rho eq " + str(fgrad), verbose=0)
    C = 0
    if(fgrad is not None):
      C = -fgrad

    g.potential['tab_transformed'] = []
    for fn in range(len(g.potential['functions'])):
      g.potential['tab_transformed'].append(None)

    pk = None
    dk = None
    ek = None

    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      if(pfn['ftype'] == 1):
        g.potential['tab_transformed'][fn] = numpy.copy(g.potential['tab_for_output'][fn]) 
        pk = fn
      elif(pfn['ftype'] == 2):
        g.potential['tab_transformed'][fn] = numpy.copy(g.potential['tab_for_output'][fn])
        dk = fn
      elif(pfn['ftype'] == 3):
        g.potential['tab_transformed'][fn] = numpy.copy(g.potential['tab_for_output'][fn])
        ek = fn

    
    # Apply C
    g.potential['tab_transformed'][pk][:,1] = g.potential['tab_transformed'][pk][:,1] - 2 * C * g.potential['tab_transformed'][dk][:,1]
    g.potential['tab_transformed'][ek][:,1] = g.potential['tab_transformed'][ek][:,1] + C * g.potential['tab_transformed'][ek][:,0]


    # Apply S
    g.potential['tab_transformed'][dk][:,1] = g.potential['tab_transformed'][dk][:,1] * S
    g.potential['tab_transformed'][ek][:,0] = g.potential['tab_transformed'][ek][:,0] * S
    

    # Output and Plot
    output.potentials_transformed(g.dir['gauge_transformed'])
    plot.potentials(g.dir['gauge_transformed_plot'], source='transformed')

    #print(rho_eq)

    #density = numpy.zeros((g.configs[config_id]['count'], 100), dtype=numpy.float64)

