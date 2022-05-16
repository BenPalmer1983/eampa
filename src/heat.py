#!/bin/python3
"""
heat
"""

import numpy
import os

from g import g
from ds import ds
from std import std
from output import output
from label import label
from configs import configs
from nl import nl



class heat:

  @staticmethod
  def gauss(config_id, heat_max=0.01, rebuild_nl=True):
    count = g.configs[config_id]['count']
    mu, sigma = 0, 1.0 
    r = heat_max * (0.5 - numpy.random.normal(mu, sigma, 3 * count))
    g.configs[config_id]['coords'][:,0] = (g.configs[config_id]['coords'][:,0] +  r[:count]) % 1.0
    g.configs[config_id]['coords'][:,1] = (g.configs[config_id]['coords'][:,1] +  r[count:2*count]) % 1.0
    g.configs[config_id]['coords'][:,2] = (g.configs[config_id]['coords'][:,2] +  r[2*count:]) % 1.0
    
    if(rebuild_nl):
      nl.rebuild_nl(config_id)
