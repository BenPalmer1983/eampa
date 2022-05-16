#!/bin/python3
"""
efs
"""

import numpy
import os

from g import g
from ds import ds
from std import std
from output import output
from label import label
from nl import nl
from configs import configs


class efs:


  @staticmethod
  def run():
    
    output.log("Run EFS", verbose=0)
        
    # make additional directories
    efs.make_dirs()
        
    # run calculations
    nl.build('efs')     
    configs.calc('efs')
    configs.rss(config_tag='efs')


    # output configs
    output.configs(g.dir['efs_configs'])
      
    # output forces
    output.config_forces(g.dir['efs_forces'])

    # Output summary
    output.log("", verbose=0)
    for config_id in range(len(g.configs)):
      if(g.configs[config_id]['tag'] == 'efs'): 
        output.config_summary(config_id)

    
    output.config_nls(g.dir['efs_nls'])
    
    # exit
    output.log("End", verbose=0)
    exit()
    


  # Make additional directories for efs
  @staticmethod
  def make_dirs():    

    g.dir['efs_configs'] = os.path.join(g.dir['efs'], "configs")
    g.dir['efs_forces'] = os.path.join(g.dir['efs'], "forces")
    g.dir['efs_nls'] = os.path.join(g.dir['efs'], "nls")

    for k in g.dir.keys():
      os.makedirs(g.dir[k], exist_ok=True)




















