#!/bin/python3
"""
rss
"""

import numpy
import os
import time

from f_toolbox import math

from g import g
from ds import ds

from output import output
from configs import configs
from bp import bp


class rss:

  log = []

  @staticmethod
  def run():
    output.log("RSS", verbose=0)

    configs.calc('efs')
    efs_rss = configs.rss(config_tag='efs')
    bp_rss = bp.run_rss()
    total_rss = efs_rss + bp_rss

    
        
    output.log("EFS ", verbose=0)
    output.log("##################################################", verbose=0)
    output.log("", verbose=0)
    for config_id in range(len(g.configs)):
      if(g.configs[config_id]['tag'].lower() == 'efs'):
        output.log("{0:20s}{1:4d}".format("Config id:",config_id), verbose=0)
        output.log("{0:64s}".format(os.path.basename(g.configs[config_id]['file_path'])), verbose=0)
        output.log("{0:20s}".format("Energy"), verbose=0)
        output.log("{0:20s}{1:12.6f}".format("Energy input:",g.configs[config_id]['energy']), verbose=0)
        output.log("{0:20s}{1:12.6f}".format("Energy calc:",g.configs[config_id]['c_energy']), verbose=0)
        output.log("{0:20s}".format("RSS"), verbose=0)
        output.log("{0:20s}{1:12.6f}".format("Energy rss:",g.configs[config_id]['rss_details']['energy']), verbose=0)
        output.log("{0:20s}{1:12.6f}".format("Force rss:",g.configs[config_id]['rss_details']['forces']), verbose=0)
        output.log("{0:20s}{1:12.6f}".format("Stress rss:",g.configs[config_id]['rss_details']['stress']), verbose=0)
        output.log("{0:20s}{1:12.6f}".format("Total rss:",g.configs[config_id]['rss']), verbose=0)
        output.log("", verbose=0)
    output.log("", verbose=0)

    output.log("BP ", verbose=0)
    output.log("##################################################", verbose=0)
    output.log("", verbose=0)
    for n in range(len(g.bp)):
      output.log("{:20s}{:4d}".format("BP id:",n), verbose=0)
      output.log("{:20s}{:12.6f} ({:12.6f}|{:12.6f})".format("a0 rss:",g.bp[n]['rss_details']['a0'],g.bp[n]['a0'],g.bp[n]['c_a0']), verbose=0)
      output.log("{:20s}{:12.6f} ({:12.6f}|{:12.6f})".format("e0 rss:",g.bp[n]['rss_details']['e0'],g.bp[n]['e0'],g.bp[n]['c_e0']), verbose=0)
      output.log("{:20s}{:12.6f} ({:12.6f}|{:12.6f})".format("b0 rss:",g.bp[n]['rss_details']['b0'],g.bp[n]['b0_gpa'],g.bp[n]['c_b0_gpa']), verbose=0)
      output.log("{:20s}{:12.6f}".format("ec rss:",g.bp[n]['rss_details']['ec']), verbose=0)
      output.log("{:20s}{:12.6f}".format("rose rss:",g.bp[n]['rss_details']['rose']), verbose=0)
      output.log("{:20s}{:12.6f}".format("msp shape rss:",g.bp[n]['rss_details']['msp_shape']), verbose=0)
      output.log("{:20s}{:12.6f}".format("msp values rss:",g.bp[n]['rss_details']['msp_values']), verbose=0)
      output.log("{:20s}{:12.6f}".format("bm eos shape rss:",g.bp[n]['rss_details']['bm_eos_shape']), verbose=0)
      output.log("{:20s}{:12.6f}".format("bm eos values rss:",g.bp[n]['rss_details']['bm_eos_values']), verbose=0)
      output.log("{:20s}{:12.6f}".format("penalty rss:",g.bp[n]['rss_details']['bp_penalty']), verbose=0)
      output.log("{:20s}{:12.6f}".format("total rss:",g.bp[n]['rss']), verbose=0)
      
      output.log("", verbose=0)
    output.log("", verbose=0)

 
    output.log("Summary ", verbose=0)
    output.log("##################################################", verbose=0)
    output.log("EFS rss:     " + str(efs_rss), verbose=0)
    output.log("BP rss:      " + str(bp_rss), verbose=0)
    output.log("Total rss:   " + str(total_rss), verbose=0)
    output.log("", verbose=0)

    




