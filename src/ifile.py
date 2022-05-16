#!/bin/python3
"""
Input File
"""

import os
import numpy

from f_toolbox import math

from g import g
from ds import ds
from output import output
from std import std
from units import units

class ifile:


  @staticmethod
  def load(file_path):
    output.log("Read input file", verbose=0)

    #   Default units
    ###############################################################
    input_units = {'energy': 'ev', 'length': 'ang', 'pressure': 'gpa',}

    #   Read file
    ###############################################################

    d = std.read(file_path)


    #   Read data from the file data
    ###############################################################
    g.input = ds.ifile()
    g.runtypes = ds.runtypes()

    for line in d:
      f = line.split(" ")
      f[0] = f[0].upper()



      if(f[0] == "#CONFIGS"):
        dir = os.path.abspath(f[1])
        if(os.path.isdir(dir)):
          g.input['configs'].append(dir)
        else:
          output.log("Directory does not exist: " + dir, verbose=0)



      elif(f[0] == "#CONFIG"):
        g.one_line_configs.append(f[1:])



      elif(f[0] == "#TYPE"):
        runtype = f[1].strip().lower()
        if(runtype in g.runtypes):
          g.input['runtype'] = runtype
        else:
          output.log("Unknown runtype: " + runtype, verbose=0)




      elif(f[0] == "#NLSIZE"):
        g.input['nlsize'] = int(f[1])



      elif(f[0] == "#POTENTIAL"):
        potfile = os.path.abspath(f[1])
        if(os.path.isfile(potfile)):
          g.input['potindex'] = potfile
        else:
          output.log("Potential index does not exist: " + potfile, verbose=0)
      

      
      elif(f[0] == "#RSS"):
        for n in range(1,len(f)):
          df = f[n].split("=")
          if(len(df) == 2):
            if(df[0].lower() == "configs"):
              g.rss_w['configs'] = float(df[1])
            elif(df[0].lower() == "e"):
              g.rss_w['e'] = float(df[1])
            elif(df[0].lower() == "f"):
              g.rss_w['f'] = float(df[1])
            elif(df[0].lower() == "s"):
              g.rss_w['s'] = float(df[1])
            elif(df[0].lower() == "bp"):
              g.rss_w['bp'] = float(df[1])
            elif(df[0].lower() == "a0"):
              g.rss_w['a0'] = float(df[1])
            elif(df[0].lower() == "e0"):
              g.rss_w['e0'] = float(df[1])
            elif(df[0].lower() == "b0"):
              g.rss_w['b0'] = float(df[1])
            elif(df[0].lower() == "ec"):
              g.rss_w['ec'] = float(df[1])
            elif(df[0].lower() == "rose"):
              g.rss_w['rose'] = float(df[1])
            elif(df[0].lower() == "msp"):
              v = numpy.asarray(df[1].split(","), dtype=numpy.float64)
              g.rss_w['msp'] = v
            elif(df[0].lower() == "bm_eos"):
              v = numpy.asarray(df[1].split(","), dtype=numpy.float64)
              g.rss_w['bm_eos'] = v
      


      elif(f[0] == "#BP"):
        bp = ds.bp()
        atoms_per_cell = 1
        for n in range(1, len(f[:])):
          df = f[n].split("=")
          if(len(df) == 2):
            if(df[0].lower() == "label"):
              bp['label'] = df[1].split(",")
            elif(df[0].lower() == "structure"):
              bp['structure'] = df[1]
              if(bp['structure'].lower() == 'fcc' or bp['structure'].lower() == 'c-fcc'):
                atoms_per_cell = 4
              elif(bp['structure'].lower() == 'bcc' or bp['structure'].lower() == 'c-bcc'):
                atoms_per_cell = 2
              elif(bp['structure'].lower() == 'sc' or bp['structure'].lower() == 'c-sc'):
                atoms_per_cell = 1
            elif(df[0].lower() == "a0"):
              bp['a0'] = float(df[1])       
            elif(df[0].lower() == "uv"):
              uv_list = df[1].split(",")      
              bp['uv'] = ds.config_uv(uv_list)
            elif(df[0].lower() == "e0"):
              bp['e0'] = float(df[1]) 
            elif(df[0].lower() == "b0"):
              b0_in = float(df[1]) 
              bp['b0'] = units.convert(input_units['pressure'], 'EV/ANG3', b0_in)
              bp['b0_gpa'] = units.convert(input_units['pressure'], 'GPA', b0_in)
            elif(df[0].lower() == "ec"):
              ec_in = ds.bp_ec(df[1].split(",")) 
              bp['ec'] = ds.bp_ec()
              bp['ec_gpa'] = ds.bp_ec() 
              for i in range(6):
                for j in range(6):
                  bp['ec'][i,j] = units.convert(input_units['pressure'], 'EV/ANG3', ec_in[i, j])
                  bp['ec_gpa'][i,j] = units.convert(input_units['pressure'], 'GPA', ec_in[i, j])
        if(bp['a0'] is not None and bp['uv'] is not None): 
          bp['v0'] = math.cellvolume(bp['a0'], bp['uv']) / atoms_per_cell
        g.bp.append(bp)
      


      elif(f[0] == "#FIT"):
        step = ds.potfit_step()
        for n in range(1, len(f[:])):
          df = f[n].split("=")
          if(len(df) == 2):
            if(df[0].lower() == "type"):
              step['type'] = df[1].lower() 
            elif(df[0].lower() == "niter"):
              step['niter'] = int(df[1]) 
            elif(df[0].lower() == "titer"):
              step['titer'] = int(df[1]) 
            elif(df[0].lower() == "tstart"):
              step['tstart'] = float(df[1]) 
            elif(df[0].lower() == "tend"):
              step['tend'] = float(df[1]) 
            elif(df[0].lower() == "pfact"):
              step['pfact'] = float(df[1]) 
            elif(df[0].lower() == "pvar"):
              step['pvar'] = float(df[1]) 
            elif(df[0].lower() == "vartype"):
              step['vartype'] = str(df[1]).strip().lower() 
            elif(df[0].lower() == "gaussian"):
              step['gaussian']  = False
              if(str(df[1]).strip().lower()[0] == 't'):
                step['gaussian']  = True
            elif(df[0].lower() == "popsize"):
              step['popsize'] = int(df[1]) 
            elif(df[0].lower() == "fresh"):
              step['fresh'] = float(df[1]) 
            elif(df[0].lower() == "search"):
              step['search'] = df[1].lower().strip()
        g.potfit_steps.append(step) 



      elif(f[0] == "#MAXTIME"):
        maxtime = float(f[1].strip())
        g.maxtime = maxtime



      elif(f[0] == "#BESTSAVEPERIOD"):
        tp = float(f[1].strip())
        g.potfit_bestsaveperiod = tp







  def one_space(inp):
    out = ''
    last = None
    for c in inp:
      if(not (last == " " and c == " ")):
        out = out + c
      last = c
    return out

