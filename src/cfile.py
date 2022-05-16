#!/bin/python3
"""
Config File
"""

import numpy
import os

from g import g
from ds import ds
from configs import configs
from label import label
from std import std


class cfile:

  @staticmethod
  def read_configs(config_tag=None):
    files = []
    for dir in g.input['configs']:
      files = files + cfile.file_list(dir)    
    for file in files:
      cfile.read(file, config_tag)


  @staticmethod
  def read(file_path, config_tag=None):

    #   Read file
    ###############################################################

    d = std.read(file_path)


    #   Read data from the file data
    ###############################################################

    c = ds.config()

    c['file_path'] = file_path

    # Default units
    l_units = 'ANG'
    e_units = 'EV'
    f_units = 'EV/ANG'
    s_units = 'EV/ANG3'
    cxyz = None
    epa = None
    epa_to_e = False
    calc_type = None

    labels = []
    coords = []
    forces = []

    for line in d:
      f = line.split(" ")
      if(f[0] == "#L_UNITS"):
        l_units = f[1]
      elif(f[0] == "#E_UNITS"):
        e_units = f[1]
      elif(f[0] == "#S_UNITS"):
        e_units = f[1]
      elif(f[0] == "#F_UNITS"):
        e_units = f[1]
      elif(f[0] == "#A0" or f[0] == "#ALAT"):
        c['a0'] = float(f[1])
      elif(f[0] == "#W"):
        c['w'] =  float(f[1])
      elif(f[0] == "#C"):
        cxyz = [int(f[1]), int(f[2]), int(f[3])]
      elif(f[0] == "#X"):
        c['uv'][0,0] = float(f[1])
        c['uv'][0,1] = float(f[2])
        c['uv'][0,2] = float(f[3])
      elif(f[0] == "#Y"):
        c['uv'][1,0] = float(f[1])
        c['uv'][1,1] = float(f[2])
        c['uv'][1,2] = float(f[3])
      elif(f[0] == "#Z"):
        c['uv'][2,0] = float(f[1])
        c['uv'][2,1] = float(f[2])
        c['uv'][2,2] = float(f[3])
      elif(f[0] == "#SX"):
        if(c['stress'] == None):
          c['stress'] = numpy.zeros((3,3,),dtype=numpy.float64)
        c['stress'][0,0] = float(f[1])
        c['stress'][0,1] = float(f[2])
        c['stress'][0,2] = float(f[3])
      elif(f[0] == "#SY"):
        c['stress'][1,0] = float(f[1])
        c['stress'][1,1] = float(f[2])
        c['stress'][1,2] = float(f[3])
      elif(f[0] == "#SZ"):
        c['stress'][2,0] = float(f[1])
        c['stress'][2,1] = float(f[2])
        c['stress'][2,2] = float(f[3])
      elif(f[0] == "#E"):
        epa = float(f[1])
        epa_in = False
      elif(f[0] == "#EPA"):
        epa = float(f[1])
        epa_in = True
      elif(f[0] == "#CALCTYPE"):
        calc_type = f[1].lower().strip()
      elif(f[0] == "#RCUT"):
        c['rcut'] =  float(f[1])


      else:
        if(len(f) == 4 or len(f) == 7):
          labels.append(label.set(f[0]))
          coords.append([float(f[1]), float(f[2]), float(f[3])])
        if(len(f) == 7):
          forces.append([float(f[4]), float(f[5]), float(f[6])]) 

    if(epa is not None and epa_in == False):
      epa = epa / len(labels)

    if(cxyz is not None):
      m = numpy.zeros((3,3,),dtype=numpy.float64)
      m[0,0] = cxyz[0]
      m[1,1] = cxyz[1]
      m[2,2] = cxyz[2]
      c['uv'] = numpy.matmul(m, c['uv'])

      e_labels = []
      e_coords = []
      e_forces = []

      for x in range(cxyz[0]):
        for y in range(cxyz[1]):
          for z in range(cxyz[2]):
            for n in range(len(labels)):
              e_labels.append(labels[n])
              cx = (x + coords[n][0]) / cxyz[0]
              cy = (y + coords[n][1]) / cxyz[1]
              cz = (z + coords[n][2]) / cxyz[2]
              e_coords.append([cx, cy, cz])
              if(len(forces) > 0):
                e_forces.append([forces[n][0], forces[n][1], forces[n][2]])

      c['labels'] = ds.config_labels(e_labels)
      c['coords'] = ds.config_coords(e_coords)   
      if(len(e_forces) > 0):
        c['forces'] = ds.config_forces(e_forces)  
           
    else:
      c['labels'] = ds.config_labels(labels)
      c['coords'] = ds.config_coords(coords)   
      if(len(forces) > 0):
        c['forces'] = ds.config_forces(forces)     
      cxyz = [1,1,1]

    c['cxyz'] = numpy.asfortranarray(cxyz, dtype=numpy.int32)
    c['count'] = len(c['labels'])
    c['energy'] = epa * c['count']

    if(calc_type is not None):
      c['calc_type'] = calc_type
    
    configs.add(c, config_tag)

    return c['config_id']


  @staticmethod
  def file_list(path_in, files=[]):
    for path in os.listdir(path_in):
      path_new = os.path.join(path_in, path)
      if(os.path.isdir(path_new)):
        files = cfile.file_list(path_new, files)
      else:
        files.append(path_new)
    return files


  @staticmethod
  def one_space(inp):
    out = ''
    last = None
    for c in inp:
      if(not (last == " " and c == " ")):
        out = out + c
      last = c
    return out