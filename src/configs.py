#!/bin/python3
"""
Configs
"""

import numpy
import os
import time

from f_toolbox import atom
from f_toolbox import math
from f_toolbox import transforms

from g import g
from ds import ds
from std import std
from output import output
from label import label
from nl import nl
from units import units


class configs:

  log = []
  calc_counter = 0
  calc_time = 0.0

  @staticmethod
  def add(c, config_tag=None):    
    c['config_id'] = len(g.configs)
    c['tag'] = config_tag
    c['volume'] = c['a0']**3 * math.tripleproduct(c['uv'][:,0], c['uv'][:,1], c['uv'][:,2])

    # Add empty arrays for computed stress and forces
    c['c_forces'] = ds.config_forces_empty(c['count'])
    c['c_stress'] = ds.config_stress_empty()
    c['c_stress_gpa'] = ds.config_stress_empty()
    c['coords'] = c['coords'] % 1.0

    g.configs.append(c)
    if(c['tag'] is None):
      tag = ""
    else:
      tag = str(c['tag'])
    output.log("Add config " + str(c['config_id']) + " " + tag, verbose=1)
    return  c['config_id']


  @staticmethod
  def calc(config_tag='all'):
    configs.log = []
    if(isinstance(config_tag, int)):
      configs.calc_efs(config_tag)
    else:
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
          configs.calc_efs(config_id)
    output.log(configs.log, verbose=2)
    configs.log = []


  @staticmethod
  def calc_efs(config_id):
    start_time = time.time()

    # Get IDs
    nl_id = g.configs[config_id]['nl_id']
    if(nl_id is None):
      nl.build_nl(config_id) 
      nl_id = g.configs[config_id]['nl_id']
      

    csize = g.configs[config_id]['count']    
    calc_type = g.configs[config_id]['calc_type']

    if(calc_type is None):
      calc_type = 'e'
      if(g.configs[config_id]['forces'] is not None):
        calc_type = 'ef'
      if(g.configs[config_id]['stress'] is not None):
        calc_type = 'efs'
        
    if(calc_type == 'e'):
      g.stats['config_calc_count'] = g.stats['config_calc_count'] + 1
      atom.e(g.nl[nl_id]['labels'][:,:], 
            g.nl[nl_id]['r'][:], 
            g.configs[config_id]['labels'][:], 
            g.configs[config_id]['c_pair_embedding_energy'][:])   
      g.configs[config_id]['c_energy'] = g.configs[config_id]['c_pair_embedding_energy'][2]
      g.configs[config_id]['last_calc_type'] = 'e'
      configs.log.append("Config " + str(config_id) + "    e calc   " + str(g.configs[config_id]['c_energy']))
    
    elif(calc_type == 'ef'):
      g.stats['config_calc_count'] = g.stats['config_calc_count'] + 1
      atom.ef(g.nl[nl_id]['labels'][:,:], 
            g.nl[nl_id]['r'][:], 
            g.nl[nl_id]['rvec'][:, :], 
            g.configs[config_id]['labels'][:], 
            g.configs[config_id]['c_pair_embedding_energy'][:], 
            g.configs[config_id]['c_forces'][:,:])
      g.configs[config_id]['c_energy'] = g.configs[config_id]['c_pair_embedding_energy'][2] 
      g.configs[config_id]['c_forces_total'] = numpy.sum(numpy.absolute(g.configs[config_id]['c_forces'][:,:]))
      g.configs[config_id]['last_calc_type'] = 'ef'
      configs.log.append("Config " + str(config_id) + "    ef calc  " + str(g.configs[config_id]['c_energy']) + "  " + str(g.configs[config_id]['c_forces_total']))

    #c['volume']
    elif(calc_type == 'efs'):
      g.stats['config_calc_count'] = g.stats['config_calc_count'] + 1
      atom.efs(g.nl[nl_id]['labels'][:,:], 
            g.nl[nl_id]['r'][:], 
            g.nl[nl_id]['rvec'][:, :], 
            g.nl[nl_id]['inhalo'][:],
            g.configs[config_id]['labels'][:], 
            g.configs[config_id]['volume'], 
            g.configs[config_id]['c_pair_embedding_energy'][:], 
            g.configs[config_id]['c_forces'][:,:], 
            g.configs[config_id]['c_stress'][:,:])
      for i in range(3):
        for j in range(3): 
          g.configs[config_id]['c_stress_gpa'][i,j] =  units.convert('EV/ANG3','GPA',g.configs[config_id]['c_stress'][i,j])
      g.configs[config_id]['c_energy'] = g.configs[config_id]['c_pair_embedding_energy'][2]
      g.configs[config_id]['c_forces_total'] = numpy.sum(numpy.absolute(g.configs[config_id]['c_forces'][:,:]))
      g.configs[config_id]['last_calc_type'] = 'efs'
      configs.log.append("Config " + str(config_id) + "    efs calc " + str(g.configs[config_id]['c_energy']) + "  " + str(g.configs[config_id]['c_forces_total']))

    calc_time = (time.time() - start_time)
    g.configs[config_id]['timer'] = g.configs[config_id]['timer'] + calc_time
    g.configs[config_id]['counter'] = g.configs[config_id]['counter'] + 1

    g.calc_counter = g.calc_counter + 1
    g.calc_time = g.calc_time + calc_time
   
  @staticmethod
  def one_line_configs():
    
    for c in g.one_line_configs:
      structure = None
      label = None
      a0 = None
      uv = None
      cxyz = None
      rcut = 6.5
      tag = 'efs'
      ctype = 'e'
      setcoords = []
      changecoords = []

      for n in range(len(c)):
        df = c[n].split("=")
        if(len(df) == 2):
          if(df[0].lower() == "label"):
            label = df[1].split(",")
          elif(df[0].lower() == "structure"):
            structure = df[1]
          elif(df[0].lower() == "a0"):
            a0 = float(df[1])
          elif(df[0].lower() == "uv"):
            uv = ds.config_uv(df[1].split(","))
          elif(df[0].lower() == "c"):
            cxyz = ds.config_cxyz(df[1].split(","))
          elif(df[0].lower() == "rcut"):
            rcut = float(df[1])
          elif(df[0].lower() == "tag"):
            tag = df[1].strip().lower()
          elif(df[0].lower() == "ctype"):
            ctype = df[1].strip().lower()
          elif(df[0].lower() == "setcoord"):
            coords = df[1].split(",")
            setcoords.append([int(coords[0]), float(coords[1]) % 1.0, float(coords[2]) % 1.0, float(coords[3]) % 1.0])
          elif(df[0].lower() == "changecoord"):
            coords = df[1].split(",")
            changecoords.append([int(coords[0]), float(coords[1]), float(coords[2]), float(coords[3])])

      configs.add_common(structure=structure, labels=label, a0=a0, uv=uv, cxyz=cxyz, rcut=rcut, index_tag=tag, calc_type=ctype, setcoords=setcoords, changecoords=changecoords)

  @staticmethod
  def add_common(structure, labels, a0, uv=[1], cxyz=[1,1,1], rcut=6.5, index_tag=None, calc_type='efs', setcoords=None, changecoords=None): 
    c = ds.config()

    # Unit vector
    uv = ds.config_uv(uv)
    cxyz = ds.config_cxyz(cxyz)

    c['a0'] = a0
    c['rcut'] = rcut
    m = numpy.zeros((3,3,),dtype=numpy.float64)
    m[0,0] = cxyz[0]
    m[1,1] = cxyz[1]
    m[2,2] = cxyz[2]
    c['uv'] = numpy.matmul(m, uv)
    c['cxyz'] = numpy.asfortranarray(cxyz, numpy.int32)
    coords_prim = configs.structure(structure)

    c['count'] = cxyz[0] * cxyz[1] * cxyz[2] * len(coords_prim)
    c['labels'] = ds.config_labels_empty(c['count'])
    c['coords'] = ds.config_coords_empty(c['count'])
    

    m = 0
    for x in range(cxyz[0]):
      for y in range(cxyz[1]):
        for z in range(cxyz[2]):
          for n in range(len(coords_prim)):
            c['labels'][m] = label.set(labels[n%len(labels)])
            c['coords'][m, 0] = (x + coords_prim[n][0]) / cxyz[0]
            c['coords'][m, 1] = (y + coords_prim[n][1]) / cxyz[1]
            c['coords'][m, 2] = (z + coords_prim[n][2]) / cxyz[2]
            m = m + 1     

    c['calc_type'] = calc_type

    # Set coords
    if(setcoords is not None):
      for n in range(len(setcoords)):
        m = setcoords[n][0]
        x = setcoords[n][1]
        y = setcoords[n][2]
        z = setcoords[n][3]
        c['coords'][m, 0] = x
        c['coords'][m, 1] = y
        c['coords'][m, 2] = z

    # Change coords
    if(changecoords is not None):
      for n in range(len(changecoords)):
        m = changecoords[n][0]
        x = changecoords[n][1]
        y = changecoords[n][2]
        z = changecoords[n][3]
        c['coords'][m, 0] = (c['coords'][m, 0] + x) % 1.0
        c['coords'][m, 1] = (c['coords'][m, 1] + y) % 1.0
        c['coords'][m, 2] = (c['coords'][m, 2] + z) % 1.0



    # Add
    config_id = configs.add(c, index_tag)
    return  config_id


  @staticmethod
  def strain(config_id, distortion=None, strain=[], invert=False):

    if(distortion == None):
      return False

    if(not isinstance(strain, list)):
      strain = [strain]
    
    if(distortion.lower() == 'd1'):
      if(invert):
        uvt = transforms.d1_inv(g.configs[config_id]['uv'], strain[0])
      else:
        uvt = transforms.d1(g.configs[config_id]['uv'], strain[0])
    elif(distortion.lower() == 'd2'):
      if(invert):
        uvt = transforms.d2_inv(g.configs[config_id]['uv'], strain[0])
      else:
        uvt = transforms.d2(g.configs[config_id]['uv'], strain[0])
    elif(distortion.lower() == 'd3'):
      if(invert):
        uvt = transforms.d3_inv(g.configs[config_id]['uv'], strain[0])
      else:
        uvt = transforms.d3(g.configs[config_id]['uv'], strain[0])


    elif(distortion.lower() == 'ctd1'):
      uvt = transforms.ctd1(g.configs[config_id]['uv'], strain[0])
      



    # RFKJ 1998 strain 1
    #elif(distortion == 'd1'):
    #  d1(uv, s)
    #minverse3(a, c)
    #transforms
    nl_id = g.configs[config_id]['nl_id']
    if(nl_id is None):
      g.configs[config_id]['uv'] = uvt
    else:
      nl.change_cell(config_id, uv=uvt, scale_uv=False)
      

  @staticmethod
  def count(config_tag='all'):
    if(isinstance(config_tag, int)):
      return 1
    else:
      counter = 0
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
          counter = counter + 1
    return counter



  @staticmethod
  def rss(config_tag='all'):
    rss_total = 0.0
    if(isinstance(config_tag, int)):
      rss_total = configs.rss_inner(config_tag)
    else:
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
          rss_total = rss_total + configs.rss_inner(config_id)
    output.log("EFS rss  " + config_tag + "   " + str(rss_total), verbose=2)
    return rss_total
   


  @staticmethod
  def rss_inner(config_id):
    # Set start value
    g.configs[config_id]['rss_details'] = ds.config_rss_details()
    rss_val = 0.0
   
    # Energy
    if(g.configs[config_id]['energy'] is not None and g.configs[config_id]['c_energy'] is not None):
      g.configs[config_id]['rss_details']['energy'] = g.rss_w['e'] * (g.configs[config_id]['energy'] - g.configs[config_id]['c_energy'])**2 

    # Forces
    if(g.configs[config_id]['forces'] is not None and g.configs[config_id]['c_forces'] is not None):
      g.configs[config_id]['rss_details']['forces'] = g.rss_w['f'] * sum(sum((g.configs[config_id]['forces'] - g.configs[config_id]['c_forces'])**2)) 

    # Stress
    if(g.configs[config_id]['stress'] is not None and g.configs[config_id]['c_stress'] is not None):
      g.configs[config_id]['rss_details']['stress'] = g.rss_w['s'] * sum(sum((g.configs[config_id]['stress'] - g.configs[config_id]['c_stress'])**2)) 

    # Overall configs weighting + total
    t = 0.0
    for k in g.configs[config_id]['rss_details'].keys():
      g.configs[config_id]['rss_details'][k] = g.configs[config_id]['rss_details'][k] * g.rss_w['configs']
      t = t + g.configs[config_id]['rss_details'][k]
    g.configs[config_id]['rss_details']['config_total'] = t
    g.configs[config_id]['rss'] = t

    return g.configs[config_id]['rss']


  @staticmethod
  def duplicate(src_id, dest_id=None):
    
    c = ds.config()
    c['tag'] = g.configs[src_id]['tag'] 
    c['file_path'] = g.configs[src_id]['file_path'] 
    c['a0'] = g.configs[src_id]['a0'] 
    c['uv'] = g.configs[src_id]['uv'] 
    c['cxyz'] = g.configs[src_id]['cxyz'] 
    c['volume'] = g.configs[src_id]['volume']     
    c['w'] = g.configs[src_id]['w'] 
    c['rcut'] = g.configs[src_id]['rcut'] 
    c['count'] = g.configs[src_id]['count'] 
    c['labels'] = numpy.copy(g.configs[src_id]['labels']) 
    c['coords'] = numpy.copy(g.configs[src_id]['coords'])
    c['energy'] = g.configs[src_id]['energy'] 
    c['forces'] = numpy.copy(g.configs[src_id]['forces']) 
    c['stress'] = numpy.copy(g.configs[src_id]['stress']) 
    c['c_energy'] = g.configs[src_id]['c_energy'] 
    c['c_forces'] = numpy.copy(g.configs[src_id]['c_forces']) 
    c['c_stress'] = numpy.copy(g.configs[src_id]['c_stress']) 
    c['c_pair_embedding_energy'] = numpy.copy(g.configs[src_id]['c_pair_embedding_energy'])
    c['c_forces_total'] = g.configs[src_id]['c_forces_total'] 
    c['calc_type'] = g.configs[src_id]['calc_type'] 
    c['last_calc_type'] = g.configs[src_id]['last_calc_type'] 
    c['rss'] = g.configs[src_id]['rss'] 
    # Ensure coords are 0.0 to 1.0
    c['coords'] = c['coords'] % 1.0

    if(dest_id is None):
      # Get new config_id
      config_id = len(g.configs)
      c['config_id'] = config_id
      g.configs.append(c)
      output.log("Duplicate config: " + str(src_id) + " to " + str(config_id), verbose=3)
      # Also duplicate neighbour list
      if(g.configs[src_id]['nl_id'] is not None):
        nl_id = nl.duplicate(g.configs[src_id]['nl_id']) 
        g.configs[config_id]['nl_id'] = nl_id
        g.nl[nl_id]['config_id'] = config_id
      return config_id
    elif(dest_id < len(g.configs)):
      nl_id_dest = g.configs[dest_id]['nl_id'] 
      nl_id_src = g.configs[src_id]['nl_id'] 
      config_id = dest_id
      c['config_id'] = config_id

      g.configs[config_id] = c
      output.log("Duplicate config: " + str(src_id) + " to " + str(config_id), verbose=3)
      # Also duplicate neighbour list
      if(nl_id_src is not None):
        nl_id = nl.duplicate(nl_id_src, nl_id_dest)         
        g.configs[config_id]['nl_id'] = nl_id
        g.nl[nl_id]['config_id'] = config_id

      return config_id
    else:
      return None

  """
  @staticmethod
  def copy(src_id, dest_id=None):
    
    c = ds.config()
    c['tag'] = g.configs[src_id]['tag'] 
    c['file_path'] = g.configs[src_id]['file_path'] 
    c['a0'] = g.configs[src_id]['a0'] 
    c['uv'] = g.configs[src_id]['uv'] 
    c['cxyz'] = g.configs[src_id]['cxyz'] 
    c['volume'] = g.configs[src_id]['volume']     
    c['w'] = g.configs[src_id]['w'] 
    c['rcut'] = g.configs[src_id]['rcut'] 
    c['count'] = g.configs[src_id]['count'] 
    c['labels'] = numpy.copy(g.configs[src_id]['labels']) 
    c['coords'] = numpy.copy(g.configs[src_id]['coords'])
    c['energy'] = g.configs[src_id]['energy'] 
    c['forces'] = numpy.copy(g.configs[src_id]['forces']) 
    c['stress'] = numpy.copy(g.configs[src_id]['stress']) 
    c['c_energy'] = g.configs[src_id]['c_energy'] 
    c['c_forces'] = numpy.copy(g.configs[src_id]['c_forces']) 
    c['c_stress'] = numpy.copy(g.configs[src_id]['c_stress']) 
    c['c_pair_embedding_energy'] = numpy.copy(g.configs[src_id]['c_pair_embedding_energy'])
    c['c_forces_total'] = g.configs[src_id]['c_forces_total'] 
    c['calc_type'] = g.configs[src_id]['calc_type'] 
    c['last_calc_type'] = g.configs[src_id]['last_calc_type'] 
    c['rss'] = g.configs[src_id]['rss'] 

    c['coords'] = c['coords'] % 1.0

    if(dest_id == None):
      dest_id = len(g.configs)
      c['config_id'] = dest_id

      # Also copy neighbour list
      if(g.configs[src_id]['nl_id'] is not None):
        print(g.configs[src_id]['nl_id']) 

      g.configs.append(c)

    output.log("Copied config " + str(src_id) + " to " + str(dest_id), verbose=1)
    """  






  @staticmethod
  def structure(structure):
    if(structure=='fcc'):
      return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    elif(structure=='bcc'):
      return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    elif(structure=='c-fcc'):
      return numpy.asarray([[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]])
    elif(structure=='c-bcc'):
      return numpy.asarray([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]])
    elif(structure=='c-sc'):
      return numpy.asarray([[0.5, 0.5, 0.5]])
    elif(structure=='c-hcp'):
      return numpy.asarray([[0.3333333,0.1666667,0.2500000],[0.6666667,0.8333333,0.7500000]])


  """
  @staticmethod
  def structure(structure):
    if(structure=='fcc'):
      return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    if(structure=='bcc'):
      return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
  """
    
  @staticmethod
  def output(tag='all'):
    for config_id in range(len(g.configs)):
      if(tag == 'all' or g.configs[config_id]['tag'] == tag):
        output.config(config_id, g.dir['data_configs'], verbose=1)


    
  @staticmethod
  def index(tag='all'):
    output.log(" ", verbose=0)
    output.log("###################### ", verbose=0)
    output.log("    Config Index       ", verbose=0)
    output.log("###################### ", verbose=0)
    for config_id in range(len(g.configs)):
      if(tag == 'all' or g.configs[config_id]['tag'] == tag):
        output.log("Config " + str(config_id), verbose=0)
        output.log("    count:     " + str(g.configs[config_id]['count']), verbose=0) 
        output.log("    c_energy:  " + str(g.configs[config_id]['c_energy']), verbose=0)
        if(g.configs[config_id]['nl_id'] is not None):
          output.log("    nl id:     " + str(g.configs[config_id]['nl_id']), verbose=0)
          output.log("    conf id:   " + str(g.nl[g.configs[config_id]['nl_id']]['config_id']), verbose=0)
          output.log("    count:     " + str(g.nl[g.configs[config_id]['nl_id']]['count']), verbose=0) 

   
    output.log("###################### ", verbose=0)
    output.log(" ", verbose=0)

