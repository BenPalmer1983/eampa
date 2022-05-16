#!/bin/python3

import numpy
import time
from f_toolbox import atom

from g import g
from ds import ds
from output import output



class nl:

  log = []

  @staticmethod
  def build(config_tag='all'):
    nl.log = []
    if(isinstance(config_tag, int)):
      nl_id = nl.build_nl(config_tag)
    else:
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
          nl_id = nl.build_nl(config_id)

    output.log(nl.log, verbose=1)
    nl.log = []

    """
    config_tag = config_tag.lower()
    if(config_tag == 'all'):
      for cid in range(len(g.configs)):
        nl_id = nl.build_nl(cid)
    else:
      for cid in range(len(g.configs)): 
        if(g.configs[cid]['tag'] == config_tag):
          nl_id = nl.build_nl(cid)
    """

  @staticmethod
  def build_nl(config_id, rebuild=False):
    c = g.configs[config_id]
    nl_id = c['nl_id']

    # If already built, and not forced to rebuild, exit
    if(nl_id is not None and rebuild == False):
      nl.log.append("Build NL " + str(c['nl_id']) + " - already built, skipping")
      return None
    # Get nl id
    if(nl_id == None):
      nl_id = len(g.nl)
      g.nl.append(ds.nl())
      c['nl_id'] = nl_id
    

    # Build neighbour list
    ta = time.time()
    s = g.input['nlsize']
    lout = ds.nl_lout(s)
    r = ds.nl_r(s)
    rvec = ds.nl_rvec(s)
    inhalo = ds.nl_inhalo(s)
    nlsize = 0
    nlsize = atom.nl(c['labels'], c['coords'], c['a0'], c['uv'], c['rcut'], lout, r, rvec, inhalo)
    tb = time.time() - ta

    # Store data
    g.nl[nl_id]['nl_id'] = nl_id
    g.nl[nl_id]['config_id'] = c['config_id']
    g.nl[nl_id]['count'] = nlsize
    g.nl[nl_id]['labels'] = lout[0:nlsize,:]
    g.nl[nl_id]['r'] = numpy.asfortranarray(r[0:nlsize])
    g.nl[nl_id]['rvec'] = numpy.asfortranarray(rvec[0:nlsize,:])
    g.nl[nl_id]['inhalo'] = numpy.asfortranarray(inhalo[0:nlsize])
    tc = time.time() - ta

    nl.log.append("Config " + str(config_id) + " - build NL " + str(c['nl_id']) + " count: " + str(nlsize) + "  "
               + "time: " + "{0:.4e}".format(tb) + "  " + "{0:.4e}".format(tc))
    return nl_id



  @staticmethod
  def rebuild_nl(config_id):
    return nl.build_nl(config_id, True)


  @staticmethod
  def update_coords(config_id):
    c = g.configs[config_id]
    nl_id = c['nl_id']

    # If not built, build it
    if(nl_id is None):
      return nl.build_nl(config_id)


    ta = time.time()
    atom.nl_update_coords(c['coords'], c['a0'], c['uv'], g.nl[nl_id]['labels'], 
                          g.nl[nl_id]['inhalo'], g.nl[nl_id]['r'], g.nl[nl_id]['rvec'])

    tb = time.time() - ta
    output.log("Config " + str(config_id) + " - update NL " + str(c['nl_id']) + "  " + "time: " + "{0:.4e}".format(tb), verbose=2)

    return nl_id

 
  @staticmethod
  def change_cell(config_id, a0=None, uv=None, scale_uv=False):
    a0_new = a0
    uv_new = uv

    # Make sure nl has been built already, but don't rebuild if already there
    nl.build_nl(config_id)

    # Break out if none set
    if(a0_new is None and uv_new is None):
      return False      

    # Process new a0
    if(a0_new is None):
      a0_new = g.configs[config_id]['a0']

    # Process new uv
    if(uv_new is None):
      uv_new = g.configs[config_id]['uv']
    else:
      uv_new = ds.config_uv(uv_new)
      if(scale_uv):
        m = numpy.zeros((3,3,),dtype=numpy.float64)
        m[0,0] = g.configs[config_id]['cxyz'][0]
        m[1,1] = g.configs[config_id]['cxyz'][1]
        m[2,2] = g.configs[config_id]['cxyz'][2]
        uv_new = numpy.matmul(m, uv_new)

    # Break out if no change
    if(a0_new == g.configs[config_id]['a0'] and numpy.array_equal(uv_new, g.configs[config_id]['uv'])):
      return False      

    nl_id = g.configs[config_id]['nl_id']
    
    ta = time.time()
    atom.nl_change_cell(a0_new, uv_new, g.configs[config_id]['a0'], g.configs[config_id]['uv'], g.nl[nl_id]['r'][:], g.nl[nl_id]['rvec'][:,:])
    tb = time.time() - ta


    output.log("Config " + str(config_id) + " - change r and rvec " + str(nl_id) + " count: " 
               + "time: " + "{0:.4e}".format(tb), verbose=3)
    
    # Update a0 and uv
    g.configs[config_id]['a0'] = a0_new
    g.configs[config_id]['uv'] = uv_new


  @staticmethod
  def duplicate(src_id, dest_id=None):
    nl_d = ds.nl()
    nl_d['config_id'] = None
    nl_d['count'] = g.nl[src_id]['count']
    nl_d['labels'] = numpy.copy(g.nl[src_id]['labels'])
    nl_d['r'] = numpy.copy(g.nl[src_id]['r'])
    nl_d['rvec'] = numpy.copy(g.nl[src_id]['rvec'])
    nl_d['inhalo'] = numpy.copy(g.nl[src_id]['inhalo'])
    
    if(dest_id is None):
      nl_id = len(g.nl)
      nl_d['nl_id'] = nl_id 
      g.nl.append(nl_d)
      output.log("Duplicate NL: " + str(src_id) + " to " + str(nl_id) + " (NEW)", verbose=3)
      return nl_id
    elif(dest_id < len(g.nl)):
      nl_id = dest_id
      nl_d['nl_id'] = nl_id 
      g.nl[dest_id] = nl_d
      output.log("Duplicate NL: " + str(src_id) + " to " + str(nl_id) + "", verbose=3)
      return nl_id
    else:
      return None

    """
    nl_id = len(g.nl)
    g.nl.append(ds.nl())
    g.nl[nl_id]['nl_id'] = nl_id
    g.nl[nl_id]['config_id'] = None
    g.nl[nl_id]['count'] = g.nl[src_id]['count']
    g.nl[nl_id]['labels'] = numpy.copy(g.nl[src_id]['labels'])
    g.nl[nl_id]['r'] = numpy.copy(g.nl[src_id]['r'])
    g.nl[nl_id]['rvec'] = numpy.copy(g.nl[src_id]['rvec'])
    g.nl[nl_id]['inhalo'] = numpy.copy(g.nl[src_id]['inhalo'])
    """



  """
  @staticmethod
  def copy(src_id, dest_id=None):
    if(dest_id == None):

  nl_id = len(g.nl)
  g.nl.append(ds.nl())
  """




