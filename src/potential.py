#!/bin/python3

import numpy
import time
import os
import matplotlib.pyplot as plt

from f_toolbox import fnc as fnc
from f_toolbox import interp as interp
from f_toolbox import atom as atom

from g import g
from ds import ds
from label import label
from std import std
from output import output
from plot import plot



class potential:


  




  @staticmethod
  def load():

    potindex  = g.input['potindex']

    if(not os.path.isfile(potindex)):
      output.log("Potential index file does not exist:", verbose=0)
      output.log("   " + str(potindex), verbose=0)
      output.log("Exiting.", verbose=0)
      exit()
    
    #   Read file
    ###############################################################

    d = std.read(potindex)

    #   Read data from index
    ###############################################################

    g.potential = ds.potential()
    potdir = os.path.dirname(potindex)
    g.potential['index_file'] = potindex
    pf_id = None
    
    for line in d:
      f = line.split(" ")
      f[0] = f[0].upper()
      if(f[0] == "#POTNAME"):
        pass

      elif(f[0] == "#FILE"):
        pf_id = len(g.potential['functions'])
        g.potential['functions'].append(ds.potential_function())
        g.potential['functions'][pf_id]['file'] = os.path.join(potdir, f[1])


      elif(f[0] == "#LABEL"):
        g.potential['functions'][pf_id]['a'] = label.set(f[1])
        if(len(f) > 2):
          g.potential['functions'][pf_id]['b'] = label.set(f[2])


      elif(f[0] == "#F_ON"):  
        if(f[1][0].upper() == 'T'):
          g.potential['functions'][pf_id]['fon'] = True
        else:
          g.potential['functions'][pf_id]['fon'] = False


      elif(f[0] == "#F_TYPE"): 
        g.potential['functions'][pf_id]['ftype'] = 1
        if(f[1][0].upper() == 'D'):
          g.potential['functions'][pf_id]['ftype'] = 2
        elif(f[1][0].upper() == 'E'):
          g.potential['functions'][pf_id]['ftype'] = 3


      elif(f[0] == "#F_GROUP"): 
        g.potential['functions'][pf_id]['fgroup'] = -1
        if(len(f) == 2):
          g.potential['functions'][pf_id]['fgroup'] = potential.group(f[1].upper())
          

      elif(f[0] == "#R_CUT"): 
        g.potential['functions'][pf_id]['rcut'] = None
        if(len(f) == 2):
          g.potential['functions'][pf_id]['rcut'] = float(f[1])


      elif(f[0] == "#FIX"):  
        if(f[1][0].upper() == 'T'):
          g.potential['functions'][pf_id]['fix'] = True
        else:
          g.potential['functions'][pf_id]['fix'] = False





    #   Read data from potential files
    ###############################################################    
    for n in range(len(g.potential['functions'])):
      fname = None

      # Read function type
      if(not os.path.isfile(g.potential['functions'][n]['file'])):
        output.log("Potential files does not exist:", verbose=0)
        output.log("   " + str(g.potential['functions'][n]['file']), verbose=0)
        output.log("Check potential index file and function files.", verbose=0)
        output.log("Exiting.", verbose=0)
        exit()
    
      fh = open(g.potential['functions'][n]['file'])
      for line in fh:
        line = std.one_space(line.strip())
        f = line.split(" ")
        if(f[0].upper() == "#TYPE"):
          fname = f[1]
          break
      fh.close()
      g.potential['functions'][n]['fname'] = fname.lower()

      # Load either tabulated or analytic
      if(g.potential['functions'][n]['fname'] == "tab"):
        g.potential['functions'][n]['tab'] = potential.read_tab(g.potential['functions'][n]['file'])
      else:
        g.potential['functions'][n]['p'], g.potential['functions'][n]['pf'] = potential.read_analytic(g.potential['functions'][n]['file'])


    #   Make tabulated versions where rcut is set
    ###############################################################   
    potential.make_tab()


    #   Output data and plot
    ###############################################################   

    # Plot potentials
    plot.potentials()
    output.potentials()


  @staticmethod
  def make_tab():
    output.log("Make tabulated where rcut has been set", verbose=2)
    for n in range(len(g.potential['functions'])):
      if(g.potential['functions'][n]['fname'] != 'tab' and g.potential['functions'][n]['rcut'] is not None):
        # Make array
        g.potential['functions'][n]['tab'] = ds.tab()
        # x-values
        g.potential['functions'][n]['tab'][:, 0] = numpy.linspace(0.0, g.potential['functions'][n]['rcut'], len(g.potential['functions'][n]['tab']))
        # fx-values
        g.potential['functions'][n]['tab'][:, 1] = fnc.fv(g.potential['functions'][n]['fname'], g.potential['functions'][n]['tab'][:, 0], g.potential['functions'][n]['p'], g.potential['functions'][n]['pf']) 
        # f'x values 5 point interp
        x = numpy.copy(g.potential['functions'][n]['tab'][:, 0])
        y = numpy.copy(g.potential['functions'][n]['tab'][:, 1])
        interp.fill(x[:], y[:], 5, g.potential['functions'][n]['tab'])
        # fx-values
        #g.potential['functions'][n]['tab'][:, 2] = fnc.fgradv(g.potential['functions'][n]['fname'], g.potential['functions'][n]['tab'][:, 0], g.potential['functions'][n]['p'], g.potential['functions'][n]['pf']) 
        
        output.log("...tab created for " + str(n) + " " + g.potential['functions'][n]['fname'], verbose=2)
      elif(g.potential['functions'][n]['fname'] != 'tab' and g.potential['functions'][n]['rcut'] is None):
        output.log("...tab not created for " + str(n) + " " + g.potential['functions'][n]['fname'] + " no rcut set", verbose=2)


  @staticmethod
  def tab_for_output(rcut_none_set=None):
    if(rcut_none_set is None):
      rcut_none_set = 10.0
    g.potential['tab_for_output'] = []
    for fn in range(len(g.potential['functions'])):
      g.potential['tab_for_output'].append(None)
    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      if(g.potential['functions'][fn]['fname'] == 'tab'):
        g.potential['tab_for_output'][fn] = pfn['tab']
      else:
        # Make tabulated
        tab = ds.tab(tlen=10001)
        # rcut
        rcut = g.potential['functions'][fn]['rcut']
        if(rcut is None):
          rcut = rcut_none_set
        # x-values
        tab[:, 0] = numpy.linspace(0.0, rcut, len(tab))
        # fx-values        
        tab[:, 1] = fnc.fv(g.potential['functions'][fn]['fname'], tab[:, 0], g.potential['functions'][fn]['p'], g.potential['functions'][fn]['pf']) 
        # fx'-values
        x = numpy.copy(tab[:, 0])
        y = numpy.copy(tab[:, 1])
        interp.fill(x[:], y[:], 5, tab)
        #tab[:, 2] = fnc.fgradv(g.potential['functions'][fn]['fname'], tab[:, 0], g.potential['functions'][fn]['p'], g.potential['functions'][fn]['pf']) 
        g.potential['tab_for_output'][fn] = tab


  @staticmethod
  def read_tab(file_path):
    d = []
    fh = open(file_path, 'r')
    for line in fh: 
      line = line.strip()
      if(line[0] != "#"):
        line = std.one_space(line)
        f = line.split(" ")
        if(len(f) > 1):
          d.append(f)
    fh.close()  
    d = numpy.asarray(d, dtype=numpy.float64, order='F')
    tab = ds.tab()
    interp.fill(d[:,0], d[:,1], 5, tab)
    return tab


  @staticmethod
  def read_analytic(file_path):
    p = None
    pf = None

    # Read file in
    raw = []
    fh = open(file_path, "r")
    for line in fh:
      raw.append(line.strip())
    fh.close()

    # Remove comments
    for n in range(len(raw)):
      t = raw[n].split("//")
      raw[n] = t[0]
      t = raw[n].split("!")
      raw[n] = t[0]

    temp = []
    # Read it in line by line
    cont = False
    for n in range(len(raw)):
      cl = False
      line = raw[n].strip()
      if(len(line) > 0 and line[-1] == "&"):
        cl = True 
        line = line[:-1].strip() + " "
      if(len(temp) > 0 and cont):
        temp[-1] = temp[-1] + line
        cont = False
      else:  
        temp.append(line)
      if(cl):
        cont = True

    file_data = []
    for line in temp:
      if(line.strip() != ""):
        line = line.replace(",", " ")
        line = line.replace("\t", " ")
        file_data.append(std.one_space(line).strip())

    # Default
    p = numpy.zeros((1,),)
    pf = numpy.zeros((1,),)

    # Read the data
    for line in file_data:
      f = line.split(" ")
      f[0] = f[0].upper()
      if(f[0] == "#P"):
        p = numpy.asarray(f[1:], dtype=numpy.float64)
      elif(f[0] == "#PF"):
        pf = numpy.asarray(f[1:], dtype=numpy.float64)
        
    return p, pf



  @staticmethod
  def group(inp):
    inp = inp.strip().upper()
    for n in range(len(g.groups)):
      if(g.groups[n] == inp):
        return n
    g.groups.append(inp)
    return (len(g.groups) - 1)




  # Parameter array for fitting
  ##############################################################

  @staticmethod
  def update_p(p_in):
    p = ds.potential_p(p_in)    

    a = 0
    for n in range(len(g.potential['functions'])):
      if(g.potential['functions'][n]['fon'] 
        and g.potential['functions'][n]['p'] is not None 
        and g.potential['functions'][n]['fix'] == False):
        b = a + len(g.potential['functions'][n]['p'])
        g.potential['functions'][n]['p'][:] = p[a:b]
        a = b

    # Remake tabulated where needed
    potential.make_tab()

    # Set potential in Fortran
    potential.set()



  #   Get parameter array
  ###############################################################   
  @staticmethod
  def get_p():
    p = []
    for n in range(len(g.potential['functions'])):
      if(g.potential['functions'][n]['fon'] 
        and g.potential['functions'][n]['p'] is not None 
        and g.potential['functions'][n]['fix'] == False):
        for m in range(len(g.potential['functions'][n]['p'])):
          p.append(g.potential['functions'][n]['p'][m]) 
    return ds.potential_p(p)



  # Set to Fortran
  ##############################################################

  @staticmethod
  def set():
    output.log("Set potentials", verbose=2)
    for n in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][n]
      if(pfn['tab'] is not None):
        if(pfn['fortran_pkey'] is None):
          ftype = pfn['ftype']
          a = pfn['a']
          b = pfn['b']
          fgroup = pfn['fgroup']
          g.potential['functions'][n]['fortran_pkey'] = atom.add_tabulated(ftype, a, b, fgroup, pfn['tab'])    
        else:
          atom.update_tabulated(pfn['tab'], g.potential['functions'][n]['fortran_pkey'])
      else:
        if(pfn['fortran_pkey'] is None):
          ftype = pfn['ftype']
          a = pfn['a']
          b = pfn['b']
          fgroup = pfn['fgroup']
          fname = pfn['fname']
          p = pfn['p']
          pf = pfn['pf']
          g.potential['functions'][n]['fortran_pkey'] = atom.add_analytic(ftype, a, b, fgroup, fname, p, pf)
        else:
          fname = pfn['fname']
          p = pfn['p']
          pf = pfn['pf']
          atom.update_analytic(fname, p, pf, g.potential['functions'][n]['fortran_pkey'])





  @staticmethod
  def summary_fortran():
    output.log("", verbose=2)
    output.log("Potential Summary - Fortran", verbose=2)
    output.log("###########################", verbose=2)
    if(atom.fcount == 0):    
      output.log("No functions loaded")
    else: 
      output.log("Function Count: " + str(atom.fcount), verbose=2)      
      output.log("PAIR", verbose=2)
      for a in range(atom.lgmax):
        for m in range(atom.lgmax):
          if(atom.pair_index_a[a, m] == -1):
            break
          else:
            b = atom.pair_index_a[a, m]
            output.log(label.get(a) + " (" + str(a) + ") " + label.get(b) + " (" + str(b) + ")    key: " + str(atom.pair_index_b[atom.pair_key(a, b)]), verbose=2)
      
      output.log("DENS", verbose=2)
      for a in range(atom.lgmax):
        for m in range(atom.lgmax):
          if(atom.dens_index_a[a, m] == -1):
            break
          else:
            fgroup = atom.dens_index_a[a, m]
            output.log(label.get(a) + " (" + str(a) + ") " + g.groups[fgroup] + " (" + str(fgroup) + ")    key: " + str(atom.dens_index_b[a, fgroup]), verbose=2)
      
      output.log("EMBE", verbose=2)
      for a in range(atom.lgmax):
        for m in range(atom.lgmax):
          if(atom.embe_index_a[a, m] == -1):
            break
          else:
            fgroup = atom.embe_index_a[a, m]
            output.log(label.get(a) + " (" + str(a) + ") " + g.groups[fgroup] + " (" + str(fgroup) + ")    key: " + str(atom.embe_index_b[a, fgroup]), verbose=2)
      output.log("", verbose=2)


