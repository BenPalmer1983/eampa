#!/bin/python3
"""
Output
"""

import os
import time
import numpy
import shutil

from datetime import datetime

from g import g
from label import label


class output:


  """
  Log
  ######################################################"""

  @staticmethod
  def log(inp, fh = None, verbose=0):
    if(not isinstance(inp, list) and "\n" in inp):
      inp = inp.split("\n")
    t = time.time() - g.start
    t = "{0:.4e}".format(t)
    while(len(t) < 10):
      t = " " + t
    if(verbose <= g.verbose['print']):
      if(isinstance(inp, list)):
        for inpn in inp:
          print(inpn)
      else:
        print(inp)
    if(isinstance(inp, list)):
      tofile = ''
      for inpn in inp:
        tofile = tofile + t + "   " + str(inpn) + '\n'
    else:
      tofile = t + "   " + str(inp) + '\n'

    
    if(verbose <= g.verbose['file']):
      g.fh.write(tofile)
      g.fh.flush()

    # Always save if a specific file provided
    if(fh != None):
      fh.write(tofile)


  """
  Configs
  ######################################################"""




  @staticmethod
  def config_summary(config_id):
    log = []

 
    atom_count = g.configs[config_id]['count']
    nl_id = g.configs[config_id]['nl_id']
    nl_count = g.nl[nl_id]['count']
    energy = g.configs[config_id]['c_energy']
    epa = g.configs[config_id]['c_energy'] / g.configs[config_id]['count']
    config_timer = g.configs[config_id]['timer']
    config_counter = g.configs[config_id]['counter']
    config_time = (config_timer / (1.0 * config_counter))
    calc_type = g.configs[config_id]['calc_type']
    if(calc_type is None):
      calc_type = ""

    log.append("{0:20s}{1:8d}".format("Config ID:", config_id))
    log.append("{0:20s}{1:8d}".format("NL ID:", nl_id))
    log.append("{0:20s}{1:8d}".format("Atoms:", atom_count))
    log.append("{0:20s}{1:8d}".format("NL Count:", nl_count))
    log.append("{0:20s}{1:8s}".format("Calc type:", calc_type))
    log.append("{0:20s}{1:12.6f}eV ({2:12.6f} eV per atom)".format("Energy:", energy, epa))
    log.append("{0:20s}{1:12.6f} {2:8d} {3:12.6f})".format("C Calc:", config_timer, config_counter, config_time))
    log.append("")

    # Output summary
    output.log(log, verbose=0)
    """
    for config_id in range(len(g.configs)):
      if(g.configs[config_id]['tag'] == 'efs'): 
        output.log("Config:   " + str(config_id), verbose=0)
        output.log("Count:   " + str(g.configs[config_id]['count']), verbose=0)
        output.log("Energy:   " + str(g.configs[config_id]['c_energy']) + "     " + str(g.configs[config_id]['c_energy'] / g.configs[config_id]['count']), verbose=0)
        output.log("", verbose=0)
    """




  @staticmethod
  def configs(dir, config_tag='all'):
    if(isinstance(config_tag, int)):
      output.config(dir, config_tag)
    else:
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
         output.config(dir, config_id)


  @staticmethod
  def config(dir, c_id):
    
    nl_id = g.configs[c_id]['nl_id']

    c_id_str = str(c_id)
    while(len(c_id_str) < 6):
      c_id_str = "0" + c_id_str
    fn = "config_" + c_id_str + ".out"
    path = os.path.join(dir, fn)

    d =     "CONFIG ID:    " + str(g.configs[c_id]['config_id']) + "\n"
    d = d + "NL ID:        " + str(g.configs[c_id]['nl_id']) + "\n"
    d = d + "TAG:          " + str(g.configs[c_id]['tag']) + "\n"
    d = d + "\n"
    d = d + "File:         " + str(g.configs[c_id]['file_path']) + "\n"
    d = d + "\n"
    d = d + "a0:           " + str(g.configs[c_id]['a0']) + "\n"
    d = d + "UV:           " + str(g.configs[c_id]['uv'][0,0]) + " " + str(g.configs[c_id]['uv'][0,1]) + " " + str(g.configs[c_id]['uv'][0,2]) +  "\n"
    d = d + "              " + str(g.configs[c_id]['uv'][1,0]) + " " + str(g.configs[c_id]['uv'][1,1]) + " " + str(g.configs[c_id]['uv'][1,2]) +  "\n"
    d = d + "              " + str(g.configs[c_id]['uv'][2,0]) + " " + str(g.configs[c_id]['uv'][2,1]) + " " + str(g.configs[c_id]['uv'][2,2]) +  "\n"
    d = d + "Volume:       " + str(g.configs[c_id]['volume']) + "\n"
    d = d + "Energy:       " + str(g.configs[c_id]['energy']) + "\n"
    if(g.configs[c_id]['stress'] is None):
      d = d + "Stress:       " + str(g.configs[c_id]['stress']) + "\n"
    else:
      d = d + "Stress:       " + str(g.configs[c_id]['stress'][0,0]) + " " + str(g.configs[c_id]['stress'][0,1]) + " " + str(g.configs[c_id]['stress'][0,2]) +  "\n"
      d = d + "              " + str(g.configs[c_id]['stress'][1,0]) + " " + str(g.configs[c_id]['stress'][1,1]) + " " + str(g.configs[c_id]['stress'][1,2]) +  "\n"
      d = d + "              " + str(g.configs[c_id]['stress'][2,0]) + " " + str(g.configs[c_id]['stress'][2,1]) + " " + str(g.configs[c_id]['stress'][2,2]) +  "\n"

    d = d + "\n"
    d = d + "Coords  (" + str(g.configs[c_id]['count']) + ")\n"
    d = d + output.lcfblock(g.configs[c_id]['labels'], g.configs[c_id]['coords'])
    d = d + "\n"

    d = d + "Forces  (" + str(g.configs[c_id]['count']) + ")\n"
    if(g.configs[c_id]['forces'] is None):
      d = d + str(g.configs[c_id]['forces'])  + "\n"
    else:
      d = d + output.datablock(g.configs[c_id]['forces'])
      d = d + "\n"
    
    uv = g.configs[c_id]['a0'] * g.configs[c_id]['uv'][:,:]
    cr = numpy.copy(g.configs[c_id]['coords'][:,:])

    for n in range(len(cr)):
      cr[n,:] = numpy.matmul(uv, cr[n,:])

    d = d + "\n"
    d = d + "Coords Real (" + str(g.configs[c_id]['count']) + ")\n"
    d = d + output.lcfblock(g.configs[c_id]['labels'], cr)
    d = d + "\n"



    if(nl_id is not None):
      d = d + "\n"
      d = d + "\n"
      d = d + "Calculated\n"
      d = d + "============================================================\n"
      d = d + "\n"
      d = d + "NL size:      " + str(g.nl[nl_id]['count']) + "\n"
      d = d + "rcut:         " + str(g.nl[nl_id]['count']) + "\n"
      d = d + "\n"
      d = d + "Energy:       " + str(g.configs[c_id]['c_energy']) + "\n" 
      d = d + "EPA:          " + str(g.configs[c_id]['c_energy'] / g.configs[c_id]['count']) + "\n" 
      d = d + "\n"
      if(g.configs[c_id]['c_stress_gpa'] is None):
        d = d + "Stress:       " + str(g.configs[c_id]['c_stress_gpa'])
      else:
        d = d + "Stress:       " + str(g.configs[c_id]['c_stress_gpa'][0,0]) + " " + str(g.configs[c_id]['c_stress_gpa'][0,1]) + " " + str(g.configs[c_id]['c_stress_gpa'][0,2]) +  "\n"
        d = d + "              " + str(g.configs[c_id]['c_stress_gpa'][1,0]) + " " + str(g.configs[c_id]['c_stress_gpa'][1,1]) + " " + str(g.configs[c_id]['c_stress_gpa'][1,2]) +  "\n"
        d = d + "              " + str(g.configs[c_id]['c_stress_gpa'][2,0]) + " " + str(g.configs[c_id]['c_stress_gpa'][2,1]) + " " + str(g.configs[c_id]['c_stress_gpa'][2,2]) +  "\n"
      d = d + "\n"
      d = d + "Forces\n"
      if(g.configs[c_id]['c_forces'] is None):
        d = d + str(g.configs[c_id]['c_forces'])  + "\n"
        d = d + "\n"
      else:
        d = d + output.datablock(g.configs[c_id]['c_forces'])
        d = d + "\n"     
      d = d + "RSS:         " + str(g.configs[c_id]['rss']) + "\n"
      d = d + "\n"

    fh = open(path, 'w')
    fh.write(d)
    fh.close()

    return True



  @staticmethod
  def config_forces(dir, config_tag='all'):
    if(isinstance(config_tag, int)):
      output.config_force(dir, config_tag)
    else:
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
         output.config_force(dir, config_id)

  @staticmethod
  def config_force(dir, c_id):
    forces = g.configs[c_id]['c_forces']
    if(forces is None):
      return False
    forces = forces.round(7)
    c_id_str = str(c_id)
    while(len(c_id_str) < 6):
      c_id_str = "0" + c_id_str
    fn = "config_" + c_id_str + ".forces"
    #fh = open(os.path.join(dir, fn), 'w')
    output.datafile(os.path.join(dir, fn), forces)
    return True



  @staticmethod
  def config_nls(dir, config_tag='all'):
    if(isinstance(config_tag, int)):
      output.config(dir, config_tag)
    else:
      for config_id in range(len(g.configs)):
        if(config_tag.lower() == 'all' or g.configs[config_id]['tag'] == config_tag.lower()):
         output.config_nl(dir, config_id)


  @staticmethod
  def config_nl(dir, config_id):
    nl_id = g.configs[config_id]['nl_id']
    c_str = str(config_id)
    while(len(c_str) < 6):
      c_str = "0" + c_str
    nl_str = str(nl_id)
    while(len(nl_str) < 6):
      nl_str = "0" + nl_str
    fname = c_str + "_" + nl_str + ".txt"
    nl = g.nl[nl_id]
    fh = open(os.path.join(dir, fname), 'w')
    for n in range(nl['count']):
      fh.write("{0:5d} {1:5d} {2:5d} {3:14.6f} \n".format(n, nl['labels'][n, 0], nl['labels'][n, 1], nl['r'][n]))
    fh.close()


  """
  Bulk Properties
  ######################################################"""


  @staticmethod
  def bp(dir):
    path = os.path.join(dir, "bp.txt")

    d = output.file_header('BP Calculations')

    for n in range(len(g.bp)):
       d = d + '********************************************************************************\n'
       d = d + '                              BP ' + str(n) + '\n'
       d = d + '********************************************************************************\n'
       d = d + '\n'
       d = d + 'Known:\n'
       d = d + '################################################################################\n'
       d = d + '\n'
       d = d + 'EoS:\n'
       d = d + 'a0:        ' + str(g.bp[n]['a0']) + '\n'
       d = d + 'uv:        ' + str(g.bp[n]['uv'][0,0]) + '  ' + str(g.bp[n]['uv'][0,1]) + '  ' + str(g.bp[n]['uv'][0,2]) + '  \n'
       d = d + '           ' + str(g.bp[n]['uv'][1,0]) + '  ' + str(g.bp[n]['uv'][1,1]) + '  ' + str(g.bp[n]['uv'][1,2]) + '  \n'
       d = d + '           ' + str(g.bp[n]['uv'][2,0]) + '  ' + str(g.bp[n]['uv'][2,1]) + '  ' + str(g.bp[n]['uv'][2,2]) + '  \n'
       d = d + 'v0:        ' + str(g.bp[n]['v0']) + '\n'
       d = d + 'e0:        ' + str(g.bp[n]['e0']) + '\n'
       d = d + 'b0:        ' + str(g.bp[n]['b0']) + '  (' + str(g.bp[n]['b0_gpa']) + ')\n'
       d = d + '\n'
       d = d + 'Elastic Constants\n'
       d = d + output.ecblock(g.bp[n]['ec_gpa'], 'ec (gpa):  ', '           ')
       d = d + '\n'
       d = d + 'Computed:\n'
       d = d + '################################################################################\n'
       d = d + '\n'
       d = d + 'EoS:\n'
       d = d + 'a0:        ' + str(g.bp[n]['c_a0']) + '\n'
       d = d + 'uv:        ' + str(g.bp[n]['c_uv'][0,0]) + '  ' + str(g.bp[n]['c_uv'][0,1]) + '  ' + str(g.bp[n]['c_uv'][0,2]) + '  \n'
       d = d + '           ' + str(g.bp[n]['c_uv'][1,0]) + '  ' + str(g.bp[n]['c_uv'][1,1]) + '  ' + str(g.bp[n]['c_uv'][1,2]) + '  \n'
       d = d + '           ' + str(g.bp[n]['c_uv'][2,0]) + '  ' + str(g.bp[n]['c_uv'][2,1]) + '  ' + str(g.bp[n]['c_uv'][2,2]) + '  \n'
       d = d + 'v0:        ' + str(g.bp[n]['c_v0']) + '\n'
       d = d + 'e0:        ' + str(g.bp[n]['c_e0']) + '\n'
       d = d + 'b0:        ' + str(g.bp[n]['c_b0']) + '  (' + str(g.bp[n]['c_b0_gpa']) + ')\n'
       d = d + '\n'
       d = d + 'Elastic Constants\n'
       d = d + 'Orthorhombic: Ravindran, Fast, Korzhavyi, Johansson 1998\n'
       d = d + output.ecblock(g.bp[n]['c_ec_gpa'], 'ec (gpa):  ', '           ')
       d = d + '\n'
       d = d + 'Elastic Constants\n'
       d = d + 'Cubic: Mehl, Singh, Klein, Papaconstantopoulos\n'
       d = d + output.ecblock(g.bp[n]['c_msp_ec_gpa'], 'ec (gpa):  ', '           ')
       d = d + '\n' 
       d = d + '\n'
       d = d + 'Other Computed Values:\n'
       d = d + '\n'
       d = d + 'b0 (r):     ' + str(g.bp[n]['c_b0_r']) + '  (' + str(g.bp[n]['c_b0_r_gpa']) + ')\n'
       d = d + 'b0 (v):     ' + str(g.bp[n]['c_b0_v']) + '  (' + str(g.bp[n]['c_b0_v_gpa']) + ')\n'
       d = d + 'b0 (avg):   ' + str(g.bp[n]['c_b0_avg']) + '  (' + str(g.bp[n]['c_b0_avg_gpa']) + ')\n'
       d = d + '\n'
       d = d + 'T melt:     ' + str(g.bp[n]['c_melting_temperature']) + '\n'
       d = d + '\n'
       d = d + '\n'
       d = d + 'RSS:\n'
       d = d + '################################################################################\n'
       for k in g.bp[n]['rss_details'].keys():
         d = d + "{0:24s}{1:22.8f}\n".format(k, g.bp[n]['rss_details'][k])   
       d = d + '\n'
       d = d + 'Total:        ' + str(g.bp[n]['rss']) + '\n'     
       d = d + '\n'
       d = d + '\n'
       d = d + '\n'
       d = d + '\n'



    fh = open(path, 'w')
    fh.write(d)
    fh.close()
    
    return d




  """
  Potentials
  ######################################################"""


  @staticmethod
  def potentials(dir=None):
    if(dir is None):
      dir = g.dir['data_potential']
    bn = os.path.basename(g.potential['index_file'])
    shutil.copyfile(g.potential['index_file'], os.path.join(dir, bn))


    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      file_name = os.path.basename(pfn['file'])
      file_name = file_name.split(".")
      file_name = str(file_name[0])

      if(pfn['tab'] is not None and pfn['fname'] == 'tab'):
        pot_name = file_name + ".tab"
        pot_path = os.path.join(dir, pot_name)
        output.potential_tabulated(pot_path, pfn['tab'])      
      else:
        pot_name = file_name + ".pot"
        pot_path = os.path.join(dir, pot_name)
        output.potential_analytic(pot_path, pfn)


  @staticmethod
  def potentials_tabulated(dir=None):
    if(dir is None):
      dir = g.dir['data_potential']
    bn = os.path.basename(g.potential['index_file'])
    shutil.copyfile(g.potential['index_file'], os.path.join(dir, bn))
    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      tab = g.potential['tab_for_output'][fn]
      file_name = os.path.basename(pfn['file'])
      #file_name = file_name.split(".")
      #file_name = str(file_name[0])
      #pot_name = file_name + ".tab"
      pot_path = os.path.join(dir, file_name)
      output.potential_tabulated(pot_path, tab) 


  @staticmethod
  def potentials_transformed(dir=None):
    if(dir is None):
      dir = g.dir['data_potential']
    bn = os.path.basename(g.potential['index_file'])
    shutil.copyfile(g.potential['index_file'], os.path.join(dir, bn))
    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      tab = g.potential['tab_transformed'][fn]
      file_name = os.path.basename(pfn['file'])
      pot_path = os.path.join(dir, file_name)
      output.potential_tabulated(pot_path, tab) 




  """
  Potfit Summary
  ######################################################"""

  @staticmethod
  def potfit_summary(dir=None):
    if(dir is None):
      dir = g.dir['fit']
    out =       "#########################################\n"  
    out = out + "                Summary                  \n"  
    out = out + "#########################################\n"    
    out = out + "\n"     

    for sn in range(len(g.potfit_steps)):
      t = g.potfit_steps[sn]['stats_time']
      out = out + "{:22s} {:12d}\n".format("Loop:", sn) 
      out = out + "{:22s} {:12s}\n".format("Minimiser:", g.potfit_steps[sn]['type'])  
      if(g.potfit_steps[sn]['stats_complete']):
        out = out + "{:22s} {:12s}\n".format("Complete:", "True")  
      else:
        out = out + "{:22s} {:12s}\n".format("Complete:", "False")  
      out = out + "{:22s} {:12.4f}\n".format("Time:", t)  
      try:
        out = out + "{:22s} {:12.4e}\n".format("Starting rss:", g.potfit_steps[sn]['stats_rss'][0])  
      except:
        pass
      try:
        out = out + "{:22s} {:12.4e}\n".format("Ending rss:", g.potfit_steps[sn]['stats_rss'][1])  
      except:
        pass
      out = out + "{:22s} {:12d} {:12.4f}\n".format("Calculation count:", g.potfit_steps[sn]['counter'], (g.potfit_steps[sn]['counter']/t)) 
      out = out + "{:22s} {:12d} {:12.4f}\n".format("Config count:", g.potfit_steps[sn]['stats_counter'], (g.potfit_steps[sn]['stats_counter']/t)) 
      out = out + "\n"    

    t = g.potfit['potfit_end_time'] - g.potfit['potfit_start_time']
    out = out + "\n"    
    out = out + "Stats\n"
    out = out + "\n"    
    out = out + "{:12s} {:16.2f}\n".format("Time: ", t)
    out = out + "{:12s} {:16d} {:16.2f}\n".format("Steps: ", g.potfit['counter'], g.potfit['counter']/ t)
    out = out + "{:12s} {:16d} {:16.2f}\n".format("Configs: ", g.calc_counter, g.calc_counter/ t)
    out = out + "{:12s} {:16.2f}\n".format("Best RSS: ", g.potfit['rss_best'])
    out = out + "\n"    
    

    fhpr = open(os.path.join(dir, 'potfit_results.txt'), 'w')
    output.log(out, fh = fhpr, verbose=1)
    fhpr.close()











  """
  Useful Functions
  ######################################################"""


  @staticmethod
  def file_header(title):
    now = datetime.now()
    out = '############################################################\n'
    for i in range((60 - len(title))//2):
      out = out + " "
    out = out + title + "\n"
    out = out + now.strftime("%H:%M:%S  %d/%m/%Y") + "\n"
    out = out + '############################################################\n\n'
    return out


  @staticmethod
  def potential_analytic(path, pfn):
    fh = open(path, 'w')
    
    fh.write("#TYPE " + str(pfn['fname']) + "\n")
    fh.write("#P &\n")
    for n in range(len(pfn['p'])):
      fh.write(str(pfn['p'][n]))
      if(n > 0 and (n+1) % 5 == 0 and n < (len(pfn['p']) - 1)):
        fh.write("&\n")
      else:
        fh.write("  ")
    fh.write("\n")
    fh.write("#PF &\n")
    for n in range(len(pfn['pf'])):
      fh.write(str(pfn['pf'][n]))
      if(n > 0 and (n+1) % 5 == 0 and n < (len(pfn['pf']) - 1)):
        fh.write("&\n")
      else:
        fh.write("  ")
    fh.write("\n")
    fh.close()


  @staticmethod
  def potential_tabulated(path, tab_in):
    fh = open(path, 'w')
    fh.write("#TYPE tab\n")
    for n in range(len(tab_in)):
      fh.write("{:16.8e}   {:16.8e}\n".format(tab_in[n, 0], tab_in[n, 1]))
    fh.close()


  @staticmethod
  def potential_data(path, tab_in):
    fh = open(path, 'w')
    for n in range(len(tab_in)):
      fh.write("{:16.8e}   {:16.8e}   {:16.8e}\n".format(tab_in[n, 0], tab_in[n, 1], tab_in[n, 2]))
    fh.close()


  @staticmethod
  def datafile(path, data_in, delim=None):
    fh = open(path, 'w')
    d = output.datablock(data_in, delim)
    fh.write(d)
    fh.close()



  @staticmethod
  def datablock(data_in, delim=None):
    d = ''
    for n in range(len(data_in)):
      for m in range(len(data_in[0,:])):      
        if(delim is None):
          v = "{0:.7e}".format(data_in[n, m])
          if(v[0] != "-"):
            v = " " + v
          d = d + " " + v + " "
        else:
          v = "{0:.7e}".format(data_in[n, m])
          if(m < len(data_in[0,:])-1):
            v = v + delim
          d = d + v     
      d = d + "\n"
    return d


  @staticmethod
  def lcfblock(labels, coords, forces=None, delim=None):
    d = ''
    for n in range(len(labels)):
      d = d + output.padr(label.get(labels[n]) + " [" + str(labels[n]) + "] ", 16) + " " 
      for m in range(len(coords[0,:])):   
        if(delim is None):
          v = "{0:.7e}".format(coords[n, m])
          if(v[0] != "-"):
            v = " " + v
          d = d + " " + v + " "
        else:
          v = "{0:.7e}".format(coords[n, m])
          if(m < len(coords[0,:])-1):
            v = v + delim
          d = d + v     
        if(forces is not None):
          for m in range(len(forces[0,:])):   
            if(delim is None):
              v = "{0:.7e}".format(forces[n, m])
              if(v[0] != "-"):
                v = " " + v
              d = d + " " + v + " "
            else:
              v = "{0:.7e}".format(forces[n, m])
              if(m < len(forces[0,:])-1):
                v = v + delim
              d = d + v     
      d = d + "\n"
    return d

  @staticmethod
  def ecblock(ec, indent_0='', indent_rest=''):
    out = ''
    for i in range(6):
      if(i == 0):
        out = out + indent_0
      else:
        out = out + indent_rest
      for j in range(6):
        #if(ec[i, j] == 0.0):
        #  out = out + "     0.0    "
        #else:
        out = out + "{0:10.2f}".format(ec[i, j]) + "  "
      out = out + "\n"
    return out
    
  @staticmethod
  def padl(inp, plen=17):
    inp = str(inp)
    while(len(inp) < plen):
       inp = " " + inp
    return inp
    
  @staticmethod
  def padr(inp, plen=17):
    inp = str(inp)
    while(len(inp) < plen):
       inp = inp + " "
    return inp
