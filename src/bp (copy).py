#!/bin/python3

import numpy
import time
import os

from f_toolbox import atom
from f_toolbox import math
from f_toolbox import transforms

from g import g
from ds import ds
from configs import configs
from nl import nl
from output import output
from plot import plot
from relax import relax
from eos import eos
from units import units


class bp: 

  log = []


  @staticmethod
  def run():
    output.log("Run BP Calculation", verbose=0)

    bp.make_dirs()
    bp.build()
     
    configs.index()
    exit() 
    bp.calc()
    bp.rss()
    bp.output()

    output.log("End", verbose=0)
    exit() 

  @staticmethod
  def output():
    bp.make_dirs()
    d = output.bp(g.dir['bp'])
    plot.plot_eos(g.dir['bp_plots'])
    plot.plot_ec(g.dir['bp_plots'])
    output.log(d, verbose=0)


  # Make additional directories for bp
  @staticmethod
  def make_dirs():    

    g.dir['bp_plots'] = os.path.join(g.dir['bp'], "plots")

    for k in g.dir.keys():
      os.makedirs(g.dir[k], exist_ok=True)


  """
  Calculate RSS
  ######################################################"""

  @staticmethod
  def rss():  
    total_rss = 0.0
    output.log("Calculate BP RSS", verbose=3)
    for n in range(len(g.bp)):
      rss_v = bp.rss_inner(n)
      output.log("BP " + str(n) + "  " + str(rss_v), verbose=2)
      total_rss = total_rss + rss_v
    output.log("BP total rss " + str(total_rss), verbose=2)
    return total_rss


  @staticmethod
  def rss_inner(n):    
    d = g.bp[n]
    rss_total = 0.0
    if(d['a0'] is not None and d['c_a0'] is not None):
      rss_total = rss_total + g.rss_w['a0'] * (d['a0'] - d['c_a0'])**2
    if(d['e0'] is not None and d['c_e0'] is not None):
      rss_total = rss_total + g.rss_w['e0'] * (d['e0'] - d['c_e0'])**2
    if(d['b0'] is not None and d['c_b0'] is not None):
      rss_total = rss_total + g.rss_w['b0'] * (d['b0'] - d['c_b0'])**2
    if(d['ec'] is not None and d['c_ec'] is not None):
      rss_total = rss_total + g.rss_w['ec'] * sum(sum((d['ec'] - d['c_ec'])**2))
    g.bp[n]['rss'] = g.rss_w['bp'] * rss_total
    return g.bp[n]['rss']



  """
  Build configurations
  ######################################################"""

  @staticmethod
  def build():
    output.log("Create BP configurations", verbose=3)
    for n in range(len(g.bp)):
      bp.build_inner(n)   

  @staticmethod
  def build_inner(n):
    if(g.bp[n]['built']):
      return False
    structure = g.bp[n]['structure']
    label = g.bp[n]['label']
    a0 = g.bp[n]['a0']
    uv = g.bp[n]['uv']
    rcut = g.bp_settings['rcut']
    cell_size = g.bp_settings['cell_size']

    # Create original unedited structure
    original_id = configs.add_common(structure, label, a0, uv, [cell_size], rcut, 'bp', 'e')
    g.bp[n]['original_id'] = original_id
    # Build NL
    nl.build_nl(g.bp[n]['original_id'])
    # Create Relax
    #g.bp[n]['relax_id'] = configs.add_common(structure, label, a0, uv, [cell_size], rcut, 'bp', 'e')
    g.bp[n]['relax_id'] = configs.duplicate(original_id)
 
    configs.index()
    exit()
    eos_strain = g.bp_settings['eos_strain']
    eos_steps = g.bp_settings['eos_steps']
    ec_strain = g.bp_settings['ec_strain']
    ec_steps = g.bp_settings['ec_steps']

    # Array sizes
    eoss = 2 * eos_steps + 1
    ecs = ec_steps

    # Create EOS
    g.bp[n]['eos_ids'] = []
    g.bp[n]['eos_strains'] = ds.bp_eos_arr(eoss)
    g.bp[n]['eos_volumes'] = ds.bp_eos_arr(eoss)
    g.bp[n]['eos_energies'] = ds.bp_eos_arr(eoss)

    g.bp[n]['ec_ids'] = []
    g.bp[n]['ec_strains'] = ds.bp_ec_arr(ecs)
    g.bp[n]['ec_volumes'] = ds.bp_ec_arr(ecs)
    g.bp[n]['ec_energies'] = ds.bp_ec_arr(ecs)

    g.bp[n]['c_ec'] = ds.bp_ec()
    g.bp[n]['c_ec_gpa'] = ds.bp_ec()

    k = 0
    for i in range(-eos_steps, eos_steps+1):
      s = eos_strain * (i / eos_steps)
      g.bp[n]['eos_strains'][k] = s
      g.bp[n]['eos_ids'].append(configs.add_common(structure, label, a0, uv, [cell_size], rcut, 'bp', 'e'))
      k = k + 1

    k = 0
    for di in range(9):
      g.bp[n]['ec_ids'].append([])
      for i in range(ec_steps):
        s = ec_strain * (i / (ec_steps - 1))
        g.bp[n]['ec_strains'][di, i] = s
        g.bp[n]['ec_ids'][di].append(configs.add_common(structure, label, a0, uv, [cell_size], rcut, 'bp', 'e'))
        k = k + 1
   
    return True



  """
  Compute EoS and Elastic Constants
  ######################################################"""

  @staticmethod
  def calc():

    output.log("Calculate BP", verbose=3)
    for n in range(len(g.bp)):
      bp.calc_inner(n) 


  @staticmethod
  def calc_inner(n):
    output.log("BP " + str(n), verbose=3)


    # Relax
    #################################

    relax_id = g.bp[n]['relax_id'] 

    
    # Relax
    relax.run_relax(relax_id, pmin='e', pvar=g.bp_settings['relax'])

    a0 = g.configs[relax_id]['a0']
    uv = g.configs[relax_id]['uv']

    g.bp[n]['c_a0'] = a0
    g.bp[n]['c_uv'] = uv


    v0 = math.cellvolume(a0, uv) / (1.0 * g.configs[relax_id]['count'])
    g.bp[n]['volrelaxed'] =  v0


    # Compute EoS
    #################################
    
    for eos_n in range(len(g.bp[n]['eos_ids'])):
      # Get id
      config_id = g.bp[n]['eos_ids'][eos_n]

      # Make transformation
      s = g.bp[n]['eos_strains'][eos_n]
      uvt = transforms.orthorhombic(uv, s)

      # Change cell and calculate energy
      nl.change_cell(config_id, a0=a0, uv=uvt, scale_uv=False)
      configs.calc_efs(config_id)


      # Save data
      g.bp[n]['eos_volumes'][eos_n] = math.cellvolume(a0, uvt) / (1.0 * g.configs[config_id]['count'])
      g.bp[n]['eos_energies'][eos_n] = g.configs[config_id]['c_energy'] / (1.0 * g.configs[config_id]['count'])

    eos_p = eos.fit(g.bp[n]['eos_volumes'], g.bp[n]['eos_energies'])
    g.bp[n]['c_v0'] = eos_p['v0']
    g.bp[n]['c_e0'] = eos_p['e0']
    g.bp[n]['c_b0'] = eos_p['b0']
    g.bp[n]['c_b0_gpa'] = eos_p['b0_gpa']
    g.bp[n]['c_b0p'] = eos_p['b0p']


    # Compute EC
    #################################
    for di in range(len(g.bp[n]['ec_ids'])):
      for i in range(len(g.bp[n]['ec_ids'][0])):
        # Get id
        config_id = g.bp[n]['ec_ids'][di][i]

        # Make transformation
        s = g.bp[n]['ec_strains'][di, i]
        uvt = transforms.dn(uv, s, di + 1)

        # Change cell and calculate energy
        nl.change_cell(config_id, a0=a0, uv=uvt, scale_uv=False)
        configs.calc_efs(config_id)


        # Save data
        g.bp[n]['ec_volumes'][di, i] = a0**3 * math.tripleproduct(uvt[:,0], uvt[:,1], uvt[:,2]) / g.configs[config_id]['count']
        g.bp[n]['ec_energies'][di, i] = g.configs[config_id]['c_energy'] / g.configs[config_id]['count']

    # Distortion data
    d = ds.bp_d()
    for di in range(len(g.bp[n]['ec_ids'])):
      p = numpy.polyfit(g.bp[n]['ec_strains'][di, :], g.bp[n]['ec_energies'][di, :], 2)
      d[di] = p[0] 

    # Calculate elastic constants
    g.bp[n]['c_ec'][0,0] = (2 / v0) * d[0]
    g.bp[n]['c_ec'][1,1] = (2 / v0) * d[1]
    g.bp[n]['c_ec'][2,2] = (2 / v0) * d[2]
    g.bp[n]['c_ec'][3,3] = (1 / (2 * v0)) * d[3]
    g.bp[n]['c_ec'][4,4] = (1 / (2 * v0)) * d[4]
    g.bp[n]['c_ec'][5,5] = (1 / (2 * v0)) * d[5]
    g.bp[n]['c_ec'][0,1] = (g.bp[n]['c_ec'][0,0] + g.bp[n]['c_ec'][1,1]) / 2 - (d[6] / v0)
    g.bp[n]['c_ec'][0,2] = (g.bp[n]['c_ec'][0,0] + g.bp[n]['c_ec'][2,2]) / 2 - (d[7] / v0)
    g.bp[n]['c_ec'][1,2] = (g.bp[n]['c_ec'][1,1] + g.bp[n]['c_ec'][2,2]) / 2 - (d[8] / v0)
    g.bp[n]['c_ec'][1,0] = g.bp[n]['c_ec'][0,1]
    g.bp[n]['c_ec'][2,0] = g.bp[n]['c_ec'][0,2]
    g.bp[n]['c_ec'][2,1] = g.bp[n]['c_ec'][1,2]

    # Calculate in GPA
    for i in range(6):
      for j in range(6):
        g.bp[n]['c_ec_gpa'][i,j] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_ec'][i,j])

    # Compliance tensor
    g.bp[n]['c_sc'] = math.minverse(g.bp[n]['c_ec'])

    # Other properties
    g.bp[n]['c_b0_r'] = bp.bulk_r(g.bp[n]['c_sc'])
    g.bp[n]['c_b0_v'] = bp.bulk_v(g.bp[n]['c_ec'])
    g.bp[n]['c_b0_avg'] = 0.5 * (g.bp[n]['c_b0_r'] + g.bp[n]['c_b0_v'])
    g.bp[n]['c_g_r'] = bp.shear_g_r(g.bp[n]['c_sc'])
    g.bp[n]['c_g_v'] = bp.shear_g_v(g.bp[n]['c_ec'])
    g.bp[n]['c_g_avg'] = 0.5 * (g.bp[n]['c_g_r'] + g.bp[n]['c_g_v'])
    g.bp[n]['c_e'] = bp.youngs_e(g.bp[n]['c_b0_avg'], g.bp[n]['c_g_avg'])
    g.bp[n]['c_melting_temperature'] = bp.melting_temperature(g.bp[n]['c_ec'])


    # GPA
    g.bp[n]['c_b0_r_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_b0_r'])
    g.bp[n]['c_b0_v_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_b0_v'])
    g.bp[n]['c_b0_avg_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_b0_avg'])
    g.bp[n]['c_g_r_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_g_r'])
    g.bp[n]['c_g_v_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_g_v'])
    g.bp[n]['c_g_avg_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_g_avg'])
    g.bp[n]['c_e_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_e'])




    output.log("", verbose=3)
    output.log("BP Results   " + str(n), verbose=3)
    output.log("===========================================", verbose=3)
    output.log("a0      " + str(g.bp[n]['c_a0']), verbose=3)
    output.log("uv      " + str(g.bp[n]['c_uv'][0,0]) + " " + str(g.bp[n]['c_uv'][0,1]) + " " + str(g.bp[n]['c_uv'][0,2]), verbose=3)
    output.log("        " + str(g.bp[n]['c_uv'][1,0]) + " " + str(g.bp[n]['c_uv'][1,1]) + " " + str(g.bp[n]['c_uv'][1,2]), verbose=3)
    output.log("        " + str(g.bp[n]['c_uv'][2,0]) + " " + str(g.bp[n]['c_uv'][2,1]) + " " + str(g.bp[n]['c_uv'][2,2]), verbose=3)
    output.log("v0      " + str(g.bp[n]['c_v0']), verbose=3)
    output.log("e0      " + str(g.bp[n]['c_e0']), verbose=3)
    output.log("b0      " + str(g.bp[n]['c_b0']), verbose=3)
    output.log("b0gpa   " + str(g.bp[n]['c_b0_gpa']), verbose=3)
    output.log("", verbose=3)
    output.log("", verbose=3)


  @staticmethod
  def bulk_r(sc):
    return 1.0 / (sc[0,0] + sc[1,1] + sc[2,2] + 2.0 * (sc[0,1] + sc[0,2] + sc[1,2]))

  @staticmethod
  def bulk_v(ec):
    return ((1.0 / 9.0) * (ec[0,0] + ec[1,1] + ec[2,2]) + (2.0 / 9.0) * (ec[0,1] + ec[0,2] + ec[1,2]))


  @staticmethod
  def shear_g_r(sc):
    return (15.0 / (4*(sc[0,0]+sc[1,1]+sc[2,2]) - 4*(sc[0,1]+sc[0,2]+sc[1,2]) + 3* (sc[3,3]+sc[4,4]+sc[5,5])))

  @staticmethod
  def shear_g_v(ec):
    return (1.0/15.0) * (ec[0,0]+ec[1,1]+ec[2,2]-ec[0,1]-ec[0,2]-ec[1,2]) + (1.0/5.0) * (ec[3,3] + ec[4,4] + ec[5,5])

  @staticmethod
  def shear_g_vec(sc):
    gvec = numpy.zeros((3,), dtype=numpy.float64)
    gvec[0] = 1.0 / (2.0 * sc[3,3])
    gvec[1] = 1.0 / (2.0 * sc[4,4])
    gvec[2] = 1.0 / (2.0 * sc[5,5])
    return gvec

  @staticmethod
  def youngs_e(b0_avg, g_avg):
    return (9 * b0_avg * g_avg) / (3 * b0_avg + g_avg)

  @staticmethod
  def shear_e_vec(sc):
    evec = numpy.zeros((3,), dtype=numpy.float64)
    evec[0] = 1.0 / sc[0,0]
    evec[1] = 1.0 / sc[1,1]
    evec[2] = 1.0 / sc[2,2]
    return evec

  @staticmethod
  def melting_temperature(ec): 
    return 598.0 + 6.66 * (ec[0,0]+ec[1,1]+ec[2,2]) - 0.003 * (ec[0,0]+ec[1,1]+ec[2,2])**2




"""


calc_v(bp_id) = (3.0D0 * calc_b0_avg(bp_id) * calc_g_avg(bp_id)) / ( 3 * calc_b0_avg(bp_id) + calc_g_avg(bp_id))


calc_v_vec(bp_id, 1, 1) = 1.0D0
calc_v_vec(bp_id, 2, 1) = -(calc_sc(bp_id, 1, 2) * calc_e_vec(bp_id, 2))
calc_v_vec(bp_id, 3, 1) = -(calc_sc(bp_id, 1, 3) * calc_e_vec(bp_id, 3))
calc_v_vec(bp_id, 1, 2) = -(calc_sc(bp_id, 2, 1) * calc_e_vec(bp_id, 1))
calc_v_vec(bp_id, 2, 2) = 1.0D0
calc_v_vec(bp_id, 3, 2) = -(calc_sc(bp_id, 2, 3) * calc_e_vec(bp_id, 3))
calc_v_vec(bp_id, 1, 3) = -(calc_sc(bp_id, 3, 1) * calc_e_vec(bp_id, 1))
calc_v_vec(bp_id, 2, 3) = -(calc_sc(bp_id, 3, 2) * calc_e_vec(bp_id, 2))
calc_v_vec(bp_id, 3, 3) = 1.0D0

b0 = 160.230732254D0 * calc_b0_avg(bp_id)
g = 160.230732254D0 * calc_g_avg(bp_id)


calc_vl(bp_id) = sqrt((b0 + (4.0D0/3.0D0) * g) / calc_rho(bp_id))
calc_vt(bp_id) = sqrt(g / calc_rho(bp_id))
calc_vm(bp_id) = ((1.0D0/3.0D0) * (2.0D0 / (calc_vt(bp_id)**3)) + 1.0D0 / (calc_vl(bp_id)**3))**(-(1.0D0/3.0D0))

calc_melting(bp_id) =   598.0D0 &
                      + 6.66D0 * (calc_ec_gpa(bp_id, 1, 1) + calc_ec_gpa(bp_id, 2, 2) + calc_ec_gpa(bp_id, 3, 3)) &
                      - 0.003D0 * (calc_ec_gpa(bp_id, 1, 1) + calc_ec_gpa(bp_id, 2, 2) + calc_ec_gpa(bp_id, 3, 3))**2  
       
h = 6.62607004D-34
k = 1.38064852D-23
pi = 3.1415926535898D0
na = 6.02214076D23
calc_debye(bp_id) = ((h / k) * &
                    ((( 3.0D0 * known_atoms_per_crystal(bp_id)) / (4.0D0 * pi)) * &
                    ((na * calc_rho(bp_id)) / (known_amu_per_crystal(bp_id) * 1.66D-27))**(1.0D0/3.0D0)) &
                    * calc_vm(bp_id))


calc_cubic_cauchy_pressure(bp_id) = calc_ec(bp_id, 1, 2) - calc_ec(bp_id, 4, 4)
calc_cubic_tetragonal_shear(bp_id) = (calc_ec(bp_id, 1, 1) - calc_ec(bp_id, 1, 2)) / 2
calc_cubic_shear_modulus(bp_id) = calc_ec(bp_id, 4, 4) 

calc_cubic_cauchy_pressure_gpa(bp_id) = calc_ec_gpa(bp_id, 1, 2) - calc_ec_gpa(bp_id, 4, 4)
calc_cubic_tetragonal_shear_gpa(bp_id) = (calc_ec_gpa(bp_id, 1, 1) - calc_ec_gpa(bp_id, 1, 2)) / 2
calc_cubic_shear_modulus_gpa(bp_id) = calc_ec_gpa(bp_id, 4, 4) 


IF((calc_ec(bp_id, 1, 1) + 2 * calc_ec(bp_id, 1, 2)) > 0 .AND. &
   (calc_ec(bp_id, 1, 1) - calc_ec(bp_id, 1, 2)) > 0 .AND. &
   (calc_ec(bp_id, 4, 4)) > 0)THEN
  calc_cubic_stability(bp_id) = 1
ELSE
  calc_cubic_stability(bp_id) = 0
END IF


IF((calc_ec(bp_id, 1, 1) > 0) .AND. &
   (calc_ec(bp_id, 2, 2) > 0) .AND. &
   (calc_ec(bp_id, 3, 3) > 0) .AND. &
   (calc_ec(bp_id, 4, 4) > 0) .AND. &
   (calc_ec(bp_id, 5, 5) > 0) .AND. &
   (calc_ec(bp_id, 6, 6) > 0) .AND. &
   ((calc_ec(bp_id, 1, 1) + calc_ec(bp_id, 2, 2) - 2 * calc_ec(bp_id, 1, 2)) > 0) .AND. &
   ((calc_ec(bp_id, 1, 1) + calc_ec(bp_id, 3, 3) - 2 * calc_ec(bp_id, 1, 3)) > 0) .AND. &
   ((calc_ec(bp_id, 2, 2) + calc_ec(bp_id, 3, 3) - 2 * calc_ec(bp_id, 2, 3)) > 0) .AND. &
   ((calc_ec(bp_id, 1, 1) + calc_ec(bp_id, 2, 2) + calc_ec(bp_id, 3, 3) + &
   2 * calc_ec(bp_id, 1, 2) + 2 * calc_ec(bp_id, 1, 3) + 2 * calc_ec(bp_id, 2, 3))> 0))THEN
  calc_stability(bp_id) = 1
ELSE
  calc_stability(bp_id) = 0
END IF









!# SHEAR MODULUS
calc_g_r(bp_id) = 15.0D0 / (&
                    4.0D0 * (calc_sc(bp_id, 1, 1) + calc_sc(bp_id, 2, 2) + calc_sc(bp_id, 3, 3)) &
                  - 4.0D0 * (calc_sc(bp_id, 1, 2) + calc_sc(bp_id, 1, 3) + calc_sc(bp_id, 2, 3)) &
                  + 3.0D0 * (calc_sc(bp_id, 4, 4) + calc_sc(bp_id, 5, 5) + calc_sc(bp_id, 6, 6)) &
                  )
calc_g_v(bp_id) = (1.0D0 / 15.0D0) * (&
                  calc_ec(bp_id, 1, 1) + calc_ec(bp_id, 2, 2) + calc_ec(bp_id, 3, 3) &
                  - calc_ec(bp_id, 1, 2) - calc_ec(bp_id, 1, 3) - calc_ec(bp_id, 2, 3) &
                  ) &
                  + (1.0D0/5.0D0) * (&
                  calc_ec(bp_id, 4, 4) + calc_ec(bp_id, 5, 5) + calc_ec(bp_id, 6, 6) &
                  )
calc_g_avg(bp_id) = 0.5D0 * (calc_g_r(bp_id) + calc_g_v(bp_id))


calc_g_vec(bp_id, 1) = 1.0D0 / (2.0D0 * calc_sc(bp_id, 4, 4))
calc_g_vec(bp_id, 2) = 1.0D0 / (2.0D0 * calc_sc(bp_id, 5, 5))
calc_g_vec(bp_id, 3) = 1.0D0 / (2.0D0 * calc_sc(bp_id, 6, 6))




calc_e(bp_id) = (9 * calc_b0_avg(bp_id) * calc_g_avg(bp_id)) / (3.0D0 * calc_b0_avg(bp_id) + calc_g_avg(bp_id))

calc_e_vec(bp_id, 1) = 1.0D0 / calc_sc(bp_id, 1, 1)
calc_e_vec(bp_id, 2) = 1.0D0 / calc_sc(bp_id, 2, 2)
calc_e_vec(bp_id, 3) = 1.0D0 / calc_sc(bp_id, 3, 3)





calc_melting(bp_id) =   598.0D0 &
                      + 6.66D0 * (calc_ec_gpa(bp_id, 1, 1) + calc_ec_gpa(bp_id, 2, 2) + calc_ec_gpa(bp_id, 3, 3)) &
                      - 0.003D0 * (calc_ec_gpa(bp_id, 1, 1) + calc_ec_gpa(bp_id, 2, 2) + calc_ec_gpa(bp_id, 3, 3))**2  
















"""


