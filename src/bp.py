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
from potential import potential
from nl import nl
from output import output
from plot import plot
from relax import relax
from eos import eos
from units import units
from rose import rose
from msp import msp

class bp: 

  log = []


  @staticmethod
  def run():
    output.log("Run BP Calculation", verbose=0)

    bp.make_dirs()
    bp.build()     
    bp.calc()
    bp.rss()
    bp.output()
    output.log("End", verbose=0)
    exit() 



  @staticmethod
  def run_rss():
    bp.build()     
    bp.calc()
    return bp.rss()


  @staticmethod
  def output():
    bp.make_dirs()

    # Output summary
    d = output.bp(g.dir['bp'])
    output.log(d, verbose=0)

    # EoS, EC, Rose, MSP plots
    plot.plot_eos(g.dir['bp_plots'])
    plot.plot_ec(g.dir['bp_plots'])
    plot.plot_rose(g.dir['bp_plots'])
    plot.plot_msp(g.dir['bp_plots'])

    # Potential data and plots
    potential.tab_for_output()
    output.potentials(g.dir['bp_pot'])
    output.potentials_tabulated(g.dir['bp_pot_data'])
    plot.potentials(g.dir['bp_pot_plots'])


  # Make additional directories for bp
  @staticmethod
  def make_dirs():    

    g.dir['bp_plots'] = os.path.join(g.dir['bp'], "plots")
    g.dir['bp_pot'] = os.path.join(g.dir['bp'], "pot")
    g.dir['bp_pot_plots'] = os.path.join(g.dir['bp'], "pot_plots")
    g.dir['bp_pot_data'] = os.path.join(g.dir['bp'], "pot_data")

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
    g.bp[n]['rss_details'] = ds.bp_rss_details()
    d = g.bp[n]


    # RSS to specific values (EoS and 9 strains)
    if(d['a0'] is not None and d['c_a0'] is not None):
      g.bp[n]['rss_details']['a0'] = g.rss_w['a0'] * (d['a0'] - d['c_a0'])**2
    if(d['e0'] is not None and d['c_e0'] is not None):
      g.bp[n]['rss_details']['e0'] = g.rss_w['e0'] * (d['e0'] - d['c_e0'])**2
    if(d['b0'] is not None and d['c_b0'] is not None):
      g.bp[n]['rss_details']['b0'] = g.rss_w['b0'] * (d['b0'] - d['c_b0'])**2
    if(d['ec'] is not None and d['c_ec'] is not None):
      g.bp[n]['rss_details']['ec'] = g.rss_w['ec'] * sum(sum((d['ec'] - d['c_ec'])**2))

    # Match to Rose eos over a wider parameter range 
    if(d['rose_eos'] and d['rose_e'] is not None and d['c_rose_e'] is not None):
      g.bp[n]['rss_details']['rose'] = g.rss_w['rose'] * sum((d['rose_e'] - d['c_rose_e'])**2)

    # Match to MSP, MKP cubic elastic constants over strain
    if(d['msp_ec'] and d['msp_e'] is not None and d['c_msp_e'] is not None):
      msp_e_zeroed = d['msp_e'][:,:] - d['msp_e'][0,:]
      c_msp_e_zeroed = d['c_msp_e'][:,:] - d['c_msp_e'][0,:]
      g.bp[n]['rss_details']['msp_shape'] = g.rss_w['msp'][0] * sum(sum((msp_e_zeroed - c_msp_e_zeroed)**2))
      g.bp[n]['rss_details']['msp_values'] = g.rss_w['msp'][1] * sum(sum((d['msp_e'] - d['c_msp_e'])**2))
    
    # Fitting to BM EoS points
    if(d['bm_eos'] and d['bm_eos_e_known_z'] is not None and d['bm_eos_e_potential_z'] is not None):
      g.bp[n]['rss_details']['bm_eos_shape'] = g.rss_w['bm_eos'][0] * sum((d['bm_eos_e_known_z'] - d['bm_eos_e_potential_z'])**2) 
    if(d['bm_eos'] and d['bm_eos_e_known'] is not None and d['bm_eos_e_potential'] is not None):
      g.bp[n]['rss_details']['bm_eos_values'] = g.rss_w['bm_eos'][1] * sum((d['bm_eos_e_known'] - d['bm_eos_e_potential'])**2)


    
    # Penalties for negative values and instabilities
    penalty = 0.0 
    if(d['c_b0'] < 0.0):
      penalty = penalty + 1.0e3
    for i in range(6):
      for j in range(6):
        if(d['c_ec'][i,j] < 0.0):
          penalty = penalty + 1.0e2
    if(not d['c_cubic_stability']):
      penalty = penalty + 1.0e3
    g.bp[n]['rss_details']['bp_penalty'] = penalty    


    # Weight
    t = 0.0
    for k in g.bp[n]['rss_details'].keys():
      g.bp[n]['rss_details'][k] = g.bp[n]['weight'] * g.rss_w['bp'] * g.bp[n]['rss_details'][k]
      t = t + g.bp[n]['rss_details'][k]  
    g.bp[n]['rss_details']['bp_total'] = t
    g.bp[n]['rss'] = t
    #print(g.bp[n]['rss_details'])
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
      bp.reset()
      return False

    # Get details for structure
    structure = g.bp[n]['structure']
    label = g.bp[n]['label']
    a0 = g.bp[n]['a0']
    uv = g.bp[n]['uv']

    # Create original unedited structure and build neighbour list
    rcut = g.bp_settings['eos_rcut']
    cell_size = g.bp_settings['eos_cell_size']
    eos_bm_original_id = configs.add_common(structure, label, a0, uv, [cell_size], rcut, 'bp', 'e')
    g.bp[n]['eos_bm_original_id'] = eos_bm_original_id
    nl.build_nl(g.bp[n]['eos_bm_original_id'])


    """
    Configs to relax and find a0 for the potential
    ######################################################"""

    # Create Relax
    g.bp[n]['relax_id'] = configs.duplicate(eos_bm_original_id)


    """
    BM EoS configs to calculate values E0 V0 B0 for 
    this potential
    ######################################################"""

    eos_strain = g.bp_settings['eos_strain']
    eos_steps = g.bp_settings['eos_steps']
    eoss = 2 * eos_steps + 1

    g.bp[n]['eos_ids'] = []
    g.bp[n]['eos_strains'] = ds.bp_eos_arr(eoss)
    g.bp[n]['eos_volumes'] = ds.bp_eos_arr(eoss)
    g.bp[n]['eos_energies'] = ds.bp_eos_arr(eoss)

    k = 0
    for i in range(-eos_steps, eos_steps+1):
      s = eos_strain * (i / eos_steps)
      g.bp[n]['eos_strains'][k] = s
      g.bp[n]['eos_ids'].append(configs.duplicate(eos_bm_original_id))
      k = k + 1


    """
    ! FOR COMPUTE VALUE CALCULATIONS
    ! Ravindran, Fast, Korzhavyi, Johansson 1998
    ! Density functional theory for calculation of elastic properties of 
    ! orthorhombic crystals: application to TiSi2
    ! 
    ######################################################"""

    rcut = g.bp_settings['ec_rcut']
    cell_size = g.bp_settings['ec_cell_size']
    ec_strain = g.bp_settings['ec_strain']
    ec_steps = g.bp_settings['ec_steps']

    # Create original unedited structure and build neighbour list
    ec_rfkj_original_id = configs.add_common(structure, label, a0, uv, [cell_size], rcut, 'bp', 'e')
    g.bp[n]['ec_rfkj_original_id'] = ec_rfkj_original_id
    nl.build_nl(g.bp[n]['ec_rfkj_original_id'])

    ecs = 2 * ec_steps + 1

    g.bp[n]['ec_ids'] = []
    g.bp[n]['ec_strains'] = ds.bp_ec_arr(ecs)
    g.bp[n]['ec_volumes'] = ds.bp_ec_arr(ecs)
    g.bp[n]['ec_energies'] = ds.bp_ec_arr(ecs)

    g.bp[n]['c_ec'] = ds.bp_ec()
    g.bp[n]['c_ec_gpa'] = ds.bp_ec()
 
    for di in range(9):
      k = 0
      g.bp[n]['ec_ids'].append([])
      for i in range(-ec_steps, ec_steps+1):
        s = ec_strain * (i / ec_steps)
        g.bp[n]['ec_strains'][di, k] = s
        g.bp[n]['ec_ids'][di].append(configs.duplicate(ec_rfkj_original_id))
        k = k + 1




    """
    ! FOR COMPUTE VALUE CALCULATIONS
    ! Ravindran, Fast, Korzhavyi, Johansson 1998
    ! Density functional theory for calculation of elastic properties of 
    ! orthorhombic crystals: application to TiSi2
    ! 
    ######################################################"""






    """
    Rose equation of state
    Generate predicted energy vs lattice parameter data
    points for fitting
    (G. Bonny exact parameter fitting code)
    ######################################################"""
    
    if(g.bp[n]['rose_eos']):
      # Create rose configs
      rcut = g.bp_settings['rose_rcut']
      steps = g.bp_settings['rose_steps']
      cell_size = g.bp_settings['rose_cell_size']
      a0_start = a0 + g.bp_settings['rose_a0_d1']
      a0_mid = a0 - g.bp_settings['rose_a0_d1']
      a0_end = a0 + g.bp_settings['rose_a0_d2']
      xa = int(2/3 * steps)
      xb = steps - xa
      g.bp[n]['rose_a'] = ds.bp_rose_a_arr(steps)
      g.bp[n]['rose_e'] = ds.bp_rose_e_arr(steps)
      g.bp[n]['c_rose_e'] = ds.bp_rose_e_arr(steps)
      g.bp[n]['rose_ids'] = []
      g.bp[n]['rose_a'][:xa] = numpy.linspace(a0_start, a0_mid, xa)
      g.bp[n]['rose_a'][xa:] = numpy.linspace(a0_mid, a0_end, xb + 1)[1:]

      for i in range(steps):
        r_a0 = g.bp[n]['rose_a'][i]
        rose_id = configs.add_common(structure, label, r_a0, uv, [cell_size], rcut, 'bp', 'e')
        g.bp[n]['rose_ids'].append(rose_id)
        nl.build_nl(rose_id)

      # Fill in rose computed e
      known_a0 = g.bp[n]['a0']
      known_e0 = g.bp[n]['e0']
      known_v0 = g.bp[n]['v0']
      known_b0 = g.bp[n]['b0']

      # Calculate energy based on b0, a0, ecoh, vol
      g.bp[n]['rose_e'] = rose.calc(g.bp[n]['rose_a'], known_a0, known_b0, known_e0, known_v0)

    
    """
    ! For Fitting Points
    ! Mehl Singh Papaconstantopoulos strains 1993
    ! Properties of ordered intermetallic alloys: first-principles
    ! and approximate methods
    !
    ! Mehl Klein Papaconstantopoulos 1993
    ! First principles calculations of elastic properties of metals
    ######################################################"""

    if(g.bp[n]['msp_ec']):
      g.bp[n]['msp_ids'] = [[],[]]

      rcut = g.bp_settings['msp_rcut']
      cell_size = g.bp_settings['msp_cell_size']

      # set strains
      g.bp[n]['msp_strains'] = numpy.linspace(0.0, g.bp_settings['msp_ec_strain'], g.bp_settings['msp_ec_steps'])

      # set energy arrays
      g.bp[n]['msp_e'] = ds.bp_msp_e_arr(g.bp_settings['msp_ec_steps'])
      g.bp[n]['c_msp_e'] = ds.bp_msp_e_arr(g.bp_settings['msp_ec_steps'])
      
      # C11-C12 strains
      for i in range(g.bp_settings['msp_ec_steps']):
        uv_ortho = ds.uv_msp_orthorhombic(g.bp[n]['msp_strains'][i])
        msp_id = configs.add_common(structure, label, a0, uv_ortho, [cell_size], rcut, 'bp', 'e')
        g.bp[n]['msp_ids'][0].append(msp_id)

      # C44 strains
      for i in range(g.bp_settings['msp_ec_steps']):
        uv_mono = ds.uv_msp_monoclinic(g.bp[n]['msp_strains'][i])
        msp_id = configs.add_common(structure, label, a0, uv_mono, [cell_size], rcut, 'bp', 'e')
        g.bp[n]['msp_ids'][1].append(msp_id)

      
      # Calculate predicted energy-strain for C11-C12 and C44
      known_e0 = g.bp[n]['e0']
      known_v0 = g.bp[n]['v0']
      known_c11 = g.bp[n]['ec'][0,0]
      known_c12 = g.bp[n]['ec'][0,1]
      known_c44 = g.bp[n]['ec'][3,3]

      # Calculate energy based on b0, a0, ecoh, vol
      g.bp[n]['msp_e'][:,0] = msp.calc_c11_c12(g.bp[n]['msp_strains'], known_e0, known_v0, known_c11, known_c12)
      g.bp[n]['msp_e'][:,1] = msp.calc_c44(g.bp[n]['msp_strains'], known_e0, known_v0, known_c44)



    
    """
    ! Ravindran, Fast, Korzhavyi, Johansson 1998
    ! Density functional theory for calculation of elastic properties of 
    ! orthorhombic crystals: application to TiSi2
    ! 
    ######################################################"""

    




    # Mark as built
    g.bp[n]['built'] = True
   
    #configs.index()
    return True

  """
  Reset BP Configurations
  ######################################################"""
  @staticmethod
  def reset():
    output.log("Reset BP configurations", verbose=3)
    for n in range(len(g.bp)):
      configs.duplicate(g.bp[n]['eos_bm_original_id'], g.bp[n]['relax_id'])
      for k in range(len(g.bp[n]['eos_ids'])):
        configs.duplicate(g.bp[n]['eos_bm_original_id'], g.bp[n]['eos_ids'][k])
      for di in range(9):
        for k in range(len(g.bp[n]['ec_ids'][di])):
          configs.duplicate(g.bp[n]['ec_rfkj_original_id'], g.bp[n]['ec_ids'][di][k])


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
    ##################################################################
    
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
   
    # Make BM EoS points used to fit to the correct curve
    if(g.bp[n]['bm_eos']):
      # Volume arrays
      g.bp[n]['bm_eos_v'] = numpy.linspace(g.bp[n]['eos_volumes'][0], g.bp[n]['eos_volumes'][-1], g.bp_settings['bm_eos_points'])

      # Known arrays
      g.bp[n]['bm_eos_e_known'] = numpy.zeros((g.bp_settings['bm_eos_points'], 3,), dtype=numpy.float64)
      p = [g.bp[n]['v0'], g.bp[n]['e0'], g.bp[n]['b0'], 2.0]
      g.bp[n]['bm_eos_e_known'] = eos.bm_calc(p, g.bp[n]['bm_eos_v'])

      # Based on potential
      g.bp[n]['bm_eos_e_potential'] = numpy.zeros((g.bp_settings['bm_eos_points'], 3,), dtype=numpy.float64)
      p = [g.bp[n]['c_v0'], g.bp[n]['c_e0'], g.bp[n]['c_b0'], 2.0]
      g.bp[n]['bm_eos_e_potential'] = eos.bm_calc(p, g.bp[n]['bm_eos_v'])



      # Known arrays
      g.bp[n]['bm_eos_v_known_z'] = numpy.linspace(0.98 * g.bp[n]['v0'], 1.02 * g.bp[n]['v0'], g.bp_settings['bm_eos_points'])
      p = [g.bp[n]['v0'], g.bp[n]['e0'], g.bp[n]['b0'], 2.0]
      g.bp[n]['bm_eos_e_known_z'] = eos.bm_calc(p, g.bp[n]['bm_eos_v_known_z'])
      g.bp[n]['bm_eos_e_known_z'] = g.bp[n]['bm_eos_e_known_z'] - g.bp[n]['e0']
      g.bp[n]['bm_eos_v_known_z'] = g.bp[n]['bm_eos_v_known_z'] - g.bp[n]['v0']


      # Known arrays
      g.bp[n]['bm_eos_v_potential_z'] = numpy.linspace(0.98 * g.bp[n]['c_v0'], 1.02 * g.bp[n]['c_v0'], g.bp_settings['bm_eos_points'])
      p = [g.bp[n]['c_v0'], g.bp[n]['c_e0'], g.bp[n]['c_b0'], 2.0]
      g.bp[n]['bm_eos_e_potential_z'] = eos.bm_calc(p, g.bp[n]['bm_eos_v_potential_z'])
      g.bp[n]['bm_eos_e_potential_z'] = g.bp[n]['bm_eos_e_potential_z'] - g.bp[n]['c_e0']
      g.bp[n]['bm_eos_v_potential_z'] = g.bp[n]['bm_eos_v_potential_z'] - g.bp[n]['c_v0']
       





    # Compute EC
    ##################################################################


    if(g.bp[n]['rfkj_ec']):
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
   
      A = numpy.zeros((9,9,), dtype=numpy.float64)
      A[0,0] = v0/2
      A[1,1] = v0/2
      A[2,2] = v0/2
      A[3,3] = 2*v0
      A[4,4] = 2*v0
      A[5,5] = 2*v0
      A[6,6] = -1.0*v0
      A[7,7] = -1.0*v0
      A[8,8] = -1.0*v0
      A[6,0] = v0 / 2
      A[6,1] = v0 / 2
      A[7,0] = v0 / 2
      A[7,2] = v0 / 2
      A[8,1] = v0 / 2
      A[8,2] = v0 / 2
      C = numpy.asarray(d, dtype=numpy.float64)
      A_inv = numpy.linalg.inv(A)
      B = numpy.matmul(A_inv, C)

      g.bp[n]['c_ec'][0,0] = B[0]
      g.bp[n]['c_ec'][1,1] = B[1]
      g.bp[n]['c_ec'][2,2] = B[2]
      g.bp[n]['c_ec'][3,3] = B[3]
      g.bp[n]['c_ec'][4,4] = B[4]
      g.bp[n]['c_ec'][5,5] = B[5]
      g.bp[n]['c_ec'][0,1] = B[6]
      g.bp[n]['c_ec'][1,0] = B[6]
      g.bp[n]['c_ec'][0,2] = B[7]
      g.bp[n]['c_ec'][2,0] = B[7]
      g.bp[n]['c_ec'][1,2] = B[8]
      g.bp[n]['c_ec'][2,1] = B[8]

      # Calculate elastic constants
      """
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
      """

      # Calculate in GPA
      g.bp[n]['c_ec_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_ec'])
      

      # Compliance tensor
      g.bp[n]['c_sc'] = math.minverse(g.bp[n]['c_ec'])

      # Properties from Elastic Constants
      g.bp[n]['c_b0_r'] = bp.bulk_r(g.bp[n]['c_sc'])
      g.bp[n]['c_b0_v'] = bp.bulk_v(g.bp[n]['c_ec'])
      g.bp[n]['c_b0_avg'] = 0.5 * (g.bp[n]['c_b0_r'] + g.bp[n]['c_b0_v'])
      g.bp[n]['c_g_r'] = bp.shear_g_r(g.bp[n]['c_sc'])
      g.bp[n]['c_g_v'] = bp.shear_g_v(g.bp[n]['c_ec'])
      g.bp[n]['c_g_avg'] = 0.5 * (g.bp[n]['c_g_r'] + g.bp[n]['c_g_v'])
      g.bp[n]['c_e'] = bp.youngs_e(g.bp[n]['c_b0_avg'], g.bp[n]['c_g_avg'])
      g.bp[n]['c_melting_temperature'] = bp.melting_temperature(g.bp[n]['c_ec'])
      g.bp[n]['c_cubic_stability'] = bp.cubic_stability(g.bp[n]['c_ec'])
      g.bp[n]['c_stability'] = bp.stability(g.bp[n]['c_ec']) 

      # Convert to GPA
      g.bp[n]['c_b0_r_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_b0_r'])
      g.bp[n]['c_b0_v_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_b0_v'])
      g.bp[n]['c_b0_avg_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_b0_avg'])
      g.bp[n]['c_g_r_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_g_r'])
      g.bp[n]['c_g_v_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_g_v'])
      g.bp[n]['c_g_avg_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_g_avg'])
      g.bp[n]['c_e_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_e'])



    """
    Rose equation of state
    Run calculations to compare to known 
    energy-lattice parameter points
    ######################################################"""
    

    if(g.bp[n]['rose_eos']):
      for i in range(len(g.bp[n]['rose_ids'])):
        rose_id = g.bp[n]['rose_ids'][i]
        configs.calc_efs(rose_id)
        g.bp[n]['c_rose_e'][i] = g.configs[rose_id]['c_energy'] / g.configs[rose_id]['count']
    





    """
    MSP MKP strains C11-C12, C44
    Run calculations to compare to known 
    energy-lattice parameter points
    ######################################################"""

    
    if(g.bp[n]['msp_ec']):
      # C11-C12
      for i in range(len(g.bp[n]['msp_ids'][0])):
        msp_id = g.bp[n]['msp_ids'][0][i]
        configs.calc_efs(msp_id)
        g.bp[n]['c_msp_e'][i,0] = g.configs[msp_id]['c_energy'] / g.configs[msp_id]['count']
      # C44
      for i in range(len(g.bp[n]['msp_ids'][1])):
        msp_id = g.bp[n]['msp_ids'][1][i]
        configs.calc_efs(msp_id)
        g.bp[n]['c_msp_e'][i,1] = g.configs[msp_id]['c_energy'] / g.configs[msp_id]['count']    

      # Array to store temporary results
      ec_temp = numpy.zeros((3,), dtype=numpy.float64)

      # Fit polynomial for c11_c12
      s_c11_c12 = g.bp[n]['msp_strains'][:]
      e_c11_c12 = g.bp[n]['c_msp_e'][:,0]
      p_c11_c12 = numpy.polyfit(s_c11_c12, e_c11_c12, 2)  

      # Fit polynomial for c44
      s_c44 = g.bp[n]['msp_strains'][:]
      e_c44 = g.bp[n]['c_msp_e'][:,1]
      p_c44 = numpy.polyfit(s_c44, e_c44, 2)

      # Get bulk modulus and
      b0 = g.bp[n]['c_b0']
      v0 = g.bp[n]['c_v0']
      if(v0 == 0.0):
        v0 = 1.0e-20

      # Solve 
      a_mat= numpy.zeros((2,2,), dtype=numpy.float64)
      a_mat[0,0] = 1.0
      a_mat[0,1] = -1.0
      a_mat[1,0] = 1.0/3.0
      a_mat[1,1] = 2.0/3.0
      b_mat= numpy.zeros((2,), dtype=numpy.float64)
      b_mat[0] = p_c11_c12[0] / v0
      b_mat[1] = b0
      c_mat = numpy.linalg.solve(a_mat, b_mat)
     
      ec_temp[0] = c_mat[0]
      ec_temp[1] = c_mat[1]
      ec_temp[2] = (2 * p_c44[0]) / v0

      g.bp[n]['c_msp_ec'] = ds.bp_ec(ec_temp)
      g.bp[n]['c_msp_ec_gpa'] = units.convert('EV/ANG3', 'GPA', g.bp[n]['c_msp_ec'])


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


  @staticmethod
  def cubic_stability(ec): 
    if((ec[0,0] + 2 * ec[0,1]) > 0 and (ec[0,0] - ec[0,1]) > 0 and ec[3,3] > 0):
      return True
    return False


  @staticmethod
  def stability(ec): 
    if((ec[0,0] > 0.0) and (ec[1,1] > 0.0) and (ec[2,2] > 0.0) and 
       (ec[3,3] > 0.0) and (ec[4,4] > 0.0) and (ec[5,5] > 0.0) and
       (ec[0,0] + ec[1,1] - 2 * ec[0,1]) > 0.0 and
       (ec[0,0] + ec[2,2] - 2 * ec[0,2]) > 0.0 and
       (ec[1,1] + ec[2,2] - 2 * ec[1,2]) > 0.0 and
       (ec[0,0] + ec[1,1] + ec[2,2] + 2 * ec[0,1] + 2 * ec[0,2] + 2 * ec[1,2]) > 0):
      return True
    return False



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


