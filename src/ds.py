#!/bin/python3
"""
Dictionaries and Lists
"""

import numpy

from f_toolbox import transforms

class ds:


  @staticmethod
  def runtypes():
    return ['fit','rss','efs','bp','gauge','test']


  @staticmethod
  def units():
    return {
    'energy': 'ev', 
    'length': 'ang', 
    'pressure': 'gpa',
    }


  @staticmethod
  def stats():
    return {
    'config_calc_count': 0,
    'nl_build_count': 0,
    }
    



  """#####################################
  INPUT FILE
  #####################################"""

  @staticmethod
  def ifile():
    return {
    'runtype': 'rss',
    'configs': [],
    'nlsize': 100000,
    'potindex': None,
    }

  
  """#####################################
  CONFIG
  #####################################"""

  # Converted to ang, ev, ev/ang, ev/ang3 when loaded
  @staticmethod
  def config():
    return {
    'config_id': None,
    'nl_id': None,
    'tag': None,
    'file_path': None,  
    'a0': 0.0,  
    'uv': ds.config_uv(None),
    'cxyz': None,
    'volume': 0.0,
    'w': 1.0,
    'rcut': 0.0,
    'count': 0,
    'labels': None,
    'coords': None,
    'energy': None,
    'forces': None,
    'stress': None,
    'c_energy': None,
    'c_forces': None,
    'c_stress': None,
    'c_stress_gpa': None,
    'c_pair_embedding_energy': numpy.zeros((3,),dtype=numpy.float64, order='F'),
    'c_forces_total': None,
    'calc_type': None, # None, e, ef, efs  (force a calculation type if needed)
    'last_calc_type': None,
    'rss': None,
    'rss_details': None,
    'timer': 0.0,
    'counter': 0,
    }    

  @staticmethod
  def config_uv(uv_in=None):
    if(uv_in is None):
      return numpy.zeros((3,3,),dtype=numpy.float64, order='F')
    uv_in = numpy.asarray(uv_in, dtype=numpy.float64, order='F')

    if(uv_in.size == 1):       
      uv = numpy.zeros((3,3,),dtype=numpy.float64, order='F')
      uv[0,0] = uv_in[0]
      uv[1,1] = uv_in[0]
      uv[2,2] = uv_in[0]
      return uv
    elif(uv_in.size == 3):
      uv = numpy.zeros((3,3,),dtype=numpy.float64, order='F')
      uv[0,0] = uv_in[0]
      uv[1,1] = uv_in[1]
      uv[2,2] = uv_in[2]
      return uv
    elif(uv_in.size == 6):
      uv = numpy.zeros((3,3,),dtype=numpy.float64, order='F')
      uv[0,0] = uv_in[0]
      uv[1,1] = uv_in[1]
      uv[2,2] = uv_in[2]
      uv[0,1] = uv_in[5]
      uv[1,0] = uv_in[5]
      uv[0,2] = uv_in[4]
      uv[2,0] = uv_in[4]
      uv[1,2] = uv_in[3]
      uv[2,1] = uv_in[3]
      return uv
    elif(uv_in.size == 9 and uv_in.ndim == 1):
      uv = numpy.zeros((3,3,),dtype=numpy.float64, order='F')
      uv[0,:] = uv_in[0:3] 
      uv[1,:] = uv_in[3:6] 
      uv[2,:] = uv_in[6:9] 
      return uv
    elif(uv_in.size == 9 and uv_in.ndim == 2):
      return uv_in


  @staticmethod
  def config_labels(labels):
    return numpy.asarray(labels, dtype=numpy.int32, order='F') 

  @staticmethod
  def config_labels_empty(s):
    return numpy.zeros((s), dtype=numpy.int32, order='F') 

  @staticmethod
  def config_coords(coords):
    return numpy.asarray(coords, dtype=numpy.float64, order='F') 

  @staticmethod
  def config_coords_empty(s):
    return numpy.zeros((s, 3), dtype=numpy.float64, order='F') 

  @staticmethod
  def config_forces(forces):
    return numpy.asarray(forces, dtype=numpy.float64, order='F') 

  @staticmethod
  def config_forces_empty(s):
    return numpy.zeros((s, 3), dtype=numpy.float64, order='F')

  @staticmethod
  def config_stress(stress):
    return numpy.asarray(stress, dtype=numpy.float64, order='F') 

  @staticmethod
  def config_stress_empty():
    return numpy.zeros((3, 3), dtype=numpy.float64, order='F')


  @staticmethod
  def config_cxyz(cxyz):
    cxyz = numpy.asarray(cxyz, dtype=numpy.int32, order='F')
    if(cxyz.size == 3):
      return cxyz
    if(cxyz.size == 1):
      cxyz_out = numpy.zeros((3), dtype=numpy.int32, order='F') 
      cxyz_out[0] = cxyz[0]
      cxyz_out[1] = cxyz[0]
      cxyz_out[2] = cxyz[0]
      return cxyz_out

  @staticmethod
  def config_rss_details():
    return {
    'config_total': 0.0, 
    'energy': 0.0,
    'forces': 0.0,
    'stress': 0.0,
    }


  @staticmethod
  def uv_msp_orthorhombic(s):
    uv = numpy.zeros((3,3,),dtype=numpy.float64, order='F')
    uv = transforms.msp_orthorhombic(uv, s)
    return uv


  @staticmethod
  def uv_msp_monoclinic(s):
    uv = numpy.zeros((3,3,),dtype=numpy.float64, order='F')
    uv = transforms.msp_monoclinic(uv, s)
    return uv





  """
  @staticmethod
  def config_cxyz(cxyz):
    cxyz = numpy.asarray(cxyz, dtype=numpy.int32, order='F')
    if(cxyz.size == 3):
      return cxyz
    if(cxyz.size == 1):
      cxyz_out = numpy.zeros((3), dtype=numpy.int32, order='F') 
      cxyz_out[0] = cxyz[0]
      cxyz_out[1] = cxyz[0]
      cxyz_out[2] = cxyz[0]
      return cxyz_out
  """




  """#####################################
  NEIGHBOUR LIST
  #####################################"""

  @staticmethod
  def nl():
    return {
    'nl_id': None,
    'config_id': None,
    'count': 0,
    'labels': None,
    'r': None,
    'rvec': None,
    'inhalo': None,
    'timer': 0.0,
    'counter': 0,
    }

  @staticmethod
  def nl_lout(s):
    return numpy.zeros((s, 4), dtype=numpy.int32, order='F')

  @staticmethod
  def nl_r(s):
    return numpy.zeros((s), dtype=numpy.float64, order='F')

  @staticmethod
  def nl_rvec(s):
    return numpy.zeros((s, 3), dtype=numpy.float64, order='F')

  @staticmethod
  def nl_inhalo(s):
    return numpy.zeros((s), dtype=numpy.int32, order='F')




  """#####################################
  POTENTIAL
  #####################################"""

  @staticmethod
  def potential():
    return {
    'index_file': '',
    'functions': [],
    'tab_for_output': [],
    'tab_transformed': [],
    }


  @staticmethod
  def potential_function():
    return {
    'fon': None,  
    'file': None,
    'a': -1,
    'b': -1,
    'rcut': None,
    'ftype': None,            # 1=PAIR, 2=DENS, 3=EMBE       
    'fgroup': -1,
    'fname': None,
    'p': None,
    'pf': None,
    'tab': None,
    'fortran_pkey': None,
    'fix': False,  
    }

  @staticmethod
  def tab(tlen=1001):
    return numpy.zeros((tlen, 3), dtype=numpy.float64, order='F')


  @staticmethod
  def potential_p(p_in):
    if(len(p_in) == 0):
      return None
    p_arr = numpy.asarray(p_in, dtype=numpy.float64, order='F')    
    return numpy.copy(p_arr, order='F')




   
  """#####################################
  RSS
  #####################################"""

  @staticmethod
  def rss_w():
    return {
    'configs': 0.001,
    'e': 1.0,
    's': 1.0,
    'f': 1.0,
    'bp': 1.0,
    'a0': 10.0,
    'e0': 10.0,
    'b0': 100.0,
    'ec': 25.0,
    'rose': 1.0,
    'msp': numpy.asarray([10.0, 10.0]),         # shape, differences
    'bm_eos': numpy.asarray([10.0, 10.0]),      # shape, differences
    }

  @staticmethod
  def bp_rss_details():
    return {
    'bp_total': 0.0, 
    'a0': 0.0,
    'e0': 0.0,
    'b0': 0.0,
    'ec': 0.0,
    'rose': 0.0,
    'msp_shape': 0.0,
    'msp_values': 0.0,
    'bm_eos_shape': 0.0, 
    'bm_eos_values': 0.0, 
    'bp_penalty': 0.0, 
    }


  """#####################################
  BULK PROPERTY CALCULATIONS
  #####################################"""


  @staticmethod
  def bp():
    return {
    'built': False,
    'label': None,
    'structure': None,
    'weight': 1.0,
    'rss': None,
    'rss_details': None,
    'a0': None,
    'volrelaxed': None,
    'uv': None,
    'v0': None,        # volume per atom
    'e0': None,
    'b0': None,
    'b0_gpa': None,
    'ec': None,
    'ec_gpa': None,
    'original_id': None,   # remove
    'eos_bm_original_id': None,
    'ec_rfkj_original_id': None,
    'relax_id': None,
    'eos_ids': None,
    'eos_strains': None,
    'eos_volumes': None,
    'eos_energies': None,
    'ec_ids': None,
    'ec_strains': None,
    'ec_volumes': None,
    'ec_energies': None,
    'c_a0': None,
    'c_uv': None,
    'c_v0': None,
    'c_e0': None,
    'c_b0': None,
    'c_b0_gpa': None,
    'c_b0p': None,
    'c_ec': None,        # stiffness tensor - elastic constants   (calculated)
    'c_ec_gpa': None,
    'c_sc': None,        # compliance tensor - inverse ec         (calculated)
    'c_b0_r': None,
    'c_b0_v': None,
    'c_b0_avg': None,
    'c_g_r': None,
    'c_g_v': None,
    'c_g_avg': None, 
    'c_e': None, 
    'c_melting_temperature': None,  
    'c_cubic_stability': None,  
    'c_stability': None,  
    'c_b0_r_gpa': None,
    'c_b0_v_gpa': None,
    'c_b0_avg_gpa': None,
    'c_g_r_gpa': None,
    'c_g_v_gpa': None,
    'c_g_avg_gpa': None, 
    'c_e_gpa': None, 
    # Relax

    #BM EoS fitting
    'bm_eos': True,
    'bm_eos_v': None,
    'bm_eos_e_known': None,
    'bm_eos_e_potential': None,
    'bm_eos_v_known_z': None,
    'bm_eos_e_known_z': None,        # Set e0 to 0.0
    'bm_eos_v_potential_z': None,
    'bm_eos_e_potential_z': None,    # Set e0 to 0.0
    # ROSE equation for fitting
    'rose_eos': True,
    'rose_ids': None,
    'rose_a': None,
    'rose_e': None,
    'c_rose_e': None,
    # MSP cubic strains for fitting    
    'msp_ec': True,
    'msp_ids': None,
    'msp_strains': None,
    'msp_e': None,
    'c_msp_e': None,
    'c_msp_ec': None,         # 6x6 elastic constants here
    'c_msp_ec_gpa': None,     # 6x6 elastic constants here
    # RFKJ orthorhombic strains for fitting  
    'rfkj_ec': True,
    'rfkj_ids': None,
    'rfkj_strains': None,
    'rfkj_e': None,
    'c_rfkj_e': None,
    
    }

  @staticmethod
  def bp_settings():
    return {
    'cell_size': 2,
    'rcut': 7.5,
    # RELAX TYPE
    'relax': 'a0',   # 'a0' or 'uva0'  
    # EOS BM
    'eos_strain': 0.03,
    'eos_steps': 4,        #  2 * eos_steps + 1
    'eos_cell_size': 2,
    'eos_rcut': 7.5,
    # BM EC
    'ec_strain': 0.03,
    'ec_steps': 4,         #  2 * ec_steps + 1
    'ec_cell_size': 2,
    'ec_rcut': 7.5,
    # BM EOS COMPARISON
    'bm_eos_points': 20,
    #ROSE
    'rose_rcut': 7.5,
    'rose_cell_size': 4,
    'rose_a0_d1': -1.0,
    'rose_a0_d2': 4.0,
    'rose_steps': 21,
    #MSP
    'msp_rcut': 7.5, 
    'msp_cell_size': 2,  
    'msp_ec_strain': 0.03,
    'msp_ec_steps': 11,
    }


  @staticmethod
  def bp_eos_arr(s):
    return numpy.zeros((s,), dtype=numpy.float64)


  @staticmethod
  def bp_ec_arr(s):
    return numpy.zeros((9, s,), dtype=numpy.float64)


  @staticmethod
  def bp_d():
    return numpy.zeros((9,), dtype=numpy.float64)

  @staticmethod
  def bp_ec(ec=None):
    if(ec is None):
      return numpy.zeros((6, 6,), dtype=numpy.float64, order='F')
    ec = numpy.asarray(ec, dtype=numpy.float64, order='F')
    ec_out = numpy.zeros((6, 6,), dtype=numpy.float64, order='F')
    if(ec.ndim == 1 and ec.size == 3):
      ec_out[0,0] = ec[0]
      ec_out[1,1] = ec[0]
      ec_out[2,2] = ec[0]
      ec_out[0,1] = ec[1]
      ec_out[0,2] = ec[1]
      ec_out[1,2] = ec[1]
      ec_out[1,0] = ec[1]
      ec_out[2,0] = ec[1]
      ec_out[2,1] = ec[1]
      ec_out[3,3] = ec[2]
      ec_out[4,4] = ec[2]
      ec_out[5,5] = ec[2]
    elif(ec.ndim == 1 and ec.size == 9):
      ec_out[0,0] = ec[0]
      ec_out[1,1] = ec[1]
      ec_out[2,2] = ec[2]
      ec_out[3,3] = ec[3]
      ec_out[4,4] = ec[4]
      ec_out[5,5] = ec[5]
      ec_out[1,2] = ec[6]
      ec_out[0,2] = ec[7]
      ec_out[0,1] = ec[8]
      ec_out[2,1] = ec[6]
      ec_out[2,0] = ec[7]
      ec_out[1,0] = ec[8]
    elif(ec.ndim == 2 and ec.size == 36):
      ec_out[:,:] = ec[:,:]
    return ec_out


  @staticmethod
  def bp_rose_a_arr(s):
    return numpy.zeros((s,), dtype=numpy.float64)


  @staticmethod
  def bp_rose_e_arr(s):
    return numpy.zeros((s,), dtype=numpy.float64)


  @staticmethod
  def bp_msp_s_arr(s):
    return numpy.zeros((s,), dtype=numpy.float64)


  @staticmethod
  def bp_msp_e_arr(s):
    return numpy.zeros((s,2,), dtype=numpy.float64)

   
  """#####################################
  potfit
  #####################################"""

  @staticmethod
  def potfit():
    return {
    'counter': 0,
    'sn': 0,
    'sn_counter': 0,
    'improvement_counter': 0,
    'p_current': None,
    'efs_rss_current': None,
    'bp_rss_current': None,
    'rss_current': None,
    'bp_current': None,
    'p_best': None,
    'efs_rss_best': None,
    'bp_rss_best': None,
    'rss_best': None,
    'bp_best': None,
    'potfit_start_time': None,
    'potfit_end_time': None,
    }

  @staticmethod
  def potfit_bp():
    return {
    'a0': None,
    'v0': None,
    'e0': None,
    'b0_gpa': None,
    'ec_gpa': None,
    }


  @staticmethod
  def potfit_step():
    return {
    'type': None,
    'counter': 0,
    'niter': 1000,
    'titer': 10,
    'tstart': 10.0,
    'tend': 0.01,
    'pfact': 0.7,
    'pvar': None,
    'vartype': None,
    'gaussian': False,
    'popsize': 1000,
    'fresh': 0.2,
    'search': None,
    'minsize': None,
    'poolsize': None,
    'samplesize': None,
    'minsize': None,
    # stats
    'stats_rss': [None, None,],
    'stats_time': None,
    'stats_counter': None,
    'stats_speed': None,
    'stats_complete': None,
    # Temp    
    'stats_start_time': None,
    'stats_start_counter': None,
    }

  @staticmethod
  def potfit_step_types():
    return {
    'sa': 'Simulated Annealing',
    'ga': 'Genetic Algorithm',
    'hs': 'Hybrid Search',
    'bh': 'Basin Hopping',
    'bfgs': 'BFGS',
    'cg': 'CG',
    'rand': 'Random Search',
    'ws': 'Wide Search',
    }

   
  """#####################################
  Gauge
  #####################################"""

  @staticmethod
  def gauge():
    return {
    'label': None,
    'structures': None,
    'a0': None,
    'uv': None,
    'size': None,
    'rcut': None,
    }

