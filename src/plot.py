#!/bin/python3
"""
Output
"""

import os
import time
import numpy
import matplotlib.pyplot as plt

from f_toolbox import fnc as fnc

from g import g
from eos import eos
from rose import rose
from msp import msp



class plot:
  

  @staticmethod
  def potentials(dir = None):
    if(dir is None):
      dir = g.dir['plots_potential']
    for fn in range(len(g.potential['functions'])):
      pfn = g.potential['functions'][fn]
      plot_name = str(os.path.basename(pfn['file']))
      plot_name = plot_name.split(".")
      plot_name = plot_name[0] + ".eps"
      plot_path = os.path.join(dir, plot_name)

      if(pfn['tab'] is not None):
        plt.plot(pfn['tab'][:,0], pfn['tab'][:,2], color="tab:green", ls="dashed")
        plt.plot(pfn['tab'][:,0], pfn['tab'][:,1], color="tab:blue")
        #plt.yscale('log')
        plt.ylim(-5, 10)
        plt.savefig(plot_path)
        plt.cla()
      else:
        rcut = pfn['rcut']
        if(rcut is None):
          rcut = 10.0
        x = numpy.linspace(0.0, rcut, 1001)
        y = fnc.fv(pfn['fname'], x[:], pfn['p'], pfn['pf']) 
        dydx = fnc.fgradv(pfn['fname'], x[:], pfn['p'], pfn['pf']) 
        plt.plot(x[:], dydx[:], color="tab:green", ls="dashed")
        plt.plot(x[:], y[:], color="tab:blue")
        #plt.yscale('log')
        plt.ylim(-5, 10)
        plt.savefig(plot_path)
        plt.cla()
        

      

  @staticmethod
  def plot_eos(dir):
    for n in range(len(g.bp)):      
      path = os.path.join(dir, 'equation_of_state_bp_' + str(n) + '.eps')
      v = g.bp[n]['eos_volumes']
      e = g.bp[n]['eos_energies']
      fit = numpy.zeros((101,2), dtype=numpy.float64, order='F')
      fit[:,0] = numpy.linspace(v[0], v[-1], 101)
      p = numpy.asarray([g.bp[n]['c_v0'], g.bp[n]['c_e0'], g.bp[n]['c_b0'], g.bp[n]['c_b0p']])
      fit[:,1] = eos.bm_calc(p, fit[:,0])

      v0 = g.bp[n]['v0']
      e0 = g.bp[n]['e0']
      b0 = g.bp[n]['b0']
      vstart = min(v[0], 0.98 * v0)
      vend = max(v[-1], 1.02 * v0)
      bmv = numpy.linspace(vstart, vend, 1001)

      p = [v0, e0, b0, 2.0]
      bme2 = eos.bm_calc(p, bmv)
      p = [v0, e0, b0, 3.0]
      bme3 = eos.bm_calc(p, bmv)
      p = [v0, e0, b0, 4.0]
      bme4 = eos.bm_calc(p, bmv)

      plot.plot_eos_inner(path, v, e, fit, bmv, bme2, bme3, bme4)

      v[:] = v[:] - g.bp[n]['c_v0'] 
      e[:] = e[:] - g.bp[n]['c_e0'] 
      fit[:,0] = fit[:,0] - g.bp[n]['c_v0'] 
      fit[:,1] = fit[:,1] - g.bp[n]['c_e0'] 
      bmv[:] = bmv[:] - v0
      bme2[:] = bme2[:] - e0
      bme3[:] = bme3[:] - e0
      bme4[:] = bme4[:] - e0
      path = os.path.join(dir, 'equation_of_state_z_bp_' + str(n) + '.eps')

      plot.plot_eos_inner(path, v, e, fit, bmv, bme2, bme3, bme4)



  @staticmethod
  def plot_eos_inner(path, v, e, fit, bmv, bme2, bme3, bme4):
    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Equation of State')        
    plt.plot(v[:], e[:], color='k',  marker="x", ls='')
    plt.plot(fit[:,0], fit[:,1], color='k', ls='solid')
    plt.plot(bmv[:], bme2[:], color='r', ls='dashed')
    plt.plot(bmv[:], bme3[:], color='g', ls='dashed')
    plt.plot(bmv[:], bme4[:], color='b', ls='dashed')
    axs.set_xlabel('Volume (ang^3/atom)')
    axs.set_ylabel('Energy (ev/atom)')
    plt.savefig(path)
    plt.close('all')

      

  @staticmethod
  def plot_ec(dir):
    for n in range(len(g.bp)):      
      path = os.path.join(dir, 'elastic_strains_bp_' + str(n) + '.eps')
      fit = numpy.zeros((9,101,2), dtype=numpy.float64, order='F')
      for di in range(len(g.bp[n]['ec_ids'])):
        s = g.bp[n]['ec_strains'][di, :]
        e = g.bp[n]['ec_energies'][di, :]

        fit[di,:,0] = numpy.linspace(s[0], s[-1], 101)
        p = numpy.polyfit(s[:], e[:], 2)
        fit[di,:,1] = p[2] + p[1] * fit[di,:,0] + p[0] * fit[di,:,0]**2
      plot.plot_ec_inner(path, g.bp[n]['ec_strains'][:, :], g.bp[n]['ec_energies'][:, :], fit)


  @staticmethod
  def plot_ec_inner(path, s, e, fit):     
    plt.clf()
    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='xx-small')
    plt.rc('ytick', labelsize='xx-small')

    fig, axs = plt.subplots(3, 3, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Elastic Constant Curves')    
    
    for dn in range(9):
      
      xa = s[dn, :]
      ya = e[dn, :]
      xb = fit[dn, :, 0]
      yb = fit[dn, :, 1]
      
      axs[int(numpy.floor(dn/3)), dn % 3].plot(xa[:], ya[:], color='k',  marker="x", ls='')
      axs[int(numpy.floor(dn/3)), dn % 3].plot(xb[:], yb[:], color='k', ls='solid')
      axs[int(numpy.floor(dn/3)), dn % 3].set_title('Distortion D' + str(dn + 1))
      axs[int(numpy.floor(dn/3)), dn % 3].set_xlabel('Strain / a0')
      axs[int(numpy.floor(dn/3)), dn % 3].set_ylabel('Energy / e0')
    
    plt.savefig(path)
    plt.close('all')


        

      

  @staticmethod
  def plot_rose(dir):
    for n in range(len(g.bp)):      
      path = os.path.join(dir, 'rose_plot_bp_' + str(n) + '.eps')
      a = numpy.linspace(g.bp[n]['rose_a'][0], g.bp[n]['rose_a'][-1], 1001)
      e = g.bp[n]['rose_e']
      ca = g.bp[n]['rose_a']
      ce = g.bp[n]['c_rose_e']

      # Fill in rose computed e
      known_a0 = g.bp[n]['a0']
      known_e0 = g.bp[n]['e0']
      known_v0 = g.bp[n]['v0']
      known_b0 = g.bp[n]['b0']

      # Calculate energy based on b0, a0, ecoh, vol
      e = rose.calc(a, known_a0, known_b0, known_e0, known_v0)

      # Plot
      plot.plot_rose_inner(path, a, e, ca, ce)
    

  @staticmethod
  def plot_rose_inner(path, a, e, ca, ce):
    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle('Rose Equation of State')        
    plt.plot(ca[:], ce[:], color='k',  marker="x", ls='', label='This potential')
    plt.plot(a[:], e[:], color='k', ls='solid', label="Rose predicted energy-volume")
    plt.legend()
    axs.set_xlabel('Lattice Constant (ang)')
    axs.set_ylabel('Energy (ev/atom)')
    plt.savefig(path)
    plt.close('all')



  @staticmethod
  def plot_msp(dir):
    for n in range(len(g.bp)):   
      path_c11c12 = os.path.join(dir, 'msp_c11_c12_plot_bp_' + str(n) + '.eps')   
      path_c44 = os.path.join(dir, 'msp_c44_plot_bp_' + str(n) + '.eps') 
      path_c11c12_z = os.path.join(dir, 'msp_c11_c12_z_plot_bp_' + str(n) + '.eps')   
      path_c44_z = os.path.join(dir, 'msp_c44_z_plot_bp_' + str(n) + '.eps') 

      s = numpy.linspace(g.bp[n]['msp_strains'][0], g.bp[n]['msp_strains'][-1], 1001)
      cs = g.bp[n]['msp_strains']
      ce = g.bp[n]['c_msp_e']

      # Fill in rose computed e      
      known_e0 = g.bp[n]['e0']
      known_v0 = g.bp[n]['v0']
      known_c11 = g.bp[n]['ec'][0,0]
      known_c12 = g.bp[n]['ec'][0,1]
      known_c44 = g.bp[n]['ec'][3,3]

      # Calculate energy based on b0, a0, ecoh, vol
      e_c11c12 = msp.calc_c11_c12(s, known_e0, known_v0, known_c11, known_c12)
      e_c44 = msp.calc_c44(s, known_e0, known_v0, known_c44)

      # Plot
      plot.plot_msp_inner(path_c11c12, s, e_c11c12, cs, ce[:,0], 'MSP Cubic Elastic Constant C11-C12')
      plot.plot_msp_inner(path_c44, s, e_c44, cs, ce[:,1], 'MSP Cubic Elastic Constant C44')

      e_c11c12[:] = e_c11c12[:] - e_c11c12[0]
      e_c44[:] = e_c44[:] - e_c44[0]
      ce[:,:] = ce[:,:] - ce[0,:]

      # Plot
      plot.plot_msp_inner(path_c11c12_z, s, e_c11c12, cs, ce[:,0],'MSP Cubic Elastic Constant C11-C12 (zeroed)')
      plot.plot_msp_inner(path_c44_z, s, e_c44, cs, ce[:,1], 'MSP Cubic Elastic Constant C44 (zeroed)')
    

  @staticmethod
  def plot_msp_inner(path, a, e, ca, ce, title):
    plt.clf()    
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(12,9))
    fig.tight_layout(pad=5.0)
    fig.suptitle(title)        
    plt.plot(ca[:], ce[:], color='k',  marker="x", ls='')
    plt.plot(a[:], e[:], color='k', ls='solid')
    axs.set_xlabel('Strain')
    axs.set_ylabel('Energy (ev/atom)')
    plt.savefig(path)
    plt.close('all')