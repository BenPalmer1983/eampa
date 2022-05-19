#!/bin/python3

import numpy
import time
import os
import sys

from scipy.optimize import minimize
from scipy.optimize import basinhopping

from f_toolbox import atom
from f_toolbox import math
from f_toolbox import transforms

from g import g
from ds import ds
from display import display
from configs import configs
from nl import nl
from output import output
from plot import plot
from relax import relax
from eos import eos
from bp import bp
from units import units
from potential import potential
from rss import rss
from search_step import search_step
from simulated_annealing import simulated_annealing
from genetic_algorithm import ga
from random_search import random_search
from hybrid_search import hybrid_search
from wide_search import wide_search

class potfit:


  @staticmethod
  def run():
    # Log
    output.log("Potfit", verbose=0)
    
    # Make dirs
    potfit.make_dirs()

    # Turn off terminal log
    g.verbose['print'] = -1

    # Create storage dictionary
    potfit.set_potfit_store()

    # Parameters
    p = potential.get_p()

    
    if(p is None):
      output.log("No parameters to optimise", verbose=0)
      output.log("End", verbose=0)
      exit()

    g.potfit['potfit_start_time'] = time.time()
    for sn in range(len(g.potfit_steps)):
      if(sn > 0):
        p = ds.potential_p(g.potfit['p_best'])
      g.potfit['sn'] = sn
      g.potfit['sn_counter'] = 0

      
      g.potfit_steps[sn]['stats_rss'][0] = g.potfit['rss_best']
      g.potfit_steps[sn]['stats_start_time'] = time.time()
      g.potfit_steps[sn]['stats_start_counter'] = g.calc_counter


      # Run
      if(g.potfit_steps[sn]['type'].lower() == 'rand'):
        p = random_search.run(potfit.get_rss, p, niter=g.potfit_steps[sn]['niter'])
      elif(g.potfit_steps[sn]['type'].lower() == 'ws'):
        p = wide_search.run(potfit.get_rss, p)
      elif(g.potfit_steps[sn]['type'].lower() == 'sa'):
        p = simulated_annealing.run(potfit.get_rss, p, 
                                    niter=g.potfit_steps[sn]['niter'], titer=g.potfit_steps[sn]['titer'], tstart=g.potfit_steps[sn]['tstart'], tend=g.potfit_steps[sn]['tend'], pfact=g.potfit_steps[sn]['pfact'], pvar=g.potfit_steps[sn]['pvar'], vartype=g.potfit_steps[sn]['vartype'], gaussian=g.potfit_steps[sn]['gaussian'])
      elif(g.potfit_steps[sn]['type'].lower() == 'ga'):
        p = ga.run(potfit.get_rss, p, niter=g.potfit_steps[sn]['niter'], popsize=g.potfit_steps[sn]['popsize'], fresh=g.potfit_steps[sn]['fresh'], threshold=1.0e20, min_final=True)
      elif(g.potfit_steps[sn]['type'].lower() == 'hs'):
        p = hybrid_search.run(potfit.get_rss, p, 
                              pool_size=g.potfit_steps[sn]['poolsize'], 
                              sample_size=g.potfit_steps[sn]['samplesize'], 
                              minima_size=g.potfit_steps[sn]['minsize'])
      elif(g.potfit_steps[sn]['type'].lower() == 'bh'):
        res = basinhopping(potfit.get_rss, p, niter=g.potfit_steps[sn]['niter'], T=10.0, stepsize=1.0)
        p = res['x']
      elif(g.potfit_steps[sn]['type'].lower() == 'bfgs'):
        res = minimize(potfit.get_rss, p, method='BFGS', options={'gtol':  1.0e-6, 'maxiter': 4, })
        p = res['x']
      elif(g.potfit_steps[sn]['type'].lower() == 'cg'):
        res = minimize(potfit.get_rss, p, method='cg', options={'gtol':  1.0e-6, 'maxiter': 4, })
        p = res['x']

      # Record stats
      g.potfit_steps[sn]['stats_rss'][1] = g.potfit['rss_best']
      g.potfit_steps[sn]['stats_time'] = time.time() - g.potfit_steps[sn]['stats_start_time']
      g.potfit_steps[sn]['stats_counter'] = g.calc_counter - g.potfit_steps[sn]['stats_start_counter']
      g.potfit_steps[sn]['stats_speed'] = g.potfit_steps[sn]['stats_counter'] / g.potfit_steps[sn]['stats_time']
      g.potfit_steps[sn]['stats_complete'] = True


    # End fitting
    potfit.end()



  @staticmethod
  def end(interupt=False): 
    g.potfit['potfit_end_time'] = time.time()

    if(interupt):
      sn = g.potfit['sn']
      g.potfit_steps[sn]['stats_rss'][1] = g.potfit['rss_best']
      g.potfit_steps[sn]['stats_time'] = time.time() - g.potfit_steps[sn]['stats_start_time']
      g.potfit_steps[sn]['stats_counter'] = g.calc_counter - g.potfit_steps[sn]['stats_start_counter']
      g.potfit_steps[sn]['stats_speed'] = g.potfit_steps[sn]['stats_counter'] / g.potfit_steps[sn]['stats_time']
      g.potfit_steps[sn]['stats_complete'] = False


    output.log("##################################################", verbose=0)
    output.log("Ending of Potential Fit", verbose=0)
    output.log("##################################################", verbose=0)
    output.log("", verbose=0)
    if(interupt):
      output.log("Reached maximum time limit", verbose=0)
    else:
      output.log("Completed within maximum time limit", verbose=0)
    output.log("", verbose=0)

    # Set parameters for potential
    p = ds.potential_p(g.potfit['p_best']) 
    potential.update_p(p)

    # Run calculation
    efs_rss, bp_rss, total_rss = potfit.get_rss_calc()
    output.log("", verbose=0)


    output.log("Best Parameters", verbose=0)
    output.log("--------------------------------------------------", verbose=0)
    line = ''
    for i in range(len(g.potfit['p_best'])):
      line = line + "{:10.3e} ".format(g.potfit['p_best'][i])
      if((i + 1) % 10 == 0 and i != len(g.potfit['p_best'])-1):
        output.log(line, verbose=0)
        line = ''
    if(line != ''):   
      output.log(line, verbose=0)
 
    t = g.potfit['potfit_end_time'] - g.potfit['potfit_start_time']
    output.log("Stats", verbose=0)
    output.log("--------------------------------------------------", verbose=0)
    output.log("{:12s} {:16.2f}".format("Time: ", t), verbose=0)
    output.log("{:12s} {:16d} {:16.2f}".format("Steps: ", g.potfit['counter'], g.potfit['counter']/ t), verbose=0)
    output.log("{:12s} {:16d} {:16.2f}".format("Configs: ", g.calc_counter, g.calc_counter/ t), verbose=0)
    output.log("{:12s} {:16.2f}".format("Best RSS: ", g.potfit['rss_best']), verbose=0)
    output.log('', verbose=0)

    # Make tabulated
    potential.tab_for_output()

    # Output potentials - analytic, tab, and plots
    output.potentials(g.dir['fit_best_pot'])
    output.potentials_tabulated(g.dir['fit_best_pottab'])
    plot.potentials(g.dir['fit_best_plot'])

    # Output bp and other plots
    output.bp(g.dir['fit_best'])
    plot.plot_eos(g.dir['fit_best_plots'])
    plot.plot_ec(g.dir['fit_best_plots'])
    plot.plot_rose(g.dir['fit_best_plots'])
    plot.plot_msp(g.dir['fit_best_plots'])

    output.potfit_summary(g.dir['fit'])
    
 
    
    # Turn on terminal log
    g.verbose['print'] = 1
    exit()
    sys.exit()

  
  
  @staticmethod
  def make_dirs():    

    g.dir['fit_best'] = os.path.join(g.dir['fit'], "best")
    g.dir['fit_best_pot'] = os.path.join(g.dir['fit_best'], "pot")
    g.dir['fit_best_pottab'] = os.path.join(g.dir['fit_best'], "pot_tab")
    g.dir['fit_best_plot'] = os.path.join(g.dir['fit_best'], "pot_plot")
    g.dir['fit_best_plots'] = os.path.join(g.dir['fit_best'], "plots")

    g.dir['fit_temp_best'] = os.path.join(g.dir['fit'], "temp")

    for k in g.dir.keys():
      os.makedirs(g.dir[k], exist_ok=True)

  

  @staticmethod
  def temporary_save(): 
    tp = g.potfit_bestsaveperiod
    dirroot = g.dir['fit_temp_best']
    tn = str(int(tp * numpy.ceil((time.time() - g.start) / tp)))
    while(len(tn) < 8):
      tn = "0" + tn

    tndir = os.path.join(dirroot, tn)
    tndir_pot = os.path.join(tndir, "pot")
    tndir_pottab = os.path.join(tndir, "pot_tab")
    tndir_plot = os.path.join(tndir, "pot_plot")
    os.makedirs(tndir, exist_ok=True)
    os.makedirs(tndir_pot, exist_ok=True)
    os.makedirs(tndir_pottab, exist_ok=True)
    os.makedirs(tndir_plot, exist_ok=True)

    output.potentials(tndir_pot)
    potential.tab_for_output()
    output.potentials_tabulated(tndir_pottab)
    plot.potentials(tndir_plot)
    output.bp(tndir_pot)


  # Function to minimise
  @staticmethod
  def get_rss(p_in):
    # End
    if((time.time() - g.start) > g.maxtime):
      potfit.end(interupt=True)

    # Set parameters for potential
    p = ds.potential_p(p_in) 
    potential.update_p(p)

    # Run calculation
    efs_rss, bp_rss, total_rss = potfit.get_rss_calc()
    
    # Save results
    g.potfit['counter'] = g.potfit['counter'] + 1
    g.potfit['improvement_counter'] = g.potfit['improvement_counter'] + 1
    g.potfit['sn_counter'] = g.potfit['sn_counter'] + 1
    g.potfit_steps[g.potfit['sn']]['counter'] = g.potfit_steps[g.potfit['sn']]['counter'] + 1
    g.potfit['p_current'] = numpy.copy(p)
    g.potfit['efs_rss_current'] = efs_rss
    g.potfit['bp_rss_current'] = bp_rss
    g.potfit['rss_current'] = total_rss

    for n in range(len(g.bp)):
      g.potfit['bp_current'][n]['a0'] = g.bp[n]['c_a0']
      g.potfit['bp_current'][n]['e0'] = g.bp[n]['c_e0']
      g.potfit['bp_current'][n]['v0'] = g.bp[n]['c_v0']
      g.potfit['bp_current'][n]['b0_gpa'] = g.bp[n]['c_b0_gpa']
      g.potfit['bp_current'][n]['ec_gpa'] = numpy.copy(g.bp[n]['c_ec_gpa'])
      g.potfit['bp_current'][n]['mspec_gpa'] = numpy.copy(g.bp[n]['c_msp_ec_gpa'])
      g.potfit['bp_current'][n]['rss_details'] = g.bp[n]['rss_details']

    # Save best 
    if((g.potfit['rss_best'] is None or total_rss < g.potfit['rss_best']) and not numpy.isnan(total_rss)):
      g.potfit['improvement_counter'] = 0
      g.potfit['p_best'] = numpy.copy(p)
      g.potfit['efs_rss_best'] = efs_rss
      g.potfit['bp_rss_best'] = bp_rss
      g.potfit['rss_best'] = total_rss

      # Save best values
      for n in range(len(g.bp)):
        g.potfit['bp_best'][n]['a0'] = g.bp[n]['c_a0']
        g.potfit['bp_best'][n]['e0'] = g.bp[n]['c_e0']
        g.potfit['bp_best'][n]['v0'] = g.bp[n]['c_v0']
        g.potfit['bp_best'][n]['b0_gpa'] = g.bp[n]['c_b0_gpa']
        g.potfit['bp_best'][n]['ec_gpa'] = numpy.copy(g.bp[n]['c_ec_gpa'])
        g.potfit['bp_best'][n]['mspec_gpa'] = numpy.copy(g.bp[n]['c_msp_ec_gpa'])

      potfit.temporary_save()



    # Display results
    display.show()

    # Log
    output.log("Total rss: " + str(total_rss) + "  " + str(g.potfit['rss_best']), verbose=2)
    output.log("Total rss: " + str(p), verbose=2)

    # 
    return total_rss



  @staticmethod
  def get_rss_calc():    
    efs_rss = 0.0
    if(configs.count('efs') > 0 and g.rss_w['configs'] != 0.0):
      configs.calc('efs')
      efs_rss = configs.rss(config_tag='efs')
    bp_rss = 0.0
    if(len(g.bp) > 0 and g.rss_w['bp'] != 0.0):  
      bp_rss = bp.run_rss()
    total_rss = bp_rss + efs_rss


    return efs_rss, bp_rss, total_rss



  @staticmethod
  def set_potfit_store():
    g.potfit = ds.potfit()

    # BP current
    if(g.potfit['bp_current'] is None):
      g.potfit['bp_current'] = []
      for n in range(len(g.bp)):
        g.potfit['bp_current'].append(ds.potfit_bp())

    # BP best
    if(g.potfit['bp_best'] is None):
      g.potfit['bp_best'] = []
      for n in range(len(g.bp)):
        g.potfit['bp_best'].append(ds.potfit_bp())


  
