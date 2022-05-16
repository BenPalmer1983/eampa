#!/bin/python3

import numpy
import time
import os

from g import g
from ds import ds

class display:

  w = 150

  @staticmethod
  def show():
    if(g.display == 1):
      display.show_1()
    elif(g.display == 2):
      display.show_2()


  @staticmethod
  def show_1():
    display.clear()
    if(g.potfit['rss_best'] is None):
      print("waiting for first result...")
      return False

    # Save some vars
    rss_best = 1.0e99
    bp_rss_best = 1.0e99
    efs_rss_best = 1.0e99
    rss_current = 1.0e99
    bp_rss_current = 1.0e99
    efs_rss_current = 1.0e99
    if(g.potfit['rss_best'] is not None): rss_best = g.potfit['rss_best']
    if(g.potfit['bp_rss_best'] is not None): bp_rss_best = g.potfit['bp_rss_best']
    if(g.potfit['efs_rss_best'] is not None): efs_rss_best = g.potfit['efs_rss_best']
    if(g.potfit['rss_current'] is not None): rss_current = g.potfit['rss_current']
    if(g.potfit['bp_rss_current'] is not None): bp_rss_current = g.potfit['bp_rss_current']
    if(g.potfit['efs_rss_current'] is not None): efs_rss_current = g.potfit['efs_rss_current']
      

    # Display

    display.hrr() ######################
    print(" Fitting ")
    display.hrr() ######################
    print(" {0:24s} {1:12d}     {2:20s} {3:12d}     {4:20s} {5:12d}".format(
          "Counter:", g.potfit['counter'], 
          "Since improvement:", g.potfit['improvement_counter'],
          "Configs/s:", int(g.calc_counter / g.calc_time))
         )    
    print(" {0:24s} {1:12.6f}     {2:20s} {3:12.6f}     {4:20s} {5:12d}".format(
          "Runtime:", time.time() - g.start, 
          "Max runtime:", g.maxtime,
          "P Size:", len(g.potfit['p_current']))
         )  
    print(" {0:24s} {1:12.6e}     {2:20s} {3:12.6e}     {4:20s} {5:12.6e}".format(
          "Rss:", rss_current, 
          "BP rss:", bp_rss_current,
          "EFS rss:", efs_rss_current)
         )
    print(" {0:24s} {1:12.6e}     {2:20s} {3:12.6e}     {4:20s} {5:12.6e}".format(
          "Rss (best):", rss_best, 
          "BP rss (best):", bp_rss_best,
          "EFS rss (best):", efs_rss_best)
         )
    display.hrr() ######################

    for sn in range(len(g.potfit_steps)):
      fit_type = g.potfit_step_types[g.potfit_steps[sn]['type'].lower()]
      print("{0:2d}  {1:24s}   {2:6d}".format(sn + 1, fit_type, g.potfit_steps[sn]['counter']), end="")
      if(g.potfit['sn'] == sn):
        print(" <--   ", end="")
      print()
    
    display.hrr() ######################
    print("> ",g.displaynote)
      
    display.hrr() ######################
    for n in range(len(g.bp)):
      print("")
      print(" BP " + str(n))
      print("{:15s}{:15s}{:15s}{:15s}{:10s}{:15s}{:15s}{:15s}".format('','Known','Best','Current','','Known','Best','Current'))
      display.hr() #------------------------
      print("BM EoS")
      print(" {0:5s}  {1:12.2f}   {2:12.2f}   {3:12.2e}      {4:5s}  {5:12.2f}   {6:12.2f}   {7:12.2e}".format(
            "a0", g.bp[n]['a0'], g.potfit['bp_best'][n]['a0'], g.potfit['bp_current'][n]['a0'],
            "e0", g.bp[n]['e0'], g.potfit['bp_best'][n]['e0'], g.potfit['bp_current'][n]['e0']
           ))
      print(" {:5s}  {:12.2f}   {:12.2f}   {:12.2e}      {:5s}  {:12.2f}   {:12.2f}   {:12.2e}".format(
            "b0", g.bp[n]['b0_gpa'], g.potfit['bp_best'][n]['b0_gpa'], g.potfit['bp_current'][n]['b0_gpa'],
            "v0", g.bp[n]['v0'], g.potfit['bp_best'][n]['v0'], g.potfit['bp_current'][n]['v0']
           ))

      if(g.bp[n]['msp_ec']):
        display.hr() #------------------------
        print("MSP/MKP Cubic Crystal Strain")
        print(" c11   {0:12.2f}    {1:12.2f}    {2:12.2e}".format(g.bp[n]['ec_gpa'][0,0], g.potfit['bp_best'][n]['mspec_gpa'][0,0], g.potfit['bp_current'][n]['mspec_gpa'][0,0]))
        print(" c12   {0:12.2f}    {1:12.2f}    {2:12.2e}".format(g.bp[n]['ec_gpa'][0,1], g.potfit['bp_best'][n]['mspec_gpa'][0,1], g.potfit['bp_current'][n]['mspec_gpa'][0,1]))
        print(" c44   {0:12.2f}    {1:12.2f}    {2:12.2e}".format(g.bp[n]['ec_gpa'][3,3], g.potfit['bp_best'][n]['mspec_gpa'][3,3], g.potfit['bp_current'][n]['mspec_gpa'][3,3]))
      display.hr() #------------------------
      print("RFKJ Orthorhombic Crystal Strain")
      print(" {0:5s}  {1:12.2f}   {2:12.2f}   {3:12.2e}      {4:5s}  {5:12.2f}   {6:12.2f}   {7:12.2e}".format(
            "c11", g.bp[n]['ec_gpa'][0,0], g.potfit['bp_best'][n]['ec_gpa'][0,0], g.potfit['bp_current'][n]['ec_gpa'][0,0],
            "c22", g.bp[n]['ec_gpa'][1,1], g.potfit['bp_best'][n]['ec_gpa'][1,1], g.potfit['bp_current'][n]['ec_gpa'][1,1]
           ))
      print(" {0:5s}  {1:12.2f}   {2:12.2f}   {3:12.2e}      {4:5s}  {5:12.2f}   {6:12.2f}   {7:12.2e}".format(
            "c33", g.bp[n]['ec_gpa'][2,2], g.potfit['bp_best'][n]['ec_gpa'][2,2], g.potfit['bp_current'][n]['ec_gpa'][2,2],
            "c44", g.bp[n]['ec_gpa'][3,3], g.potfit['bp_best'][n]['ec_gpa'][3,3], g.potfit['bp_current'][n]['ec_gpa'][3,3]
           ))
      print(" {0:5s}  {1:12.2f}   {2:12.2f}   {3:12.2e}      {4:5s}  {5:12.2f}   {6:12.2f}   {7:12.2e}".format(
            "c55", g.bp[n]['ec_gpa'][4,4], g.potfit['bp_best'][n]['ec_gpa'][4,4], g.potfit['bp_current'][n]['ec_gpa'][4,4],
            "c66", g.bp[n]['ec_gpa'][5,5], g.potfit['bp_best'][n]['ec_gpa'][5,5], g.potfit['bp_current'][n]['ec_gpa'][5,5]
           ))
      print(" {0:5s}  {1:12.2f}   {2:12.2f}   {3:12.2e}      {4:5s}  {5:12.2f}   {6:12.2f}   {7:12.2e}".format(
            "c12", g.bp[n]['ec_gpa'][0,1], g.potfit['bp_best'][n]['ec_gpa'][0,1], g.potfit['bp_current'][n]['ec_gpa'][0,1],
            "c13", g.bp[n]['ec_gpa'][0,2], g.potfit['bp_best'][n]['ec_gpa'][0,2], g.potfit['bp_current'][n]['ec_gpa'][0,2]
           ))
      print(" {0:5s}  {1:12.2f}   {2:12.2f}   {3:12.2e}      ".format(
            "c23", g.bp[n]['ec_gpa'][1,2], g.potfit['bp_best'][n]['ec_gpa'][1,2], g.potfit['bp_current'][n]['ec_gpa'][1,2]
           ))

      display.hr() #------------------------
      print("RSS Breakdown")
      a0 = g.bp[n]['rss_details']['a0']
      e0 = g.bp[n]['rss_details']['e0']
      b0 = g.bp[n]['rss_details']['b0']
      ec = g.bp[n]['rss_details']['ec']
      rose = g.bp[n]['rss_details']['rose']
      msp_shape = g.bp[n]['rss_details']['msp_shape']
      msp_values = g.bp[n]['rss_details']['msp_values']
      bm_eos_shape = g.bp[n]['rss_details']['bm_eos_shape']
      bm_eos_values = g.bp[n]['rss_details']['bm_eos_values']


      print(" {:22s} {:12.4e} {:12.4e} {:12.4e} {:12.4e} {:12.4e}".format(
            'a0, e0, b0, ec, rose: ', a0, e0, b0, ec, rose))
      print(" {:22s} {:12.4e} {:12.4e} {:12.4e} {:12.4e}".format(
            'mspv, msps, bmv, bms ', msp_values, msp_shape, bm_eos_values, bm_eos_shape))


      display.hrr() ######################
      print("Parameters")
      #for i in range(
      print(g.potfit['p_current'])
      display.hrr() ######################




  @staticmethod
  def show_2():
    if(g.potfit['improvement_counter'] == 0):
      print(" {:10d} {:18.6e} {:18.6e}".format(
             g.potfit['counter'], g.potfit['rss_best']))







  @staticmethod
  def hr():
    for i in range(display.w):
      print("-", end="")
    print()

  @staticmethod
  def hrr():
    for i in range(display.w):
      print("#", end="")
    print()

  @staticmethod
  def clear():
    os.system('cls' if os.name == 'nt' else 'clear') 
    
