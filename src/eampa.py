#!/bin/python3

import os
import sys
import numpy
import time
import shutil

from g import g
from display import display
from ifile import ifile
from cfile import cfile
from nl import nl
from label import label
from output import output
from potential import potential
from configs import configs
from heat import heat


from efs import efs
from bp import bp
from rss import rss
from relax import relax
from test import test
from potfit import potfit


class eampa:


  @staticmethod
  def main():
    display.clear()
    g.start = time.time()
    eampa.set_dirs()
    g.fh = open(os.path.join(g.dir['out'],'log.txt'), 'w')
    output.log("EAMPA")

    eampa.load()
    eampa.run()

    g.fh.close()





  @staticmethod
  def set_dirs():
    g.dir['out'] = os.path.join(os.getcwd(), "out")
    if(os.path.isdir(g.dir['out'])):
      shutil.rmtree(g.dir['out'])

    g.dir['efs'] = os.path.join(g.dir['out'], "efs")
    g.dir['bp'] = os.path.join(g.dir['out'], "bp")
    g.dir['fit'] = os.path.join(g.dir['out'], "fit")
    g.dir['plots'] = os.path.join(g.dir['out'], "plots")
    g.dir['plots_potential'] = os.path.join(g.dir['plots'], "potential")
    g.dir['data'] = os.path.join(g.dir['out'], "data")
    g.dir['data_potential'] = os.path.join(g.dir['data'], "potential")
    g.dir['data_configs'] = os.path.join(g.dir['data'], "configs")

    for k in g.dir.keys():
      os.makedirs(g.dir[k], exist_ok=True)

  @staticmethod
  def load():

    if(len(sys.argv)<2):
      output.log("Please specifiy input file.", verbose=0)

    if(not os.path.isfile(sys.argv[1])):
      output.log("Input file does not exist.", verbose=0)
      exit()

    # Read Input File
    ifile.load(sys.argv[1])

    # Load any config files
    cfile.read_configs('efs')

    # Load one line configs
    configs.one_line_configs()

    # Load potential files
    potential.load()

    # Load into atom module
    potential.set()
    potential.set()
    potential.set()
    potential.summary_fortran()

  """
  Run types - choose which type of run to perform
  """

  @staticmethod
  def run():
    if(g.input['runtype'] == 'rss'):
      rss.run()
    elif(g.input['runtype'] == 'efs'):
      efs.run()
    elif(g.input['runtype'] == 'bp'):
      bp.run()
    elif(g.input['runtype'] == 'fit'):
      potfit.run()
    elif(g.input['runtype'] == 'relax'):
      output.log("Run Relax", verbose=0)
      relax.run()
    elif(g.input['runtype'] == 'test'):
      output.log("Run Test", verbose=0)
      test.run()




  @staticmethod
  def quit():
    g.end = time.time()
    output.log("Quitting", verbose=0)
    g.fh.close()
    exit()


if __name__ == "__main__":
  eampa.main()    

