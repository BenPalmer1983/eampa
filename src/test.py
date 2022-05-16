#!/bin/python3
"""
test
"""

import numpy
import os
import time

from f_toolbox import math

from g import g
from ds import ds
from std import std
from output import output
from label import label
from nl import nl
from configs import configs
from relax import relax
from heat import heat
from rss import rss
from bp import bp
from md import md
from potfit import potfit


class test:


  @staticmethod
  def run():
    
    output.log("Test")
    # test.common_configs()
    # test.change_nl_a0()
    #test.change_nl_uv_1()  
    #test.change_nl_uv_2()  
    #test.change_nl_uva0()
    #test.set_of_configs()
    #test.test_stress()
    #test.heat()
    #test.relax_1()
    #test.relax_2()
    #test.relax_3()
    #test.relax_4()
    #test.relax_5()
    #test.relax_6()
    #test.update_nl()
    #test.calc_rss()
    #test.failing()
    #test.calc_bp()
    #test.verlet()
    #test.potfit()
    #test.relax_7()
    #test.relax_8()
    #test.relax_9()
    #test.duplicateconfig()
    #test.copyconfiga()
    #test.copyconfigb()
    #test.calc_bp2()
    #test.calc_bp3()
    #test.calc_bp4()
    #test.openmp()
    test.strains()





  @staticmethod
  def strains():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    print(config_id)
    relax.run_relax(config_id, pmin='e', pvar='a0')
    print(g.configs[config_id]['c_energy'])
    print(g.configs[config_id]['a0'])
    print(g.configs[config_id]['uv'])
    print(g.configs[config_id]['volume'])
    relaxed_id = config_id

    
    a = numpy.linspace(0.0, 0.01, 11)
    e = numpy.zeros((11,), dtype=numpy.float64)

    for n in range(len(a)):
      config_id = configs.duplicate(relaxed_id)
      configs.strain(config_id, 'ctd1', a[n])
      configs.calc(config_id)
      e[n] = g.configs[config_id]['c_energy']
      print(a[n], e[n])
    


    for n in range(len(a)):
      config_id = configs.duplicate(relaxed_id)
      configs.strain(config_id, 'd1', a[n])
      configs.calc(config_id)
      e[n] = g.configs[config_id]['c_energy']
      print(a[n], e[n])


  @staticmethod
  def openmp():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.00], [2], 7.5, 'test', 'e')
    a = time.time()
    for n in range(1000):
      configs.calc(config_id)
    print(time.time() - a)
    print(g.configs[config_id]['c_energy'])



  @staticmethod
  def calc_bp4():
    print("calc_bp")
    for i in range(10):
      rss = bp.run_rss()
      print(rss)


  @staticmethod
  def calc_bp3():
    print("calc_bp")
    bp.run()


  @staticmethod
  def calc_bp2():
    print("calc_bp")
    bp.build()
    configs.index()
    bp.reset()
    configs.index()
    #configs.index()
    bp.calc()
    total_rss = bp.rss()
    print(total_rss) 
    bp.calc()
    total_rss = bp.rss()
    print(total_rss) 
    bp.reset()
    bp.calc()
    total_rss = bp.rss()
    print(total_rss) 
    bp.reset()
    bp.calc()
    total_rss = bp.rss()
    print(total_rss) 



  @staticmethod
  def copyconfigb():
    config_a = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    config_b = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    config_c = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    configs.calc(config_a)
    configs.calc(config_b)
    configs.calc(config_c)
    configs.index()
    heat.gauss(config_b, 0.001)
    configs.calc(config_b)

    configs.index()

    configs.duplicate(config_b, config_c)
    configs.index()

  @staticmethod
  def copyconfiga():
    config_a = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    config_b = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 6.5, 'test', 'e')
    config_c = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')

    configs.index()
    configs.duplicate(config_a, config_b)
    configs.calc(config_a)
    configs.calc(config_b)
    configs.index()



  @staticmethod
  def duplicateconfig():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    configs.calc(config_id)
    print(g.configs[config_id]['c_energy'])
    configs.index()
    config_id = configs.duplicate(config_id)
    configs.index()
    heat.gauss(config_id, 0.001)
    configs.calc(config_id)
    configs.index()
    config_id = configs.duplicate(config_id)
    configs.calc(config_id)
    configs.index()


  @staticmethod
  def relax_9():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    relax.run_relax(config_id, pmin='e', pvar='a0')
    

  @staticmethod
  def relax_8():
    config_id = 0
    configs.calc(config_id)
    print(g.configs[config_id]['c_energy'])
    relax.run_relax(config_id, pmin='e', pvar='coords')
    print(g.configs[config_id]['c_energy'])
    print(g.configs[config_id]['coords'])


  @staticmethod
  def relax_7():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    heat.gauss(config_id, 0.001)
    configs.calc(config_id)
    print(g.configs[config_id]['c_energy'])
    relax.run_relax(config_id, pmin='e', pvar='coords')
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.00], [4], 7.5, 'test', 'e')
    configs.calc(config_id)
    print(g.configs[config_id]['c_energy'])


  @staticmethod
  def potfit():
    print("Testing potfit")
    potfit.run()


  @staticmethod
  def verlet():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [0.99], [4], 7.5, 'test', 'e')
    heat.gauss(config_id, 0.001)
    configs.calc(config_id)
    md.run(config_id)


  @staticmethod
  def calc_bp():
    print("calc_bp")
    bp.build()
    bp.calc()


  @staticmethod
  def failing():
    # al_100_SCF0000.out
    # ec_6_005.out
    nl.build('efs')
    configs.calc('efs')
    print(g.configs)


  @staticmethod
  def calc_rss():
    output.log("Run test RSS", verbose=0)
    nl.build('efs')
    configs.calc('efs')
    print(g.rss_w)
    total_rss = configs.rss('efs')
    print(total_rss)
    

  @staticmethod
  def update_nl():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.calc(config_id)
    nl.rebuild_nl(config_id)
    nl.update_coords(config_id)
    configs.calc(config_id)

    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'efs')
    heat.gauss(config_id, 0.005)
    nl.rebuild_nl(config_id)
    configs.calc(config_id)
    print(g.configs[config_id]['c_forces'][-10:,:])
    nl.update_coords(config_id)
    configs.calc(config_id)
    print(g.configs[config_id]['c_forces'][-10:,:])
 

  @staticmethod
  def relax_1():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    relax.run_relax(config_id, pmin='e', pvar='a0')
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    relax.run_relax(config_id, pmin='f', pvar='a0')


  @staticmethod
  def relax_2():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [0.99], [4], 7.5, 'test', 'e')
    relax.run_relax(config_id, pmin='e', pvar='uv')


  @staticmethod
  def relax_3():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [0.99], [4], 7.5, 'test', 'e')
    relax.run_relax(config_id, pmin='e', pvar='uva0')


  @staticmethod
  def relax_4():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [0.99], [4], 7.5, 'test', 'e')
    heat.gauss(config_id, 0.005)
    configs.calc(config_id)
    relax.run_relax(config_id, pmin='e', pvar='coords')


  @staticmethod
  def relax_5():
    # Not working
    config_id = configs.add_common('fcc', ['Al'], 4.018, [1.00], [4], 7.5, 'test', 'e')
    heat.gauss(config_id, 0.005)
    configs.calc(config_id)
    relax.run_relax(config_id, pmin='ef', pvar='coords')


  @staticmethod
  def relax_6():
    config_id = configs.add_common('fcc', ['Al'], 4.018, [1.00], [4], 7.5, 'test', 'e')
    


  @staticmethod
  def heat():
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.calc(config_id)
    heat.gauss(config_id, 0.01)
    configs.calc(config_id)
    heat.gauss(config_id, 0.01)
    configs.calc(config_id)
 


  @staticmethod
  def set_of_configs():
    output.log("Configs 6.5 cutoff", verbose=0)
    configs.add_common('fcc', ['Al'], 4.00, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.01, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.01013367, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.02, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.03, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 6.5, 'test', 'e')
    configs.calc('all')

    output.log("Configs 7.5 cutoff")
    configs.add_common('fcc', ['Al'], 4.00, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.01, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.01013367, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.02, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.03, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.calc('all')

    output.log("Configs 8.5 cutoff")
    configs.add_common('fcc', ['Al'], 4.00, [1.0], [4], 8.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.01, [1.0], [4], 8.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.01013367, [1.0], [4], 8.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.02, [1.0], [4], 8.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.03, [1.0], [4], 8.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 8.5, 'test', 'e')
    configs.calc('all')


  @staticmethod
  def common_configs():
    output.log("Test Common Configs")
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.02], [4], 6.5, 'test', 'ef')
    configs.add_common('fcc', ['Al'], 4.04, [1.04], [4], 6.5, 'test', 'efs')
    configs.add_common('fcc', ['Al'], 4.05, [1.02], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.05, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.06, [1.0], [4], 6.5, 'test', 'e')
    configs.add_common('bcc', ['Al'], 4.06, [1.0], [4], 6.5, 'test', 'e')



  @staticmethod
  def change_nl_a0():
    output.log("Build NL Configs")
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.02], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.04], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.05, [1.02], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.038, [0.995], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.05, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.06, [1.0], [4], 7.5, 'test', 'e')

    output.log("")
    output.log("Change a0 Configs")
    configs.calc('all')
    nl.change_cell(1, 4.05, [1.0])
    configs.calc(1)
    nl.change_cell(1, 4.06, [1.0])
    configs.calc(1)
    nl.change_cell(1, 4.04, [1.0])
    configs.calc(1)
    




  @staticmethod
  def change_nl_uv_1():
    output.log("Build NL Configs")
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.02], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.04], [4], 7.5, 'test', 'e')

    output.log("")
    output.log("Change UV Configs")
    configs.calc('all')
    nl.change_cell(1, 4.04, [1.0])
    configs.calc(1)
    nl.change_cell(1, 4.04, [1.02])
    configs.calc(1)
    nl.change_cell(1, 4.04, [1.04])
    configs.calc(1)
    nl.change_cell(1, 4.04, [1.0])
    configs.calc(1)


  @staticmethod
  def change_nl_uv_2():
    output.log("Build NL Configs")
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.calc('all')
    print(config_id)
    nl.change_cell(config_id, 4.04, [1.0])
    configs.calc(config_id)
    nl.change_cell(config_id, 4.04, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    configs.calc(config_id)
    nl.change_cell(config_id, 4.04, [[1.0, 0.01, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    configs.calc(config_id)
    config_id = configs.add_common('fcc', ['Al'], 4.04, [[1.0, 0.01, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [4], 7.5, 'test', 'e')
    configs.calc(config_id)



  @staticmethod
  def change_nl_uva0():
    output.log("Build NL Configs")
    configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.03, [1.02], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.04, [1.04], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.05, [1.02], [4], 7.5, 'test', 'e')
    configs.add_common('fcc', ['Al'], 4.038, [0.995], [4], 7.5, 'test', 'e')

    output.log("")
    output.log("Change a0 and UV Configs")
    configs.calc('all')
    nl.change_cell(1, 4.03, [1.02])
    configs.calc(1)
    nl.change_cell(1, 4.04, [1.04])
    configs.calc(1)
    nl.change_cell(1, 4.05, [1.02])
    configs.calc(1)
    nl.change_cell(1, 4.038, [0.995])
    configs.calc(1)
    nl.change_cell(1, 4.04, [1.0])
    configs.calc(1)

    configs.output('all')
    #configs.calc_efs()



  @staticmethod
  def test_stress():
    output.log("Build NL Configs")
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.0], [4], 7.5, 'test', 'efs')
    configs.calc(config_id)
    print(g.configs[config_id]['c_stress'])
    config_id = configs.add_common('fcc', ['Al'], 4.04, [1.01, 1.02, 1.012,0.01,0.0,0.0], [4], 7.5, 'test', 'efs')
    configs.calc(config_id)
    print(g.configs[config_id]['c_stress'])
    



  @staticmethod
  def minverse():

    uv = numpy.zeros((3,3,), dtype=numpy.float64)
    uv[0,0] = 4.0
    uv[0,1] = 1.0
    uv[1,0] = 15
    uv[1,1] = 4.0
    uv[2,2] = 4.0
    uv[1,2] = -0.4
    uv[2,1] = 0.4

    s = time.time()
    for i in range(10000):
      uvinv = numpy.linalg.inv(uv)
    print(uvinv)
    print(time.time() - s)
    print()

    s = time.time()
    for i in range(10000):
      uvinva = math.minverse(uv)
    print(uvinva)
    print(time.time() - s)
    print()

    s = time.time()
    for i in range(10000):
      uvinva = math.minverse3(uv)
    print(uvinva)
    print(time.time() - s)
    print()



  @staticmethod
  def vecvol():
    uv = numpy.zeros((3,3,), dtype=numpy.float64)

    uv[0,0] = 4.0
    uv[0,1] = 1.0
    uv[1,0] = 1.0
    uv[1,1] = 4.0
    uv[2,2] = 4.0
    uv[1,2] = 0.4
    uv[2,1] = 0.4

    vol = math.vecvol(uv)
    print(uv)
    print(vol)
    
    a = uv[0,:]
    b = uv[1,:]
    c = uv[2,:]
    vol = numpy.dot(c, numpy.cross(a, b))
    print(vol)
    