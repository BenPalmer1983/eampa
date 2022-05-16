#!/bin/python3
import numpy
"""
units
"""

class units:

  @staticmethod
  def convert(conv_from, conv_to, value):
    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

    if(type(value) == int or type(value) == str):
      try:
        value = float(value)
      except:
        return None
      return units.convert_inner(conv_from, conv_to, value)
    elif(type(value) == float or type(value) == numpy.float64):
      return units.convert_inner(conv_from, conv_to, value)

    elif(type(value) == numpy.ndarray):
      output = numpy.copy(value)
      if(value.ndim == 1):
        for i in range(value.shape[0]):
          output[i] = units.convert_inner(conv_from, conv_to, value[i])
      elif(value.ndim == 2):
        for i in range(value.shape[0]):
          for j in range(value.shape[1]):
            output[i,j] = units.convert_inner(conv_from, conv_to, value[i,j])


      return output

  @staticmethod
  def convert_inner(conv_from, conv_to, value_in):

    # LENGTH METERS
    length = {
    'MILES': 0.000621373, 
    'KM': 0.001,
    'M': 1.0,
    'FT': 3.28084,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }
    
    # AREA METERS SQUARED
    area = {
    'M2': 1.0,
    'CM2': 1E4,
    'MM2': 1E6,
    'UM2': 1E12,
    'NM2': 1E18,
    'ANG2': 1E20,
    }
    
    # VOLUME METERS CUBED
    volume = {
    'M2': 1.0,
    'CM2': 1E6,
    'MM2': 1E9,
    'UM2': 1E18,
    'NM2': 1E27,
    'ANG2': 1E30,
    }

    # ENERGY J
    energy = {
    'J': 1.0,
    'EV': 6.2415E18,
    'RY': 4.5874E17,
    }

    # FORCE N
    force = {
    'N': 1.0,
    'RY/BOHR': 2.4276e7,
    'EV/ANG':6.2414E8,    
    }
    
    # VELOCITY
    velocity = {
    'M/S': 1.0,
    'MPH': 2.25,    
    }
    
    # PRESSURE
    pressure = {
    'PA': 1.0,
    'GPA': 1.0E-9,    
    'BAR': 1.0E-5,    
    'ATMOSPHERE': 9.8692E-6,    
    'PSI': 1.45038E-4, 
    'KBAR': 1.0E-8,   
    'RY/BOHR3': 6.857E-14,   
    'EV/ANG3': 6.241E-12
    }
    
    # CHARGE DENSITY (UNIT CHARGE PER VOLUME - ANG^3)
    charge_density = {
    'ANG-3': 1.0,
    'BOHR-3': 0.14812,    
    }
    
    # TEMPERATURE
    
    
    unit_list = [length, area, volume, energy, force, velocity, pressure, charge_density]
    
    for l in unit_list:
      if(conv_from in l.keys() and conv_to in l.keys()):
        return round((l[conv_to] / l[conv_from]) * float(value_in),9)
  


  @staticmethod
  def test():
    print ("Units Testing")

    f = 'miles'
    t = 'cm'
    x = 2.5
    y = units.convert(conv_from=f, conv_to=t, value=x)
    print(x, f, " = ", y, t)

    aa = numpy.random.rand(10,5)
    bb = units.convert(conv_from=f, conv_to=t, value=aa)
    print(aa)
    print(bb)
    #print(type(aa))



if __name__ == "__main__":
  units.test()    

