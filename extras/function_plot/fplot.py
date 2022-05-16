import sys
import numpy
import matplotlib.pyplot as plt
import os
from f_toolbox import fnc
from f_toolbox import interp
from scipy.optimize import minimize
from scipy.optimize import basinhopping



class ff:

  d = None
  dfull = None
  name = None
  pf = None
  

  def main():
    print("Function Plot")

    if(len(sys.argv)<2):
      print("Please specifiy potential file.")
      exit()

    os.makedirs('plots', exist_ok=True)
    #os.makedirs('pots', exist_ok=True)
    
    potfile = sys.argv[1] 
    potfilebn = os.path.basename(potfile)
    plotname = "plots/" + potfilebn + '.eps'
    plotname_zoom1 = "plots/" + potfilebn + '_zoom1.eps'
    plotname_zoom2 = "plots/" + potfilebn + '_zoom2.eps'

    start_d = 0
    if(len(sys.argv) >= 3):
      start_d = int(sys.argv[2])
    
    if(not os.path.isfile(potfile)):
      print("Pot file does not exist.")
      exit()

    pfunc = function_file.read(potfile)


    if(pfunc['tab'] is None):
      ff.name = pfunc['type']
      ff.p = pfunc['p']
      ff.pf = pfunc['pf']
    
      x = numpy.linspace(0.0, 10.0, 10001)
      y = ff.f(ff.name, x[:], ff.p, ff.pf)
    else:
      x = pfunc['tab'][:,0]
      y = pfunc['tab'][:,1]

    fig, ax = plt.subplots()
    ax.plot(x[:], y[:])
    plt.savefig(plotname)


    fig, ax = plt.subplots()
    ax.plot(x[:], y[:])
    ax.set_ylim(-4, 10)
    plt.savefig(plotname_zoom1)

    fig, ax = plt.subplots()
    ax.plot(x[:], y[:])
    ax.set_ylim(-4, 1000)
    plt.savefig(plotname_zoom2)

      

  @staticmethod
  def f(fname, x, p, pf):
    return fnc.fv(fname, x, p, pf)

  @staticmethod
  def fv(p):
    return fnc.fv(ff.name, ff.d[:,0], p, ff.pf)

  @staticmethod
  def rss_opt(p):
    y = fnc.fv(ff.name, ff.d[:,0], p, ff.pf)
    return sum((y[:] - ff.d[:,1])**2)


class function_file:
 
  @staticmethod
  def read(path):

    fh = open(path, "r")
    for line in fh:
      line = function_file.one_space(line.strip())
      f = line.split(" ")
      if(f[0].upper() == "#TYPE"):
        fname = f[1].lower()
      break
    fh.close()

    if(fname == 'tab'):
      pfunc = function_file.read_tabulated(path)
    else:
      pfunc = function_file.read_analytic(path)
    
    return pfunc

  
  @staticmethod
  def read_analytic(path):

    # Read file in
    raw = []
    fh = open(path, "r")
    for line in fh:
      raw.append(line.strip())
    fh.close()

    # Remove comments
    for n in range(len(raw)):
      t = raw[n].split("//")
      raw[n] = t[0]
      t = raw[n].split("!")
      raw[n] = t[0]


    temp = []
    # Read it in line by line
    cont = False
    for n in range(len(raw)):
      cl = False
      line = raw[n].strip()
      if(len(line) > 0 and line[-1] == "&"):
        cl = True 
        line = line[:-1].strip() + " "
      if(len(temp) > 0 and cont):
        temp[-1] = temp[-1] + line
        cont = False
      else:  
        temp.append(line)
      if(cl):
        cont = True
    
    file_data = []
    for line in temp:
      if(line.strip() != ""):
        line = line.replace(",", " ")
        file_data.append(function_file.one_space(line).strip())

    pfunc = {
              'type': None,
              'p': None,
              'pf': None,
              'pl': None,
              'pu': None,
              'vr': None,
              'tab': None,
            }

    for line in file_data:   
      f = line.split(" ")  
      if(f[0].upper() == "#TYPE"):
        pfunc['type'] = f[1]
      elif(f[0].upper() == "#P"):
        pfunc['p'] = numpy.asarray(f[1:], dtype=numpy.float64)
      elif(f[0].upper() == "#PF"):
        pfunc['pf'] = numpy.asarray(f[1:], dtype=numpy.float64)
      elif(f[0].upper() == "#PL"):
        pfunc['pl'] = numpy.asarray(f[1:], dtype=numpy.float64)
      elif(f[0].upper() == "#PU"):
        pfunc['pu'] = numpy.asarray(f[1:], dtype=numpy.float64)
      elif(f[0].upper() == "#VR"):
        pfunc['vr'] = float(f[1])
    return pfunc


  @staticmethod
  def read_tabulated(path):  

    pfunc = {
              'type': None,
              'p': None,
              'pf': None,
              'pl': None,
              'pu': None,
              'vr': None,
              'tab': None,
            }

    
    # Read file in
    raw = []
    fh = open(path, "r")
    for line in fh:
      raw.append(line.strip())
    fh.close()

    d = []
    for line in raw[1:]:
      line = function_file.one_space(line.replace("\t", " ").replace(",", " ")).strip()
      if(line != ""):
        d.append(line.split(" "))
    pfunc['tab'] = numpy.asarray(d, dtype=numpy.float64)
    return pfunc

  @staticmethod
  def one_space(line, sep=" "):
    out = ''   
    indata = 0
    last_char = None
    for char in line:
      if(indata == 1 and char != "'" and last_char != "\\"):
        out = out + char
      elif(indata == 1 and char == "'" and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 2 and char != '"' and last_char != "\\"):
        out = out + char
      elif(indata == 2 and char == '"' and last_char != "\\"):
        out = out + char
        indata = 0
      elif(indata == 0 and not (char == " " and last_char == " ")):
        out = out + char
      last_char = char
    return out 

if __name__ == "__main__":
  ff.main()    







