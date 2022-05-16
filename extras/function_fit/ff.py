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
    print("Function Fit")
    os.makedirs('plots', exist_ok=True)
    if(len(sys.argv)<2):
      print("Please specifiy data file.")
      exit()
    if(len(sys.argv)<3):
      print("Please specifiy potential file.")
      exit()

    datafile = sys.argv[1]
    potfile = sys.argv[2] 
    potfileout = potfile.replace(".", "_fit.")
    potfilebn = os.path.basename(potfile)
    plotname = "plots/" + potfilebn + '.eps'
    plotname_zoom = "plots/" + potfilebn + '_zoom.eps'

    start_d = 0
    if(len(sys.argv) >= 4):
      start_d = int(sys.argv[3])
    
    if(not os.path.isfile(datafile)):
      print("Data file does not exist.")
      exit()
    if(not os.path.isfile(potfile)):
      print("Pot file does not exist.")
      exit()

    ff.d = data_file.read(datafile)
    ff.dfull = numpy.zeros((len(ff.d), 4,),dtype=numpy.float64,order='F')
    interp.fill(ff.d[:,0], ff.d[:,1], 5, ff.dfull)

    """
    xi = 0.5
    yi = interp.interpolate_row(xi, ff.dfull[:,:], 5)    
    print(yi)
    
    
    xi = 1.5
    yi = interp.interpolate_row(xi, ff.dfull[:,:], 5)    
    print(yi)
    
    
    xi = 2.5
    yi = interp.interpolate_row(xi, ff.dfull[:,:], 5)    
    print(yi)
    
    
    xi = 3.5
    yi = interp.interpolate_row(xi, ff.dfull[:,:], 5)    
    print(yi)
    
    
    
    xi = 5.5
    yi = interp.interpolate_row(xi, ff.dfull[:,:], 5)    
    print(yi)
    exit()
    """

    ff.d = ff.d [start_d:,:]
    pfunc = function_file.read(potfile)
    ff.name = pfunc['type']
    ff.p = pfunc['p']
    ff.pf = pfunc['pf']
    p = numpy.copy(ff.p)

    print("Start",ff.rss_opt(p))

    # Fit
    """
    p = sa.run(ff.rss_opt, p, tempiter=20, tmax=1000)
    print("SA",ff.rss_opt(p))

    res = basinhopping(ff.rss_opt, p, niter=50, T=1.0, stepsize=0.5)
    p = res['x']
    print(ff.rss_opt(p))
    print("BH",ff.rss_opt(p))

    res = minimize(ff.rss_opt, p, method='Nelder-Mead', options={'maxiter': 1000, 'disp': False})
    p = res['x']
    print("NM",ff.rss_opt(p))
    """

    p = sa.run(ff.rss_opt, p, tempiter=10, tmax=100)
    print("SA",ff.rss_opt(p))

    res = minimize(ff.rss_opt, p, method='bfgs', options={'gtol': 1e-10, 'maxiter': 1000, 'disp': False})
    p = res['x']
    print("BFGS",ff.rss_opt(p))
    

    y_fit = ff.f(ff.name, ff.d[:,0], p, ff.pf)
    plt.plot(ff.d[:,0], ff.d[:,1])
    plt.plot(ff.d[:,0], y_fit)
    #plt.ylim(-10, 20)
    plt.savefig(plotname)

    y_fit = ff.f(ff.name, ff.d[:,0], p, ff.pf)
    plt.plot(ff.d[:,0], ff.d[:,1])
    plt.plot(ff.d[:,0], y_fit)
    plt.ylim(-4, 10)
    plt.savefig(plotname_zoom)

    fh = open(potfileout, 'w')
    fh.write("#TYPE " + str(ff.name) + "\n")
    fh.write("#P ")
    for pn in p:
      fh.write(str(pn) + " ")    
    fh.write("\n")    
    fh.write("#PF ")
    for pfn in ff.pf:
      fh.write(str(pfn) + " ")    
    fh.write("\n")
      

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

class sa:
  
  def run(f, p, niter=1000, tempiter=20, tmax=100, tmin=0.01):
    
    p_best = numpy.copy(p)
    rss = f(p)
    rss_best = rss
    p_new = numpy.zeros(len(p))
    t = numpy.linspace(tmax, tmin, tempiter)
    pvar = 0.9
    v = 10.0    
    for n in range(tempiter):
      p = numpy.copy(p_best)
      for m in range(niter):
        p_new[:] = p[:] + pvar**n * v * (0.5 - numpy.random.rand(len(p)))
        rss_new = f(p_new)
        if(rss_new < rss or numpy.random.uniform() < numpy.exp((rss-rss_new) / t[n])):
            p = numpy.copy(p_new)
            rss = rss_new
            if(rss_new < rss_best):
              p_best = numpy.copy(p_new)
              rss_best = rss
    return p_best

class data_file:

  @staticmethod
  def read(file_name):
    d = []
    fh = open(file_name, 'r')
    for row in fh:
      row = row.strip()
      if(row != '' and row[0] != '#'):
        row = data_file.one_space(row)
        f = row.split(" ")        
        line = []
        for fn in f:
          line.append(float(fn))        
        d.append(line)
    d = numpy.asarray(d)
    return d

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

class function_file:
 
  @staticmethod
  def read(path):

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







