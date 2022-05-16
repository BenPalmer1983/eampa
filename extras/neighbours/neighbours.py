import numpy
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import sys

class neighbours:

  configs = []
  rcut = 6.5
  bins = 200
  din = None
  dout = None

  def main():
    if(len(sys.argv)<2):
      print("Please specifiy input db directory.")
      exit()
    if(len(sys.argv)<3):
      print("Please specifiy output directory.")
      exit()

    din = sys.argv[1]
    dout = sys.argv[2]

    if(not os.path.isdir(din)):
      print("Input directory does not exist.")
      exit()
 
    os.makedirs(dout, exist_ok=True)
    neighbours.din = din
    neighbours.dout = dout

    files = neighbours.files(din)
    for f in files:
      neighbours.read_config(f)
    neighbours.make()
    neighbours.tally()
   

  def files(thepath, file_arr = []):
    files = os.listdir(thepath)
    for f in files:
      filepath = thepath + "/" + f
      if(os.path.isdir(filepath)):    
        file_arr = neighbours.files(filepath, file_arr)
      elif(os.path.isfile(filepath)):
        file_arr.append(filepath)
    return file_arr
  

  def read_config(f):
    uv = numpy.zeros((3,3),)
    a0 = 0.0
    coords_crystal = []
    coords = []

    fh = open(f, 'r')
    for line in fh:
      line = line.replace("\t", " ")
      line = line.strip()
      line = line.split("//")
      line = neighbours.one_space(line[0])
      f = line.split(" ")
      if(f[0] == "/*"):
        pass
      elif(f[0].upper() == "#ALAT"):
        a0 = float(f[1])
      elif(f[0].upper() == "#X"):
        uv[0,0] = float(f[1])
        uv[0,1] = float(f[2])
        uv[0,2] = float(f[3])
      elif(f[0].upper() == "#Y"):
        uv[1,0] = float(f[1])
        uv[1,1] = float(f[2])
        uv[1,2] = float(f[3])
      elif(f[0].upper() == "#Z"):
        uv[2,0] = float(f[1])
        uv[2,1] = float(f[2])
        uv[2,2] = float(f[3])
      elif(f[0][0] != "#" and (len(f) == 4 or len(f) == 7)):
        coords_crystal.append([f[0], float(f[1]), float(f[2]), float(f[3])])
    fh.close()

    tr = a0 * uv
    neighbours.configs.append({'a0': a0, 'uv': uv, 'tr': tr, 'coords_crystal': coords_crystal, 'coords': None, 'gcoords': None, 'nl': None,})


  def make():
    for cn in range(len(neighbours.configs)):
      neighbours.make_coords(cn)
      neighbours.make_gcoords(cn)
    for cn in range(len(neighbours.configs)):
      neighbours.make_nl(cn)

  def make_coords(cn):
    neighbours.configs[cn]['coords'] = []
    tr = neighbours.configs[cn]['tr']
    for n in range(len(neighbours.configs[cn]['coords_crystal'])):
      label = neighbours.configs[cn]['coords_crystal'][n][0]
      v = numpy.asarray(neighbours.configs[cn]['coords_crystal'][n][1:4])
      v = numpy.matmul(tr, v)
      neighbours.configs[cn]['coords'].append([n, label, v])

  def make_gcoords(cn):
    neighbours.configs[cn]['gcoords'] = []
    tr = neighbours.configs[cn]['tr']
    for i in range(-1, 2):
      for j in range(-1, 2):
        for k in range(-1, 2):
          for n in range(len(neighbours.configs[cn]['coords_crystal'])):
            label = neighbours.configs[cn]['coords_crystal'][n][0]
            v = numpy.asarray(neighbours.configs[cn]['coords_crystal'][n][1:4])
            v[0] = v[0] + i
            v[1] = v[1] + j
            v[2] = v[2] + k
            v = numpy.matmul(tr, v)
            neighbours.configs[cn]['gcoords'].append([n, label, v])

  def make_nl(cn):
    rcutsq = neighbours.rcut ** 2
    c = neighbours.configs[cn]['coords']
    g = neighbours.configs[cn]['gcoords']
    neighbours.configs[cn]['nl'] = []
    for i in range(len(c)):
      for j in range(len(g)):
        if(c[i][0] < g[j][0]): 
          rdsq = sum((c[i][2][:] - g[j][2][:])**2)
          if(rdsq <= rcutsq):
            rd = numpy.sqrt(rdsq)
            neighbours.configs[cn]['nl'].append([c[i][0], c[i][1], g[j][0], g[j][1], rd])
  
  def tally():
    x = numpy.linspace(0.0, neighbours.rcut, neighbours.bins+1)
    neighbours.td = numpy.zeros((neighbours.bins,),dtype=numpy.int32)
    nl_plot = []
    for cn in range(len(neighbours.configs)):
      nl = neighbours.configs[cn]['nl']
      for n in range(len(nl)):
        nl_plot.append(nl[n][4])
        #bin = int(numpy.floor((nl[n][4] / neighbours.rcut) * neighbours.bins))
        #neighbours.td[bin] = neighbours.td[bin] + 1
    mlx = MultipleLocator(0.1)
    majlx = MultipleLocator(1)
    mly = MultipleLocator(1)
    plt.axes().xaxis.set_major_locator(majlx)
    plt.axes().xaxis.set_minor_locator(mlx)
    plt.axes().yaxis.set_minor_locator(mly)
    plt.grid(color='grey', linestyle='-', linewidth=1)
    count, bins, ignored = plt.hist(nl_plot, neighbours.bins, density=True)
    plt.savefig(neighbours.dout + '/' + neighbours.din + '_neighbours.eps')




  
  def one_space(inp):
    out = ''
    last = None
    for c in inp:
      if(not (last == " " and c == " ")):
        out = out + c
      last = c
    return out

neighbours.main()
