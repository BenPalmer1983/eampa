
import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib.ticker import MultipleLocator

"""
python3 neighbours.py structure=FCC a0=3.42 title="FCC Iron: Neighbours to Atom at 0,0,0" out="fe_fcc.eps"

python3 neighbours.py structure=BCC a0=2.80 title="BCC Iron: Neighbours to Atom at 0,0,0" out="fe_bcc.eps"

"""

class neighbours:

  def run():
    s = {
         'structure': None,
         'a0': None,
         'size': None,
         'rcut': 10.0,
         'title': '',
         'out': 'neighbours.eps',
        }
    
    
    for i in range(1, len(sys.argv)):
      f = sys.argv[i].split("=")    
      s[f[0]] = f[1]

    if(s['structure'] == None or s['a0'] == None):
      print("Example Useage")
      print("python3 neighbours.py structure=FCC a0=4.04")
      print("python3 neighbours.py structure=FCC a0=4.04 title=\"FCC Aluminium: Neighbours to Atom at 0,0,0\" out=\"al.eps\"")
      exit()

    s['a0'] = float(s['a0'])
    if(s['size'] is None):
      xmax = s['rcut'] + 1.0
      s['size'] = numpy.ceil(s['rcut'] / s['a0']) + 1
    else:
      xmax = s['size'] -1 * s['a0']   
    s['size'] = int(s['size'])



    mx = 0.0
    my = 0.0
    mz = 0.0 
   
    d = {}

    for i in range(-s['size'], s['size']):
      for j in range(-s['size'], s['size']):
        for k in range(-s['size'], s['size']):
          a = neighbours.getcoords(s['structure'])
          for n in range(len(a)):
            c = [a[n][0] + i, a[n][1] + j, a[n][2] + k] 
            if(not (mx == c[0] and my == c[1] and mz == c[2])):
              r = s['a0'] * numpy.sqrt((mx - c[0])**2 + (my - c[1])**2 + (mz - c[2])**2)
              if(r<= s['rcut']):
                if(r not in d.keys()):
                  d[r] = 0
                d[r] = d[r] + 1
    
    d_arr = numpy.zeros((len(d.keys()),2),)
    n = 0

    for k in d.keys():
      d_arr[n,0] = k
      d_arr[n,1] = d[k]
      n = n + 1

    mlx = MultipleLocator(0.1)
    majlx = MultipleLocator(1)
    mly = MultipleLocator(1)
    plt.clf()
    plt.figure(figsize=(12,8))    
    plt.title(s['title'])
    plt.xlabel('Separation/Angs')
    plt.ylabel('Neighbour Count')
    plt.xlim(0.0, xmax)
    plt.ylim(0.0, 1.05 * max(d_arr[:,1]))
    plt.stem(d_arr[:,0], d_arr[:,1])
    plt.axes().xaxis.set_major_locator(majlx)
    plt.axes().xaxis.set_minor_locator(mlx)
    plt.axes().yaxis.set_minor_locator(mly)
    plt.grid(color='grey', linestyle='-', linewidth=1)
    plt.savefig(s['out'], format='eps')
    plt.close('all') 


  def getcoords(structure):
    if(structure.upper() == "FCC"):
      return [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]]
    elif(structure.upper() == "BCC"):
      return [[0.0,0.0,0.0],[0.5,0.5,0.5]]
    elif(structure.upper() == "SC"):
      return [[0.0,0.0,0.0],[0.5,0.5,0.5]]



neighbours.run()


"""
python3 neighbours.py structure=FCC a0=4.04 title="FCC Aluminium: Neighbours to Atom at 0,0,0" out="al.eps"


"""
