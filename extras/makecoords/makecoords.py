import os
import numpy
import sys

"""
python3 makecoords.py Fe1,Fe2 fcc 4.04 1,1,1,0,0,0 4,4,4 fe_bulk

"""


class makecoords:

  @staticmethod
  def make():
    print("Example usage:")
    print("""python3 makecoords.py 1,2 fcc 4.04 1,0,0,0,1,0,0,0,1 4,4,4 alexample
python3 makecoords.py 1,2 fcc 4.04 cubic 4,4,4 alexample
python3 makecoords.py 1,2 c-hcp 4.04 hcp 4,4,4 alexample
python3 makecoords.py 1 hcp 2.718747 hcp 8,8,8 ruhcp
""")
    
    if(len(sys.argv) < 7):
      print("Check arguments.")
      exit()
    
    try:
      labels_prim = sys.argv[1].split(",")
      structure = sys.argv[2].lower()
      a0 = float(sys.argv[3])
      uv = sys.argv[4].split(",")
      c = sys.argv[5].split(",")
      dir_out = sys.argv[6]
      #shapes = sys.argv[7:]
    except:
      print("Check input arguments")
      exit()

        
    os.makedirs(dir_out, exist_ok=True)

    uv = makecoords.unit_vector(uv)
    c = numpy.asarray(c, dtype=numpy.int)  
    cm = numpy.zeros((3,3,), dtype=numpy.float64)
    cm[0,0] = c[0]
    cm[1,1] = c[1]
    cm[2,2] = c[2]
    uv = numpy.matmul(cm,uv) 
    box = a0 * uv
    box_inv = numpy.linalg.inv(box)
    
    print("Labels:     ", labels_prim)
    print("Structure:  ", structure)
    print("a0:         ", a0)
    print("uv:         ", uv[0,:])
    print("            ", uv[1,:])
    print("            ", uv[2,:])
    print("c:          ", c)
    print("box:        ", box[0,:])
    print("            ", box[1,:])
    print("            ", box[2,:])
    print("box inv:    ", box_inv[0,:])
    print("            ", box_inv[1,:])
    print("            ", box_inv[2,:])
    
    coords_prim = makecoords.structure(structure)
    labels = []
    coords = numpy.zeros((len(coords_prim) * c[0] * c[1] * c[2], 3,), dtype=numpy.float64)
    
    # Make coords
    m = 0
    for i in range(c[0]):
      for j in range(c[1]):
        for k in range(c[2]):
          for n in range(len(coords_prim)):
            labels.append(labels_prim[n % len(labels_prim)])
            coords[m,:] = numpy.asarray([(i+coords_prim[n,0])/c[0], (j+coords_prim[n,1])/c[1], (k+coords_prim[n,2])/c[2]])
            m = m + 1

    # Apply shape to crystal coords

            
    # Real coords
    coords_real = numpy.zeros((len(coords_prim) * c[0] * c[1] * c[2], 3,), dtype=numpy.float64)
    for n in range(len(coords)):
      coords_real[n,:] = numpy.matmul(box, coords[n,:])
    

    # Apply shape to real coords

    #a0

    # Output
    makecoords.output_coords(labels, coords, dir_out, 'coords_crystal.txt')
    makecoords.output_coords(labels, coords_real, dir_out, 'coords_real.txt')
    makecoords.output_coords_gb(labels, a0, c[0], c[1], c[2], coords, dir_out, 'coords_real_gb.txt')


  #@staticmethod
  #def structure(structure):
  #  if(structure=='fcc'):
  #    return numpy.asarray([[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]])
  #    #return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
  @staticmethod
  def structure(structure):
    if(structure=='fcc'):
      return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])
    elif(structure=='bcc'):
      return numpy.asarray([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    elif(structure=='sc'):
      return numpy.asarray([[0.0, 0.0, 0.0]])
    elif(structure=='hcp'):
      return numpy.asarray([[0.0,0.0,0.0],[0.3333333,0.6666667,0.5000000]])
    elif(structure=='c-fcc'):
      return numpy.asarray([[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]])
    elif(structure=='c-bcc'):
      return numpy.asarray([[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]])
    elif(structure=='c-sc'):
      return numpy.asarray([[0.5, 0.5, 0.5]])
    elif(structure=='c-hcp'):
      return numpy.asarray([[0.3333333,0.1666667,0.2500000],[0.6666667,0.8333333,0.7500000]])

  @staticmethod
  def unit_vector(uv_in):
    uv = numpy.zeros((3,3,), dtype=numpy.float64)
    if(isinstance(uv_in, str)):
      pass
    elif(isinstance(uv_in, list)):
      if(len(uv_in) == 1):
        if(uv_in[0] == 'cubic'):
          uv[0,0] = 1
          uv[1,1] = 1
          uv[2,2] = 1
        elif(uv_in[0] == 'hcp'):
          uv[0,0] = 1
          uv[1,0] = -0.5
          uv[1,1] = 0.866025
          uv[2,2] = 1.582300
      elif(len(uv_in) == 6):
        uv[0,0] = uv_in[0]
        uv[1,1] = uv_in[1]
        uv[2,2] = uv_in[2]
        uv[0,1] = uv_in[5]
        uv[1,0] = uv_in[5]
        uv[0,2] = uv_in[4]
        uv[2,0] = uv_in[4]
        uv[1,2] = uv_in[3]
        uv[2,1] = uv_in[3]
      elif(len(uv_in) == 9):
        uv[0,:] = uv_in[0:3]
        uv[1,:] = uv_in[3:6]
        uv[2,:] = uv_in[6:9]
    return uv
    

  @staticmethod
  def output_coords(labels, coords, path, filename):
    fh = open(os.path.join(path, filename),'w')
    for n in range(len(labels)):
      fh.write("{:10s}{:12.7f}{:12.7f}{:12.7f}\n".format(labels[n], coords[n,0], coords[n,1], coords[n,2]))  
    fh.close()


  @staticmethod
  def output_coords_gb(labels, a0, cx, cy, cz, coords, path, filename):
    fh = open(os.path.join(path, filename), 'w')
    fh.write("{:6d}\n".format(len(labels)))
    fh.write("{:7.3f} {:3.1f} {:3.1f} {:3.1f}\n".format(a0, cx, cy, cz))
    for n in range(len(labels)):
      fh.write("{:3s}{:12.7f}{:12.7f}{:12.7f}\n".format(labels[n], coords[n,0] * cx, coords[n,1] * cy, coords[n,2] * cz))  
    fh.close()


  @staticmethod
  def padl(inp, plen=17):
    while(len(inp) < plen):
       inp = " " + inp
    return inp
    
  @staticmethod
  def padr(inp, plen=17):
    while(len(inp) < plen):
       inp = inp + " "
    return inp
    
makecoords.make()



