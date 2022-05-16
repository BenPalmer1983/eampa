import os
import numpy

class processpw:

  input_dir = []
  output_dir = []
  dft = {}
  configs = []
  mask = {}
  rydberg = 13.6056980659
  relaxed = 'all'
  atoms_min = 32
 
  @staticmethod
  def run():

    processpw.read_input()
    processpw.read_files()
    processpw.output() 

  @staticmethod
  def read_input():

    print("Read input file")
    section = None
    fh = open('input.in', 'r')
    inp = []
    for line in fh:
      line = line.strip()

      # Which section of the input file
      if(processpw.starts_with(line, "#")):
        section = None
      if(line == "#INPUT"):
        section = "input"
      elif(line == "#OUTPUT"):
        section = "output"
      elif(line == "#DFT"):
        section = "dft"
      elif(line == "#MASK"):
        section = "mask"
      elif(line == "#RELAXED"):
        section = "relaxed"

      # Input
      if(section == "input" and not (processpw.starts_with(line, "#") or line == "")):
        processpw.input_dir.append(line)

      # DFT
      if(section == "dft"):
        if("=" in line):
          f = line.split("=")
          label = f[0].upper()
          f = f[1].split(",")
          processpw.dft[label] = f

      # MASK
      if(section == "mask"):
        if("=" in line):
          f = line.split("=")
          label_mask = f[0].upper()
          f = f[1].split(",")
          for fm in f:
            label = fm.upper().strip()
            if(label not in processpw.mask):
              processpw.mask[label] = label_mask

      # RELAXED
      if(section == "relaxed" and not (processpw.starts_with(line, "#") or line == "")):
        processpw.relaxed = line
    fh.close()  

  @staticmethod
  def read_files():
    for n in range(len(processpw.input_dir)):
      input_dir = processpw.input_dir[n]
      output_dir = processpw.input_dir[n] + "_processed"
      print("Dir: ", input_dir)
      if(not os.path.isdir(input_dir)):
        print("Input directory does not exist")
      else:
        files = processpw.list_files(input_dir, ext='.out')
        for file_name_in in files:
          file_name_out = output_dir + file_name_in[len(input_dir):]
          os.makedirs(os.path.dirname(file_name_out), exist_ok=True)
          processpw.read_file(file_name_in, os.path.dirname(file_name_out))



  @staticmethod
  def read_file(file_name_in, file_dir_out):
    # Read File
    f = []
    fh = open(file_name_in, 'r')
    for line in fh:
      if(line.strip() != ""):
        f.append(line[:-1])
    fh.close()

    # Check file
    job_done = False
    converged = True
    n = 0
    while(n<len(f)):
      if("JOB DONE" in f[n]):
        job_done = True
      if("convergence NOT achieved" in f[n]):
        converged = False
      n = n + 1

    # Output to User
    print(file_name_in, end=": \t")
    if(job_done):
      print("Job Done", end=" \t")
    else:
      print("Job NOT Done", end=" \t")
    if(converged):
      print("Convergence Achieved", end=" \t")
    else:
      print("#### Convergence FAILED - Skipping Failed ####")
      return None
    print()

    # Load number of atoms
    nat = 0
    n = 0
    while(n<len(f)):
      l = f[n]
      if("number of atoms/cell" in l):
        nat = int(l[35:].strip())
        n = len(f)
      n = n + 1

    a0 = []
    uv = []
    coords = []
    total_energy = []
    forces = []
    stresses = []

    n = 0
    while(n<len(f)):
      l = f[n]
      # Load all a0 & uv
      if("lattice parameter" in l):
        a0.append(float(l[34:47].strip()) * 0.529) # CONVERT TO ANGSTROM
      elif("CELL_PARAMETERS" in l):
        a0.append(float(l[22:34].strip()) * 0.529) # CONVERT TO ANGSTROM
        uv.append(numpy.zeros((3,3,),))
        for i in range(3):
          n = n + 1
          l = f[n]
          uv[-1][i,0] = float(l[2:15].strip())
          uv[-1][i,1] = float(l[16:29].strip())
          uv[-1][i,2] = float(l[30:].strip())
      elif("crystal axes: (cart. coord. in units of alat)" in l):
        n = n + 1
        uv.append(numpy.zeros((3,3,),))
        for i in range(3):
          l = f[n]
          uv[-1][i,0] = float(l[25:36].strip())
          uv[-1][i,1] = float(l[37:47].strip())
          uv[-1][i,2] = float(l[47:57].strip()) 
          n = n + 1
      elif("   Cartesian axes" in l):
        coords.append([])
        n = n + 2
        l = f[n]
        while(not "number of k points" in f[n]):
          l = f[n]
          label = l[10:25].strip().upper()
          x = float(l[40:51].strip())
          y = float(l[52:63].strip())
          z = float(l[64:75].strip())
          if(label in processpw.mask):
            label = processpw.mask[label]
          coords[-1].append({})
          coords[-1][-1]['label'] = label
          coords[-1][-1]['x'] = x
          coords[-1][-1]['y'] = y
          coords[-1][-1]['z'] = z
          n = n + 1
      elif("ATOMIC_POSITIONS (crystal)" in l):
        coords.append([])
        n = n + 1
        for i in range(nat):
          l = processpw.one_space(f[n])
          d = l.split(" ")
          label = d[0].upper()
          x = float(d[1])
          y = float(d[2])
          z = float(d[3])
          if(label in processpw.mask):
            label = processpw.mask[label]
          coords[-1].append({})
          coords[-1][-1]['label'] = label
          coords[-1][-1]['x'] = x
          coords[-1][-1]['y'] = y
          coords[-1][-1]['z'] = z
          n = n + 1

      elif("!    total energy   " in l):
        total_energy.append(float(l[32:50].strip())  * processpw.rydberg)  # In EV

      elif("Forces acting on atoms (cartesian axes, Ry/au):" in l):
        forces.append([])
        n = n + 1
        while("atom" in f[n] and "force" in f[n]):
          fx = float(f[n][34:49].strip()) * 25.710    # Ry/au -> eV/angstrom
          fy = float(f[n][50:63].strip()) * 25.710    # Ry/au -> eV/angstrom
          fz = float(f[n][64:].strip()) * 25.710    # Ry/au -> eV/angstrom
          forces[-1].append({})
          forces[-1][-1]['fx'] = fx
          forces[-1][-1]['fy'] = fy
          forces[-1][-1]['fz'] = fz          
          n = n + 1

      elif("          total   stress  (Ry/bohr**3) " in l):
        n = n + 1
        stresses.append(numpy.zeros((3,3,),)) 
        for i in range(3):
          l = processpw.one_space(f[n]).strip()
          l = l.split(" ")
          stresses[-1][i,0] = float(l[0]) *  14583.6    # Ry/bohr^3 -> gpa
          stresses[-1][i,1] = float(l[1]) * 14583.6    # Ry/bohr^3 -> gpa
          stresses[-1][i,2] = float(l[2]) * 14583.6    # Ry/bohr^3 -> gpa 

          """
          stresses[-1][i,0] = float(l[3:15].strip()) *  14583.6    # Ry/bohr^3 -> gpa
          stresses[-1][i,1] = float(l[15:28].strip()) * 14583.6    # Ry/bohr^3 -> gpa
          stresses[-1][i,2] = float(l[28:42].strip()) * 14583.6    # Ry/bohr^3 -> gpa     
          """
          n = n + 1 
      n = n + 1

    if(len(coords) == 0):
      print("#### Error: coords unreadable ####")
    elif(len(a0) == 0):
      print("#### Error: a0 unreadable ####")

    elif(len(a0) == 1):      
      processpw.get_config(file_name_in, file_dir_out, None, nat, a0[0], uv[0], coords[0], total_energy[0], forces[0], stresses[0], True)
    elif(len(a0) > 1):
      if(processpw.relaxed == 'all'):
        # first
        i = 0
        processpw.get_config(file_name_in, file_dir_out, i, nat, a0[i+1], uv[i+1], coords[i+1], total_energy[i], forces[i], stresses[i], False)

        # middle
        for i in range(1, len(a0) - 2):
          processpw.get_config(file_name_in, file_dir_out, i, nat, a0[i], uv[i], coords[i+1], total_energy[i], forces[i], stresses[i], False)

        # last
        i = len(a0) - 2
        processpw.get_config(file_name_in, file_dir_out, i, nat, a0[i+1], uv[i+1], coords[i+1], total_energy[i], forces[i], stresses[i], True)
        """
        i = 0
        istr = str(i)
        while(len(istr)<3):
          istr = "0" + istr
        fn = filename.replace(".out", "_" + istr + ".out")
        cg = processpw.get_config(file_name_in, fn, output_dir, nat, a0[i+1], uv[i+1], coords[i+1], total_energy[i], forces[i], stresses[i], False)
        processpw.configs.append(cg)
        for i in range(1, len(a0) - 2):
          istr = str(i)
          while(len(istr)<3):
            istr = "0" + istr
          fn = filename.replace(".out", "_" + istr + ".out")
          fn.replace(".out", "_" + str(i) + ".out")
          cg = processpw.get_config(file_name_in, fn, output_dir, nat, a0[i], uv[i], coords[i+1], total_energy[i], forces[i], stresses[i], False)
          processpw.configs.append(cg)
        i = len(a0) - 2
        istr = str(i)
        while(len(istr)<3):
          istr = "0" + istr
        fn = filename.replace(".out", "_" + istr + ".out")
        cg = processpw.get_config(file_name_in, fn, output_dir, nat, a0[i+1], uv[i+1], coords[i+1], total_energy[i], forces[i], stresses[i], True)
        processpw.configs.append(cg)
        """
      else:
        i = len(a0) - 2
        cg = processpw.get_config(file_name_in, fn, output_dir, nat, a0[i+1], uv[i+1], coords[i+1], total_energy[i], forces[i], stresses[i], True)


  @staticmethod
  def get_config(file_name_in, file_dir_out, counter, nat, a0, uv, coords, total_energy, forces, stresses, to_crystal=False):

    file_name_out = processpw.make_out_name(file_name_in, file_dir_out, counter)
    print(file_name_out)
    print(total_energy)
    if(total_energy is not None):
      coh = 0.0
      relaxed_dft_energy = 0.0
      for i in range(len(coords)):
        e = processpw.dft[coords[i]['label']]
        #print(e)
        if(e[2].upper() == 'EV'):
          relaxed_dft_energy = relaxed_dft_energy + float(e[1]) / float(e[0])
        elif(e[2].upper() == 'RY'):
          relaxed_dft_energy = relaxed_dft_energy + (float(e[1]) / float(e[0])) * processpw.rydberg
        if(e[4].upper() == 'EV'):
          coh = coh + float(e[3])
        elif(e[4].upper() == 'RY'):
          coh = coh + float(e[3]) * processpw.rydberg

      #print()

      energy = coh + (total_energy - relaxed_dft_energy)
      epa = energy / nat   
      #print(energy)

    uv_inverse = numpy.linalg.inv(uv)
    uv_norm = numpy.copy(uv)
    uv_norm = uv_norm / uv[0,0]
    a0_norm = a0 * uv[0,0]
    
    coords_crystal = []      
    for i in range(len(coords)):
      c = numpy.zeros((3,),)
      label = coords[i]['label']
      if(to_crystal):
        c[0] = coords[i]['x']
        c[1] = coords[i]['y']
        c[2] = coords[i]['z']
        cc = numpy.matmul(uv_inverse, c)
        coords_crystal.append({'label': label, 'x': cc[0], 'y': cc[1], 'z': cc[2],})
      else:
        coords_crystal.append(coords[i])

    cg = {
            'file_name_in': file_name_in,
            'file_name_out': file_name_out,
            'nat': nat,
            'a0': a0,
            'a0_norm': a0_norm,
            'uv': uv,
            'uv_inverse': uv_inverse,
            'uv_norm': uv_norm,
            'coords': coords,
            'coords_crystal': coords_crystal,
            'total_energy': total_energy,
            'energy': energy,
            'epa': epa,
            'forces': forces,
            'stresses': stresses,
            }
    processpw.configs.append(cg)

  @staticmethod
  def output():
    print("Output Files")
    for n in range(len(processpw.configs)):
      file_name_out = str(processpw.configs[n]['file_name_out'])

      fh = open(file_name_out, 'w')
      print(file_name_out, end="   ")
      a0 = processpw.configs[n]['a0_norm']
      uv = processpw.configs[n]['uv_norm']
      epa = processpw.configs[n]['epa']
      nat = processpw.configs[n]['nat']
      coords = processpw.configs[n]['coords_crystal']
      forces = processpw.configs[n]['forces']
      stresses = processpw.configs[n]['stresses']

      weighting = 1.0
      if("w((" in processpw.configs[n]['file_name_in'] and "))w" in processpw.configs[n]['file_name_in']):
        ww = processpw.configs[n]['file_name_in'].split("w((")
        weighting = float(ww[1].split("))w")[0])


      c = 1
      if(len(coords) < processpw.atoms_min):
        c = int(numpy.ceil((processpw.atoms_min / len(coords))**(1/3)))
        print(c, end="  ")
      print()

      fh.write("/* Converted from PWscf output file */ \n")
      fh.write("#L_UNITS ang \n")
      fh.write("#E_UNITS eV \n")
      fh.write("#S_UNITS GPA \n")
      fh.write("#F_UNITS EV/ANG \n")
      fh.write("#ALAT 	" + processpw.float_str(c * a0)+ "     // a0 \n")
      fh.write("#X 	" + processpw.float_str(uv[0,0]) + "  " + processpw.float_str(uv[0,1]) + "  " + processpw.float_str(uv[0,2]) + "     // X \n")
      fh.write("#Y 	" + processpw.float_str(uv[1,0]) + "  " + processpw.float_str(uv[1,1]) + "  " + processpw.float_str(uv[1,2]) + "     // Y \n")
      fh.write("#Z 	" + processpw.float_str(uv[2,0]) + "  " + processpw.float_str(uv[2,1]) + "  " + processpw.float_str(uv[2,2]) + "     // Z \n")
      
      fh.write("#SX 	" + processpw.float_str(stresses[0,0]) + "  " + processpw.float_str(stresses[0,1]) + "  " + processpw.float_str(stresses[0,2]) + "     // SX \n")
      fh.write("#SY 	" + processpw.float_str(stresses[1,0]) + "  " + processpw.float_str(stresses[1,1]) + "  " + processpw.float_str(stresses[1,2]) + "     // SY \n")
      fh.write("#SZ 	" + processpw.float_str(stresses[2,0]) + "  " + processpw.float_str(stresses[2,1]) + "  " + processpw.float_str(stresses[2,2]) + "     // SZ \n")
     

      fh.write("#EPA 	" + processpw.float_str(epa) + "     // Energy per Atom \n")
      fh.write("#C 	1 1 1 \n")
      fh.write("#W 	" + str(weighting) + " \n")
      fh.write("#RCUT 	6.5        // Rcut \n")

      for x in range(c):
        for y in range(c):
          for z in range(c):
            for i in range(nat):
              fh.write(coords[i]['label'])
              fh.write('   ')
              fh.write(processpw.float_str((x + coords[i]['x']) / c))
              fh.write('   ')
              fh.write(processpw.float_str((y + coords[i]['y']) / c))
              fh.write('   ')
              fh.write(processpw.float_str((z + coords[i]['z']) / c))
              fh.write('   ')
              if(type(forces) == list):
                fh.write(processpw.float_str(forces[i]['fx']))
                fh.write('   ')
                fh.write(processpw.float_str(forces[i]['fy']))
                fh.write('   ')
                fh.write(processpw.float_str(forces[i]['fz']))
                fh.write('\n')
        

      fh.close()


  @staticmethod
  def make_out_name(in_name, out_dir, n=None):
    base = os.path.basename(in_name)
    if(n == None):
      return os.path.join(out_dir, base) 
    else:
      nstr = str(n)
      while(len(nstr)<4):
        nstr = "0" + nstr
      if('.' not in base):
        return os.path.join(out_dir, base + "_SCF" + nstr)
      else:
        rbase = ''
        for c in base:
          rbase = c + rbase
        base = ''
        for c in rbase:
          if(c == '.'):
            base = '_SCF' + nstr + '.' + base
          else:
            base = c + base
        return os.path.join(out_dir, base)

  @staticmethod
  def one_space(l):
    out = ''
    last = None
    for c in l:
      if(not (c == " " and last == " ")):
        out = out + c
      last = c
    return out

  @staticmethod
  def make_dir(dir):
    dirs = dir.split("/")
    try:
      dir = ''
      for i in range(len(dirs)):
        dir = dir + dirs[i]
        if(not os.path.exists(dir) and dir.strip() != ''):
          os.mkdir(dir) 
        dir = dir + '/'
      return True
    except:
      return False



  @staticmethod
  def float_str(x):
    return str('{:10.5f}'.format(float(x)))


  @staticmethod
  def starts_with(line, char):
    if(len(line) == 0):
      return False
    if(len(line) < len(char)):
      return False
    if(line[0:len(char)] == char):
      return True
    return False


  @staticmethod
  def list_files(path, ext='', out=[]):
    l = os.listdir(path)
    le = len(ext)
    for f in l:
      fpath = path + "/" + f
      if(os.path.isdir(fpath)):
        processpw.list_files(fpath, ext, out)
      else:
        if(le > 0):
          if(len(f) > le):
            if(fpath[-le:].lower() == ext.lower()):
              out.append(fpath)
        else:
          out.append(fpath)
    return out

processpw.run()
