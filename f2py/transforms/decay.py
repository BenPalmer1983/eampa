import os
import numpy
from pz import pz
from isotopes import isotopes
import matplotlib.pyplot as plt
import copy
import hashlib

class decay:

  path_isotopes = "../data/isotopes.pz"
  loaded = False

  @staticmethod
  def set(path_isotopes):
    decay.path_isotopes = path_isotopes
    isotopes.set(path_isotopes)
    decay.load()

  @staticmethod
  def load():
    if(decay.loaded == False):
      decay.loaded = True


  @staticmethod
  def chain_isotopes(key, out=[]):
    if(not isotopes.is_valid(key)):
      return out
    if(key not in out):
      out.append(key)
    else:
      return out
    if(isotopes.is_stable(key)):
      return out
    else:
      dm = isotopes.get_decay_modes(key)
      for k in dm.keys():        
        decay.chain_isotopes(k, out)
      return out


  @staticmethod
  def make_chain(key, l=0, out=[], bf=0.0):
    if(not isotopes.is_valid(key)):
      return out
    if(l == 0):   
      decay.chains_store = []
    out.append([l, key, bf])
    if(isotopes.is_stable(key)):
      if(len(out) > 1):
        for i in range(len(out)-1,1,-1):
          if(out[i-1][0]>=out[i][0]):
            out.pop(out[i-1][0])
      return decay.chains_store.append(copy.deepcopy(out))
    else:
      l = l + 1 
      dms = isotopes.get_decay_modes(key)
      for k in dms.keys():
        bf = dms[k]['branching_factor']
        decay.make_chain(k, l, out, bf)
      return out      

  @staticmethod
  def calculate(parent, time, i_data_in, log=None, custom_chain=None):

    if(log != None):
      log_dir = decay.get_file_dir(log)
      print(log_dir)
      decay.make_dir(log_dir)

    decay.results = {
                    'tally': {},
                    'unique': None,
                    'chains': None, 
                    } 
    decay.chains_store = None

    if(custom_chain == None):
      decay.chains_store = []
      decay.make_chain(parent, 0, [])
      decay.results['chains'] = []
      cn = 0
      for chain in decay.chains_store:
        decay.results['chains'].append([])
        for iso in chain:
          k = iso[1]
          bf = iso[2]
          i_data = isotopes.get(k)
          half_life = None
          decay_constant = None
          n0 = 0.0
          w = 0.0
          if(not i_data['stable']):
            half_life = i_data['half_life'] 
            decay_constant = i_data['decay_constant']           
          if(k in i_data_in.keys()):
            n0 = i_data_in[k]['n0']
            w = i_data_in[k]['w']
          d = {
              'isotope_key': k, 
              'bf': bf,
              'w': w,
              'n0': n0,
              'nend': 0,
              'half_life': half_life,
              'decay_constant': decay_constant,
              } 
          decay.results['chains'][cn].append(d)
        cn = cn + 1
    else:
      decay.results['chains'] = custom_chain


    # Find unique and make tally
    decay.results['unique'] = []
    for chain in decay.results['chains']:
      for iso in chain:
        k = iso['isotope_key']
        if(k not in decay.results['tally'].keys()):
          decay.results['unique'].append(k)
          # Use provided data
          if(iso['half_life'] == None):
            stable = True
            half_life = None
            decay_constant = 0.0
          else:
            stable = False
            half_life = iso['half_life']
            decay_constant = numpy.log(2) / iso['half_life']
          w = iso['w']
          n0 = iso['n0'] 
          # Get proton/neutron etc from isotopes database
          i_data = isotopes.get(k)
          if(i_data is None):     
            decay.results['tally'][k] = {
                                        'printable': 'Custom',
                                        'element': 'ZZ',
                                        'protons': 999,
                                        'nucleons': 999,
                                        'metastable': 9,
                                        'stable': stable,
                                        'half_life': half_life,
                                        'decay_constant': decay_constant,
                                        'w': w,
                                        'n0': n0,
                                        'nend': 0.0,
                                        }
          else:
            decay.results['tally'][k] = {
                                        'printable': decay.pad(isotopes.get_printable_name(k), 12),
                                        'element': i_data['element'],
                                        'protons': i_data['protons'],
                                        'nucleons': i_data['nucleons'],
                                        'metastable': i_data['metastable'],
                                        'stable': stable,
                                        'half_life': half_life,
                                        'decay_constant': decay_constant,
                                        'w': w,
                                        'n0': n0,
                                        'nend': 0.0,
                                        }


    decay.results['chains_individual'] = []
    for cn in range(len(decay.results['chains'])):
      for n in range(len(decay.results['chains'][cn])):
        nc = []
        if(decay.results['chains'][cn][n]['n0']>0.0 or decay.results['chains'][cn][n]['w']>0.0):
          for j in range(n, len(decay.results['chains'][cn])):
            iso = copy.deepcopy(decay.results['chains'][cn][j])
            if(j>n):
              iso['n0'] = 0.0
              iso['w'] = 0.0
            nc.append(iso)
        if(len(nc)>0):
          decay.results['chains_individual'].append(nc)

    for cn in range(len(decay.results['chains_individual'])):
      chain = decay.results['chains_individual'][cn]
      n0 = numpy.zeros((len(chain),),)
      w = numpy.zeros((len(chain),),)
      l = numpy.zeros((len(chain),),)
      b = numpy.zeros((len(chain)-1,),)
      for n in range(len(decay.results['chains_individual'][cn])):
        k = decay.results['chains_individual'][cn][n]['isotope_key']
        n0[n] = decay.results['chains_individual'][cn][n]['n0']
        w[n] = decay.results['chains_individual'][cn][n]['w']
        l[n] = decay.results['tally'][k]['decay_constant']
        if(n>0):
          b[n-1] = decay.results['chains_individual'][cn][n]['bf']
        # Some observationally stable will still have a decay constant, so set to -1.0
        if(isotopes.is_stable(k)):
          l[n] = -1.0
      nt = decay.calculate_activity(time, l, b, w, n0) 
      for n in range(len(decay.results['chains_individual'][cn])):
        decay.results['chains_individual'][cn][n]['nend'] = nt[n]

    set = []
    for cn in range(len(decay.results['chains_individual'])):
      ckey = ''
      for n in range(len(decay.results['chains_individual'][cn])):
        k = decay.results['chains_individual'][cn][n]['isotope_key']
        ckey = ckey + str(decay.results['chains_individual'][cn][n]['isotope_key']) + "N0:" + str(decay.results['chains_individual'][cn][n]['n0']) + "W:" + str(decay.results['chains_individual'][cn][n]['w'])
        ckeyh = hashlib.md5(ckey.encode())
        ckeyh = ckeyh.hexdigest()
        if(ckeyh not in set):
          decay.results['tally'][k]['nend'] = decay.results['tally'][k]['nend'] + decay.results['chains_individual'][cn][n]['nend']
          set.append(ckeyh)


    # Log
    if(log != None):
      width = 140
      fh = open(log, 'w')
      fh.write("Unique Isotopes\n")
      fh.write(decay.hr(width) + "\n")
      for k in decay.results['unique']: 
        line = decay.results['tally'][k]['printable']
        fh.write(line + "\n")
      fh.write("\n")
      fh.write("\n")
      fh.write("Decay Chains\n")
      fh.write(decay.hr(width) + "\n")
      fh.write("\n")
      fh.write("\n")
      for cn in range(len(decay.results['chains'])):
        chain = decay.results['chains'][cn]

        fh.write(decay.pad(cn+1,6)) 
        for n in range(len(chain)):
          iso = chain[n]
          k = iso['isotope_key']
          if(n>0):
            bf = "{0:3e}".format(iso['bf'])
            fh.write(" --[" + str(bf) + "]--> ")
          fh.write(decay.pad(decay.results['tally'][k]['printable'], 12))
        fh.write("\n")
        fh.write(decay.pad("T1/2",6)) 
        for n in range(len(chain)):
          if(n>0):
            fh.write(decay.pad("",17)) 
          if(decay.results['chains'][cn][n]['half_life'] == None):
            fh.write(decay.pad("[Stable]",12))
          else:
            fh.write(decay.pad("[" + str("{0:8e}".format(decay.results['chains'][cn][n]['half_life'])) + "]",12))
        fh.write("\n")
        fh.write(decay.pad("L",6)) 
        for n in range(len(chain)):
          if(n>0):
            fh.write(decay.pad("",17)) 
          if(decay.results['chains'][cn][n]['decay_constant'] == None):
            fh.write(decay.pad("[Stable]",12))
          else:
            fh.write(decay.pad("[" + str("{0:8e}".format(decay.results['chains'][cn][n]['decay_constant'])) + "]",12))
        fh.write("\n")
        fh.write("\n")        
        fh.write("\n")


     
      fh.write("Amounts\n")
      fh.write(decay.hr(width) + "\n")
      fh.write("\n")
      fh.write(decay.hr(width) + "\n")
      line = decay.pad("Isotope", 12)
      line = line + decay.pad("T(1/2)", 18)
      line = line + decay.pad("Decay Constant", 18)
      line = line + decay.pad("W", 18)
      line = line + decay.pad("N(t=0)", 18)
      line = line + decay.pad("N(t=" + str(time) + ")", 18)
      line = line + decay.pad("A(t=0)", 18)
      line = line + decay.pad("A(t=" + str(time) + ")", 18)
      fh.write(line + "\n")
      fh.write(decay.hr(width) + "\n")


      for k in decay.results['tally'].keys():
        line = decay.pad(decay.results['tally'][k]['printable'], 12)
        if(decay.results['tally'][k]['half_life'] is None):
          line = line + decay.pad("Stable", 18)
          line = line + decay.pad("Stable", 18)
        else:
          line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['half_life'])).strip(), 18)
          line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['decay_constant'])).strip(), 18)
        line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['w'])).strip(), 18)
        line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['n0'])).strip(), 18)
        line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['nend'])).strip(), 18)
        if(decay.results['tally'][k]['half_life'] is not None):
          line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['decay_constant'] * decay.results['tally'][k]['n0'])).strip(), 18)
          line = line + decay.pad(str("{0:16e}".format(decay.results['tally'][k]['decay_constant'] * decay.results['tally'][k]['nend'])).strip(), 18)
        fh.write(line + "\n")
      fh.write(decay.hr(width) + "\n")
      fh.write("\n")
      fh.write("\n")
      
      fh.close()
  
    return decay.results

  ##########################################
  # DECAY EQUATIONS
  ##########################################
    
  @staticmethod
  def calculate_activity(t, lam, b, w, n0):
    nt = numpy.zeros((len(n0),),)
    for m in range(0,len(n0)):
      if(lam[m] <= 0):
        nt[m] = decay.calc_stable(t, m, lam, b, w, n0)
      else:
        nt[m] = decay.calc_unstable(t, m, lam, b, w, n0)
    return nt

  @staticmethod
  def calc_unstable(t, m, lam, b, w, n0):
    y = 0.0 
    k = 0
    while(k<=m):    
      y = y + decay.r(k, m, lam, b) * (decay.f_unstable(t, k, m, lam) * n0[k] + decay.g_unstable(t, k, m, lam) * w[k])
      k = k + 1
    return y

  @staticmethod
  def f_unstable(t, k, m, lam): 
    s = 0.0
    i = k
    while(i<=m):
      p = decay.lprod(lam, k, m, 3, i)
      s = s + numpy.exp(-lam[i] * t) * (1.0 / p)
      i = i + 1
    s = s * (-1)**(m-k)  
    return s


  @staticmethod
  def g_unstable(t, k, m, lam): 
    # Term a
    p = decay.lprod(lam, k, m, 1)
    a = 1.0 / p
    # Term b
    s = 0.0
    i = k
    while(i<=m):
      p = decay.lprod(lam, k, m, 3, i)  
      s = s + (1 / lam[i]) * numpy.exp(-lam[i] * t) * (1.0 / p)
      i = i + 1
    s = s * (-1)**(m-k+1)           
    return a + s    

  @staticmethod
  def calc_stable(t, m, lam, b, w, n0):
    y = n0[m] + w[m] * t 
    k = 0
    while(k<=m-1):    
      y = y + decay.r(k, m, lam, b) * (decay.f_stable(t, k, m, lam) * n0[k] + decay.g_stable(t, k, m, lam) * w[k])
      k = k + 1
    return y    

  @staticmethod
  def f_stable(t, k, m, lam):  
    mm = m - 1
    p = decay.lprod(lam, k, mm, 1)
    a = 1.0 / p

    s = 0.0
    i = k
    while(i<=mm):
      p = decay.lprod(lam, k, mm, 3, i)
      s = s + (1 / lam[i]) * numpy.exp(-lam[i] * t) * (1.0 / p)
      i = i + 1
    s = s * (-1)**(mm-k+1)  
    return a + s
   
  @staticmethod 
  def g_stable(t, k, m, lam): 
    mm = m - 1

    p = decay.lprod(lam, k, mm, 1) 
    a = (t/p)

    s = 0.0
    i = k
    while(i<=mm):
      p = decay.lprod(lam, k, mm, 4, i) 
      s = s + p
      i = i + 1
    q = decay.lprod(lam, k, mm, 2)
    b = (-s / q)

    c = 0.0
    i = k
    while(i<= mm):
      p = lam[i]**2 * decay.lprod(lam, k, mm, 3, i) 
      c = c + numpy.exp(-1.0 * lam[i] * t) / p
      i = i + 1
    c = c * (-1.0)**(mm - k)
    nd = a + b + c

    return nd       
    

  @staticmethod
  def r(k, m, lam, b):
    if(k == m):
      return 1  
    else:
      p = 1.0
      i = k
      while(i<=(m-1)):
        p = p * b[i] * lam[i]
        i = i + 1
      return p
 
  @staticmethod   
  def lprod(lam, k, m, t=1, i=None): 
    # PROD k,m lam[j]
    if(t == 1):
      p = 1.0
      j = k
      while(j<=m):
        p = p * lam[j]
        j = j + 1
      return p
    # PROD k,m lam[j]**2
    elif(t == 2):
      p = 1.0
      j = k
      while(j<=m):
        p = p * lam[j]**2
        j = j + 1
      return p
    # PROD k,m, i!=j (lam[i]-lam[j])
    elif(t == 3):
      p = 1.0
      j = k
      while(j<=m):
        if(i != j):
          p = p * (lam[i] - lam[j])
        j = j + 1
      return p
    # PROD k,m, i!=j lam[j]
    elif(t == 4):
      p = 1.0
      j = k
      while(j<=m):
        if(i != j):
          p = p * lam[j]
        j = j + 1
      return p
  
  @staticmethod
  def pad(inp, l=16):
    inp = str(inp)
    while(len(inp)<l):
      inp = inp + " "
    return inp

  @staticmethod
  def hr(l=16):
    inp = ""
    while(len(inp)<l):
      inp = inp + "="
    return inp

  @staticmethod
  def get_file_dir(file_path):
    file_path = file_path.strip()
    if(file_path[0] != "/"):
      root = os.getcwd()
      file_path = root + "/" + file_path
    file_path = file_path.split("/")
    path = ""
    for i in range(1,len(file_path) - 1):
      path = path + "/" + file_path[i]
    return path

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
  def test():
    print("Test Decay")

    # Load isotopes dictionary
    decay.set("../data/isotopes.pz")

    idata = {}
    parent = 84216
    time = 10
    idata[84216] = {'w': 0.20, 'n0': 100.0}
    idata[82212] = {'w': 0.0, 'n0': 5.0}
    idata[83212] = {'w': 0.07, 'n0': 15.0}
    idata[81208] = {'w': 0.005, 'n0': 0.0}
    # 84Po212 0 0 (default)
    idata[82208] = {'w': 0.01, 'n0': 300.0}
    decay.calculate(parent, time, idata, "testing/log_84216_new.txt")


def main():
  decay.test()

if __name__ == "__main__":
    main()    

