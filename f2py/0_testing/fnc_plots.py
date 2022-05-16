
#export PYTHONPATH=$PYTHONPATH:"/cloud/Code/python/eampa/f2py_lib"


from f_toolbox import fnc
#from eampa_lib.f_spline import spline
import numpy
import matplotlib.pyplot as plt




class make_plots:


  
  @staticmethod
  def run():
  
    xa = numpy.linspace(0.0, 7.0, 100) 
    xb = numpy.linspace(0.5, 7.0, 100) 
    xc = numpy.linspace(0.0, 2.0, 100) 
    xd = numpy.linspace(0.0, 120.0, 100)   
    xlabel = 'radius'
    ylabel = 'energy'
    
    
    # ZERO
    p = [0.0]
    pf = [0.0]           # No fixed parameters    
    make_plots.plot('zero', xa, p, pf, 'Zero', xlabel, ylabel) 
    make_plots.make_pot_file('zero', p, pf, 'zero')
    
    # BUCKINGHAM
    p = [6.0,0.5,12.0]
    pf = [0.0]           # No fixed parameters    
    make_plots.plot('buckingham', xb, p, pf, 'Buckingham', xlabel, ylabel) 
    make_plots.make_pot_file('buckingham', p, pf, 'buckingham')
    
    # LJ
    p = [2.3,3.5]
    pf = [0.0]           # No fixed parameters    
    make_plots.plot('lennard_jones', xb, p, pf, 'Lennard-Jones', xlabel, ylabel)
    make_plots.make_pot_file('lennard_jones', p, pf, 'lennard_jones') 
    
    # MORSE
    p = [4.669,1.256,2.8]
    pf = [0.0]            # No fixed parameters    
    make_plots.plot('morse', xb, p, pf, 'Morse', xlabel, ylabel) 
    make_plots.make_pot_file('morse', p, pf, 'morse') 
    
    # SLATER 4S
    p = [5.0, 1.323]
    pf = [0.0]
    make_plots.plot('slater_4s', xa, p, pf, 'Slater 4S', xlabel, ylabel) 
    make_plots.make_pot_file('slater_4s', p, pf, 'slater_4s') 
    
    # SLATER 4S
    p = [5.0, 1.323]
    pf = [4.0]
    make_plots.plot('slater_4s', xa, p, pf, 'Slater 4S Cutoff 4.0', xlabel, ylabel, 'slater_4s_cutoff') 
    make_plots.make_pot_file('slater_4s', p, pf, 'slater_4s_cutoff') 
    
    # SLATER 4S DOUBLE
    p = [5.0, 1.323, 3.0, 30.5]
    pf = [0.0]
    make_plots.plot('slater_4s', xa, p, pf, 'Slater 4S Double', xlabel, ylabel, 'slater_4s_double') 
    make_plots.make_pot_file('slater_4s', p, pf, 'slater_4s_double') 
    
    # ZBL
    p = [4.669,1.256,2.8]
    pf = [26.0, 26.0]  
    make_plots.plot('zbl', xa, p, pf, 'ZBL', xlabel, ylabel) 
    make_plots.make_pot_file('zbl', p, pf, 'zbl') 
    
    # ZBL
    p = [4.669,1.256,2.8]
    pf = [26.0, 26.0]  
    make_plots.plot('zbl', xa, p, pf, 'ZBL', xlabel, ylabel, 'zbl_zoom', 1.0, -1.0) 
    
    
    # SPLINE
    p = [-165.0, -78.49908, -78.15495, 1.8679553]
    pf =      [3, 1, 0, 0, 0] 
    pf = pf + [0, 0, 0, 0, 0] 
    pf = pf + [2.4325824, 2.86626, 3.0307584, 4.11246]      # First = power (cubic)       
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_1', 5.0, -10.0) 
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_1') 
    
    
    # SPLINE
    p = [3, -2, 1, -1.5]
    pf =      [3, 1, 26.0, 26.0, 0] 
    pf = pf + [0, 0, 0, 0, 0] 
    pf = pf + [2.4325824, 2.86626, 3.0307584, 4.11246]      # First = power (cubic)       
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_2', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_2') 
    
    
    # SPLINE
    p = [0, 0, 0, 0]
    pf =      [3, 1, 26.0, 26.0, 0] 
    pf = pf + [0, 0, 0, 0, 0] 
    pf = pf + [2.4325824, 2.86626, 3.0307584, 4.11246]      # First = power (cubic)       
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_3', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_3') 
    
    
    # SPLINE
    p = [-165.0, -78.49908, -78.15495, 1.8679553]
    pf =      [5, 1, 0, 0, 0] 
    pf = pf + [0, 0, 0, 0, 0] 
    pf = pf + [2.4325824, 2.86626, 3.0307584, 4.11246]      # First = power (quintic)
    make_plots.plot('spline', xa, p, pf, 'Quintic Spline', xlabel, ylabel, 'spline_quintic', 5.0, -10.0) 
    make_plots.make_pot_file('spline', p, pf, 'spline_quintic') 
    
    
    # SPLINE
    p =       [0, 0, 0, 0]
    pf =      [3, 1, 26.0, 26.0, 0] 
    pf = pf + [1.0, 0, 0, 0, 0] 
    pf = pf + [2.4325824, 2.86626, 3.0307584, 4.11246]      # First = power (cubic)       
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_4', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_4') 
    
    
    # SPLINE
    p = [0, 0, 0, 0]
    pf =      [3, 1, 26.0, 26.0, 6] 
    pf = pf + [1.0, 0, 0, 0, 0] 
    pf = pf + [2.4325824, 2.86626, 3.0307584, 4.11246]      # First = power (cubic)       
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_5', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_5') 
    
    
    
    # SPLINE
    p =       [0, 0, 0, 0]
    p = p +   [0, 0, 0, 0]
    pf =      [3, 1, 26.0, 26.0, 6] 
    pf = pf + [1.0, 0, 0, 0, 0] 
    pf = pf + [2.5, 2.5, 2.5, 2.5]      # First = power (cubic)   
    pf = pf + [3.5, 3.5, 3.5, 3.5]      # First = power (cubic)      
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_6', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_6') 
    
    
    # SPLINE
    p =       [0, 0, 0, 0]
    p = p +   [0, 0, 0, 0]
    p = p +   [0, 0, 0, 0]
    pf =      [3, 1, 26.0, 26.0, 6] 
    pf = pf + [1.0, 0, 0, 0, 0] 
    pf = pf + [2.8, 2.8, 2.8, 2.8]  
    pf = pf + [4.1, 4.1, 4.1, 4.1]    
    pf = pf + [4.9, 4.9, 4.9, 4.9]  
    make_plots.plot('spline', xa, p, pf, 'Cubic Spline', xlabel, ylabel, 'spline_cubic_7', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'spline_cubic_7') 
    
    
    # KNOT SPLINE
    p = [-0.5, 0.03, -0.1, 0.0, 0.2]
    p = p + [0.1, -0.1, 0.0, 0.1, 0.0 ]
    p = p + [0.1, -0.1, 0.0, 0.1, 0.0 ]
    settings = [5, 0, 0]
    settings = settings + [1, 0.0, 0.0, 0.0, 0.0] 
    settings = settings + [1, 6.5, 0, 0, 0] 
    settings = settings + [1, 0.5, 1.5, 26.0, 26.0]
    settings = settings + [0, 0]
    nodes = [2.0, 3.0, 4.0, 4.5, 5.0]    
    pf = settings + nodes       
    make_plots.plot('knot_spline', xa, p, pf, 'Knot Spline', xlabel, ylabel, 'knot_spline_1', 3.0, -3.0) 
    make_plots.make_pot_file('spline', p, pf, 'knot_spline_1') 
        
    
    # KNOT SPLINE
    p = [-0.5, 0.03, -0.1, 0.0, 0.2]
    p = p + [0.1, -0.1, 0.0, 0.1, 0.0 ]
    p = p + [0.1, -0.1, 0.0, 0.1, 0.0 ]
    settings = [5, 0, 0]
    settings = settings + [1, 0.0, 0.0, 0.0, 0.0] 
    settings = settings + [1, 6.5, 0, 0, 0] 
    settings = settings + [1, 0.5, 1.5, 26.0, 26.0]
    settings = settings + [0, 0]
    nodes = [2.0, 3.0, 4.0, 4.5, 5.0]    
    pf = settings + nodes       
    make_plots.plot('knot_spline', xa, p, pf, 'Knot Spline', xlabel, ylabel, 'knot_spline_2', 3.0, -3.0) 
    make_plots.make_pot_file('spline', p, pf, 'knot_spline_2') 
    
            
    
    # KNOT SPLINE
    p =     [1.0, -0.1, 0.03, -0.1, 0.0, 0.2]
    p = p + [-10, 0.0, -0.1, 0.0, 0.1, 0.0 ]
    p = p + [1.0, 0.1, -0.1, 0.0, 0.1, 0.0 ]
    settings = [5, 0, 0]
    settings = settings + [1, 0.0, 0.0, 0.0, 0.0] 
    settings = settings + [1, 6.5, 0, 0, 0] 
    settings = settings + [1, 0.4, 1.2, 26.0, 26.0]
    settings = settings + [0, 0]
    nodes = [1.2, 2.0, 3.0, 4.0, 4.5, 5.0]    
    pf = settings + nodes       
    make_plots.plot('knot_spline', xa, p, pf, 'Knot Spline', xlabel, ylabel, 'knot_spline_3', 3.0, -3.0) 
    make_plots.make_pot_file('spline', p, pf, 'knot_spline_3') 
    
    
    # MISHIN DENSITY
    p =  [0.5, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0]
    pf = [6.0]      
    make_plots.plot('mishin_density', xa, p, pf, 'Mishin Density', xlabel, ylabel, 'mishin_density', 3.0, -3.0) 
    make_plots.make_pot_file('mishin_density', p, pf, 'mishin_density') 
    
    
    
    xlabel = 'density'
    
    # FS EMBEDDING
    p = [1.5]
    pf = [0.0]         # No fixed parameters   
    make_plots.plot('fs_embedding', xc, p, pf, 'FS Embedding', xlabel, ylabel, None, 3.0, -3.0) 
    make_plots.make_pot_file('fs_embedding', p, pf, 'fs_embedding') 
    
          
    # MENDELEV EMBEDDING
    p = [0.6]
    pf = [0.0]         # No fixed parameters   
    make_plots.plot('mendelev_embedding', xc, p, pf, 'Mendelev Embedding', xlabel, ylabel, None, 3.0, -3.0) 
    make_plots.make_pot_file('mendelev_embedding', p, pf, 'mendelev_embedding') 
    
          
    # ACKLAND EMBEDDING
    p = [-0.5, 1.0]
    pf = [0.0]         # No fixed parameters   
    make_plots.plot('ackland_embedding', xc, p, pf, 'Ackland Embedding', xlabel, ylabel, None, 3.0, -3.0) 
    make_plots.make_pot_file('ackland_embedding', p, pf, 'ackland_embedding') 
    
    
    
    
    xlabel = 'energy'
    # Alpha Iron Ackland Mendelev
    p =       [-27.444805994228, 15.738054058489, 2.2077118733936, -2.4989799053251, 4.2099676494795]
    p = p +   [0.77361294129713, 0.80656415937789, -2.3194358924605, 2.6577406128280, -1.0260416933564]
    p = p +   [0.35018615891957, -0.058531821042271, -0.0030458824556234]
    pf =      [3, -1.0, 26.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
    pf = pf + [2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.3,3.7,4.2,4.7,5.3]   
    p_abs = numpy.absolute(p)
    pl = p[:] - 0.1 * p_abs[:]
    pu = p[:] + 0.1 * p_abs[:]
    make_plots.plot('spline', xa, p, pf, 'Fe Spline - Pair', xlabel, ylabel, 'Fe_spline_pair', 15.0, -10.0)
    make_plots.plot('spline', xa, p, pf, 'Fe Spline - Pair', xlabel, ylabel, 'Fe_spline_pair_zoom', 2.0, -2.0)
    make_plots.make_pot_file('spline', p, pf, 'Fe_spline_pair', pl, pu)
    
    # Alpha Iron Ackland Mendelev
    p =       [-27.444805994228, 15.738054058489, 2.2077118733936, -2.4989799053251, 4.2099676494795]
    p = p +   [0.77361294129713, 0.80656415937789, -2.3194358924605, 2.6577406128280, -1.0260416933564]
    p = p +   [0.35018615891957, -0.058531821042271, -0.0030458824556234]
    pf =      [1.0, 2.05, 9734.2365892908, 26.0, 26.0, 7.4122709384068, -0.64180690713367,-2.6043547961722, 0.6262539393123, 0] 
    pf = pf + [2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.3,3.7,4.2,4.7,5.3] 
    make_plots.plot('ackland_spline', xa, p, pf, 'Fe Alpha Ackland - Pair', xlabel, ylabel, 'Fe_alpha_ackland_pair', 150.0, -10.0)
    make_plots.plot('ackland_spline', xa, p, pf, 'Fe Alpha Ackland - Pair', xlabel, ylabel, 'Fe_alpha_ackland_pair_zoom', 2.0, -2.0)
    make_plots.make_pot_file('ackland_spline', p, pf, 'Fe_alpha_ackland_pair')
    make_plots.plot('ackland_spline', xa, p, pf, 'Ackland Spline', xlabel, ylabel, 'ackland_spline', 15.0, -2.0)
    make_plots.make_pot_file('ackland_spline', p, pf, 'ackland_spline')

    # Alpha Iron Ackland Mendelev
    p = [-6.7314115586063e-4, 7.6514905604792e-8]
    pf = [0.0]         # No fixed parameters   
    make_plots.plot('ackland_embedding', xa, p, pf, 'Fe Alpha Ackland - Embedding', xlabel, ylabel, 'Fe_alpha_ackland_embed',1.0,-4.0) 
    make_plots.make_pot_file('ackland_embedding', p, pf, 'Fe_alpha_ackland_embed')
  

    xlabel = 'density'
    # Alpha Iron Ackland Mendelev
    p =       [11.686859407970, -0.01471074009883, 0.47193527075943]
    pf =      [3, -1.0, 0.0, 0, 0, 0, 0, 0, 0, 0] 
    pf = pf + [2.4, 3.2, 4.2]      # First = power (cubic)       
    make_plots.plot('spline', xa, p, pf, 'Fe Alpha Ackland - Density', xlabel, ylabel, 'Fe_alpha_ackland_dens', 5.0, -10.0)
    make_plots.make_pot_file('spline', p, pf, 'Fe_alpha_ackland_dens')


    
    # SPLINE
    p =       [84.370697642971, -370.62682398816, 655.67905839695, -549.25562909352, 177.09076987937]
    p = p +   [-3.0363100686117, -3.0909218331545, -4.4211545804424, -1.9205251191921, -0.31060214601543]
    p = p +   [0.93316319033507, -3.0761272816964, 3.5391800334269, -1.7626076559210, 0.32587864688719]
    pf =      [4, -1, 0, 0, 8] 
    pf = pf + [0, 0, 0, 0, 0] 
    pf = pf + [2.9, 2.9, 2.9, 2.9, 2.9]
    pf = pf + [5.0, 5.0, 5.0, 5.0, 5.0]
    pf = pf + [6.4, 6.4, 6.4, 6.4, 6.4]
    make_plots.plot('spline', xa, p, pf, 'K Spline', xlabel, ylabel, 'spline_k', 5.0, -10.0) 
    make_plots.make_pot_file('spline', p, pf, 'spline_k') 

    
    # SPLINE
    p =       [-0.88124070590302, 6.5448285611211e-7, -2.5046137745775e-6, 2.3771143627356e-5, -4.4929936875438e-5, 2.3927653989776e-5]
    pf =      [65, 75, 90, 95, 100] 
    make_plots.plot('spline_embedding', xd, p, pf, 'Spline Embedding', xlabel, ylabel, 'spline_embedding', 5.0, -10.0) 
    make_plots.make_pot_file('spline_embedding', p, pf, 'spline_embedding') 



    

  @staticmethod
  def plot(f_name, x, p, pf, plot_title, xlabel, ylabel, plot_name=None, ymax=None, ymin=None):
    p = numpy.asarray(p, dtype=numpy.float64)
    pf = numpy.asarray(pf, dtype=numpy.float64)
    y = fnc.fv(f_name, x, p, pf)
    dy = fnc.fgradv(f_name, x, p, pf)
    
    if(ymin == None):
      ymin = max(-100, min(1.05 * numpy.amin(y[:]), 1.05 * numpy.amin(dy[:])))
    if(ymax == None):
      ymax = min(100, max(1.05 * numpy.amax(y[:]), 1.05 * numpy.amax(dy[:])))
    if(ymin == ymax):
      ymin = ymin - 0.1
      ymax = ymax + 0.1
    
    if(plot_name == None):
      plot_name = f_name
    
    plt.clf()   
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig, axs = plt.subplots(1, 1, figsize=(7,5), dpi=144)
    fig.tight_layout(pad=5.0)    
    fig.suptitle(plot_title)
    axs.plot(x, y, color='k', ls='solid', label='potential')
    axs.plot(x, dy, color='k', ls='dashed', label='gradient')
    #axs.plot(x, y_a, color='k', ls='dashed', label='gradient')
    #axs.set_ylim(-5.0,5.0)
    axs.legend()
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_ylim([ymin, ymax])
    plt.savefig('test/' + plot_name + '.png')
    plt.close(fig)
    plt.clf() 

  @staticmethod
  def make_pot_file(f_name, p, pf, pot_name=None, pl=None, pu=None):
    if(pot_name == None):
      pot_name = f_name 
      
    fh = open('test/' + pot_name + ".pot", 'w')
    fh.write('#TYPE   ' + f_name + '\n')
    fh.write('#P      ')   
    for n in range(len(p)):
      fh.write(str(p[n]) + ' ')
      if(n > 0 and n < (len(p) - 1) and (n+1) % 5 == 0):
        fh.write(' &\n        ')
    fh.write('\n')
    fh.write('#PF     ')   
    for n in range(len(pf)):
      fh.write(str(pf[n]) + ' ')
      if(n > 0 and n < (len(pf) - 1) and (n+1) % 5 == 0):
        fh.write(' &\n        ')
    fh.write('\n')
    if(pl is not None):
      fh.write('#PL     ')   
      for n in range(len(pl)):
        fh.write(str(pl[n]) + ' ')
        if(n > 0 and n < (len(pl) - 1) and (n+1) % 5 == 0):
          fh.write(' &\n        ')
      fh.write('\n')
    if(pu is not None):
      fh.write('#PU     ')   
      for n in range(len(pu)):
        fh.write(str(pu[n]) + ' ')
        if(n > 0 and n < (len(pu) - 1) and (n+1) % 5 == 0):
          fh.write(' &\n        ')
      fh.write('\n')        
    fh.close()


make_plots.run()