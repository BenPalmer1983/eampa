from ds import ds

class g:


  # directories
  dir = {'out': None, }


  # file handles
  fh = None        # g.fh is the main log file

  # Input file
  input = {}

  # Potential
  potential = {}

  # Shared lists
  labels = []
  groups = []
  nl = []
  configs = []
  one_line_configs = []
  bp = []
  bp_settings = ds.bp_settings()

  rss_w = ds.rss_w()


  # timers
  start = None
  end = None

  stats = ds.stats()

  # Verbosity
  # usual setting 0, 1
  verbose = {'print': 0, 'file': 1}
  # 0 always print
  # 1 best write to file
  # 2 detailed, write to file but don't usually log
  # 3 log for testing

  # Display (fitting)
  display = 1
  
   
  """#####################################
  Speed
  #####################################"""
  
  calc_counter = 0
  calc_time = 0.0


  """#####################################
  Potfit
  #####################################"""

  potfit = ds.potfit()
  potfit_step_types = ds.potfit_step_types()
  potfit_steps = []
  maxtime = 1.0e6
  displaynote = ''
  potfit_bp = []     # Compare fitting to exact eos and ec
  potfit_bestsaveperiod = 60
