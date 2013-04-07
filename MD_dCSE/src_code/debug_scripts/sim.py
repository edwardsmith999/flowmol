#!/usr/bin/env python

import os
import shutil
import re
import math
import struct
import time
import datetime
try:
  from numpy import *
  np = True
except ImportError:
  np = False


# ---------------------------------------------------------------
# Function average
#
# - Returns average of list of numbers
################################################################
def average(array, pick=''):
  total = 0.0
  l = len(array)
  picked = 0
  for i in range(l):
    if pick:
      if i % pick == 0:
        picked += 1
        total += array[i]
        continue
    else:
      total += array[i]
      picked += 1

  average = total/picked
  return average



# ---------------------------------------------------------------
# Function sigma
#
# - Returns standard deviation of list of numbers.
# - Takes an optional argument of the mean to save on computation.
################################################################# 
def sigma(array, mu='', pick= ''):
  total = 0.0
  picked = 0
  l = len(array)

  if not mu:
    mu = average(array)
  for i in range(l):
    elem = array[i]
    if pick:
      if i % pick == 0:
        picked += 1
        total += (elem - mu)*(elem - mu)
        continue
    else:
      picked += 1
      total += (elem - mu)*(elem - mu)
  sigma = total/math.sqrt(picked)
  return sigma



# ----------------------------------------------------------------------------
# Class Simulation
#
# Encapsulates the process of setting up, running and post 
# processing a simulation in the MD_dCSE code.
# 
# - Allows the creation of simulation objects.
#     - (default) assume code is stored in './..' 
#       e.g.    >>  sim1 = Simulation()
#     - (optional) code is stored in 'codedir'
#       e.g.    >>  sim2 = Simulation('codedir')
#
# - Sets all inputs to default values as in default.in.
#
# - Allows modification of these input values by default keywords
#       e.g.    >>  sim1.input('DENSITY') = 0.7
#               >>  sim2.input('DENSITY') = 0.9
#
# - There are also extra methods which do the same thing 
#       e.g.    >>  sim1.density(0.84)
#               >>  sim2.tplot(1)
#
# - The simulation files can then be created in two ways:
#     - (default) modify the existing MD.in file (to keep comments)
#       e.g.    >>  sim1.createInput()
#     - (optional) create new file with just the key-value pairs. 
#       e.g.    >>  sim2.createInput('filename')
#
# - The restart flag can be set (by default ./results/final_state is used
#     as the restart file):
#       e.g.    >>  sim1.restart = True
#
# - The simulation can then be run.
#     - (default) run simulation
#       e.g.    >>  sim1.run()
#     - (optional) run simulation with specific restart file
#       e.g.    >>  sim1.run('./sheared_final_state')
#
# - There are then various post-processing routines:
#     - Viscosity. Returns array of viscosity average and standard deviation.
#         Currently requires Lees-Edwards boundary conditions and pressure flag.
#         e.g.    >>  print sim1.viscosity()
#                       > [1.67, 0.01] 
#     - Macrolist. Returns array of specific macroscopic property, denoted by 
#         keyword in top of column of 'macroscopic_properties'. Requires 
#         MACRO_OUTFLAG=2 
#         e.g.    >>  sim2.macro()
#                 >>  sim2.run()
#                 >>  steps = sim2.macrolist('iter')
#                 >>  total_energy = sim2.macrolist('TE')
#     - Vslice. Optional boolean arguments to plot the data, and the analytic
#         solution. Requires velocity flag.
#         - (default) return array of arrays of [[spatial_coords], [average v 1] ...] 
#           e.g.    >>  sim1.vslice()
#         - (optional) plot the arrays. Uses analysis class, which uses gnuplot. Plots
#           are saved in './results' directory.
#           e.g.    >>  sim1.vslice(True)
#         - (optional) plot accompanying analytical solution for results.
#           e.g.    >>  sim1.vslice(True, True)
#     - Writesummary. Writes a summary of simulation time and any changed inputs, to 
#         'Simulations.txt' in the script directory. Only inputs changed using specialist
#         methods as opposed to direct alteration of keywords are printed.
#         - (default) Includes viscosity results if calculated.
#           e.g.    >>  sim1.writeSummary()
#         - (optional) Gives user the option of including a custom string
#           e.g.    >>  sim1.writeSummary('some stuff')
#
# - Can disable running of simulations i.e. to re-run a script with just post processing.
#       e.g.    >>  sim1.post = True
#     This will effectively ignore any sim.run() commands. 

class Simulation():
  def __init__(self, sourceDir='./..'):
    # Initialise directories and files
    self.initDir = os.getcwd()
    # Ensure full path is saved
    self.sourceDir = os.path.abspath(sourceDir) 
    
    self.runDir = ''
    self.restart = False
    self.restartFile = self.sourceDir + '/results/final_state'

    self.serial_exe = 'md.exe'
    self.parallel_exe = 'parallel_md.exe'
    self.platform = 'intel'
    self.debug = False

    self.inputFile = 'MD.in'
    self.simOutFile = 'md.out'

    # Don't run any simulations if this is true
    self.post = False

    # Set the default state of all the variables
    # Corresponds to values of default.in (as of 9/10/12)
    self.input={}
    self.input['DENSITY'] = 0.8
    self.input['RCUTOFF'] = 1.12246204830937
    self.input['INPUTTEMPERATURE'] = 1.0
    self.input['INITIAL_CONFIG_FLAG'] = [0, 'special_case_name']
    self.input['SPARSE_FENE'] = [20, 30.0, 1.5, 15.0, 10.0, 10.0, 4]
    self.input['DENSE_FENE'] = [8, 30.0, 1.5, 1.1, 8, 8, 8]
    self.input['INITIALNUNITS'] = [8, 8, 8]
    self.input['NSTEPS'] = 1000
    self.input['TPLOT'] = 100
    self.input['DELTA_T'] = 0.005
    self.input['INTEGRATION_ALGORITHM'] = 0
    self.input['FORCE_LIST'] = 3
    self.input['DELTA_RNEIGHBR'] = 0.3
    self.input['RESCUE_SNAPSHOT_FREQ'] = 21600
    self.input['FIXED_REBUILD_FLAG'] = [0, 20]
    #self.input['SORT_FLAG'] = [2, 1, 3]
    self.input['SORT_FLAG'] = [0, 0, 0]
    self.input['ENSEMBLE'] = 0
    self.input['PERIODIC'] = [1, 1, 1]
    self.input['BFORCE'] = [0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    self.input['DEFINE_SHEAR'] = [1, 0, 0, 1.0]
    self.input['POTENTIAL_FLAG'] = 0
    self.input['FENE_INFO'] = [8, 30.0, 1.5]
    self.input['SOLVENT_INFO'] = [0, 0, 0.0, 0.0, 0.0]
    self.input['SEED'] = [1, 2]
    self.input['WALLSLIDEV'] = [0.0, 0.0, 0.0]
    self.input['SLIDEDISTBOTTOM'] = [0.0, 0.0, 0.0]
    self.input['SLIDEDISTTOP'] = [0.0, 0.0, 0.0]
    self.input['FIXDISTBOTTOM'] = [0.0, 0.0, 0.0]
    self.input['FIXDISTTOP'] = [0.0, 0.0, 0.0]
    self.input['TETHEREDDISTBOTTOM'] = [0.0, 0.0, 0.0]
    self.input['TETHEREDDISTTOP'] = [0.0, 0.0, 0.0]
    self.input['THERMSTAT_FLAG'] = 0
    self.input['THERMSTATBOTTOM'] = [0.0, 0.0, 0.0]
    self.input['THERMSTATTOP'] = [0.0, 0.0, 0.0]
    self.input['CYLINDERS'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0]
    self.input['FILL_CYLINDER'] = [0, 0.0, 0, 0, 0]
    self.input['FREEZE'] = [0.0, 0.0, 10000]
    self.input['PROCESSORS'] = [1, 1, 1]
    self.input['VMD_OUTFLAG'] = [0, 0, "0, 100, 200, 300, 310, 431, 600, 700, 900, 970, 980, 1000"]
    self.input['MACRO_OUTFLAG'] = 1
    self.input['MASS_OUTFLAG'] = [0, 100]
    self.input['VELOCITY_OUTFLAG'] = [0, 100]
    self.input['PRESSURE_OUTFLAG'] = [0, 100]
    self.input['VISCOSITY_OUTFLAG'] = [0, 100]
    self.input['CV_CONSERVE'] = 1
    self.input['MFLUX_OUTFLAG'] = [0, 1]
    self.input['VFLUX_OUTFLAG'] = [0, 1]
    self.input['EFLUX_OUTFLAG'] = [0, 100]
    self.input['ETEVTCF_OUTFLAG'] = [0, 0]
    self.input['R_GYRATION_OUTFLAG'] = [0, 0]
    self.input['RDF_OUTFLAG'] = [0, 3.5, 100]
    self.input['STRUCT_OUTFLAG'] = [0, 1, 2, 50]
    self.changed = {}
    

    # Initialise various simulation properties
    self.finalStep = self.input['NSTEPS']
    self.simProgress = 0
    
    # Initialise post processing values
    self.visc_mu = 0.0
    self.visc_sig = 0.0
    self.Pxy = []
    self.eta = []

    # Simulation run properties
    self.runBG = False
    self.simTime = 0.0

    
   
  ############### ALTERNATIVE INPUT SETTING METHODS #####################
  # i.e. as well as using self.input['NSTEPS'] = 2000, 
  # can use the following methods.

  def density(self, rho):
    self.input['DENSITY'] = rho
    self.changed['DENSITY'] = self.input['DENSITY']
    return

  def cutoff(self, cut):
    self.input['RCUTOFF'] = cut
    self.changed['RCUTOFF'] = self.input['RCUTOFF']
    return

  def temperature(self, t):
    self.input['TEMPERATURE'] = t
    self.changed['TEMPERATURE'] = self.input['TEMPERATURE']
    return

  def cell(self, xunit, yunit, zunit):
    self.input['INITIALNUNITS'] = [xunit, yunit, zunit]
    self.changed['INITIALNUNITS'] = self.input['INITIALNUNITS']
    return

  # Determine cell size from nparticles.
  # Assume cube and round to nearest whole number.
  def nparticles(self, np):
    nc=np/4.0
    n=int(round(math.pow(nc, (1.0/3.0))))
    self.cell(n, n, n)
    self.changed['NPARTICLES'] = n
    return

  def nsteps(self, steps):
    self.input['NSTEPS'] = steps 
    self.changed['NSTEPS'] = self.input['NSTEPS']
    return

  def tplot(self, tp):
    self.input['TPLOT'] = tp
    self.changed['TPLOT'] = self.input['TPLOT']
    return

  def dt(self, delta_t):
    self.input['DELTA_T'] = delta_t
    self.changed['DELTA_T'] = self.input['DELTA_T']
    return

  def ensemble(self, ens):
    self.input['ENSEMBLE'] = ens 
    self.changed['ENSEMBLE'] = self.input['ENSEMBLE']
    return

  def pbc(self, x, y, z):
    self.input['PERIODIC'] = [x, y, z] 
    self.changed['PERIODIC'] = self.input['PERIODIC']
    return


  def shear(self, dirn, itern, rateflag, rate):
    self.input['DEFINE_SHEAR'] = [dirn, itern, rateflag, rate]
    self.changed['DEFINE_SHEAR'] = self.input['DEFINE_SHEAR']
    return


  def potential(self, pot):
    self.input['POTENTIAL_FLAG'] = pot
    self.changed['POTENTIAL_FLAG'] = self.input['POTENTIAL_FLAG']
    return

  def polymer(self, chainl):
    self.input['POTENTIAL_FLAG'] = 1
    self.input['FENE_INFO'][0] = chainl
    self.changed['POTENTIAL_FLAG'] = self.input['POTENTIAL_FLAG']
    self.changed['FENE_INFO'] = self.input['FENE_INFO']
  
  def lj(self):
    self.input['POTENTIAL_FLAG'] = 0
    self.changed['POTENTIAL_FLAG'] = self.input['POTENTIAL_FLAG']


  def solventRatio(self, sr):
    self.input['SOLVENT_INFO'][1] = sr
    if sr == 1:
      self.input['SOLVENT_INFO'][0] = 1
    self.changed['SOLVENT_INFO'] = self.input['SOLVENT_INFO']
    return

  def random(self):
    self.input['SEED'] = [1, 1]
    self.changed['SEED'] = self.input['SEED']
    return

  def freeze(self, initTemp, finalTemp, eqSteps=''):
    self.input['FREEZE'][0] = initTemp
    self.input['FREEZE'][1] = finalTemp
    if eqSteps:
      self.input['FREEZE'][2] = eqSteps
    self.changed['FREEZE'] = self.input['FREEZE']
    return

  def fill(self, on, disp=0.9, thermo=1, tethered=0, wall_thermo=0):
    if on:
      self.input['FILL_CYLINDER'] = [on, disp, thermo, tethered, wall_thermo]
    else:
      self.input['FILL_CYLINDER'][0] = 0
    return

  def nproc(self, nx, ny, nz):
    self.input['PROCESSORS'] = [nx, ny, nz]
    self.changed['PROCESSORS'] = self.input['PROCESSORS']
    return

  def vmd(self, unwrapped=False):
    if not unwrapped:
      self.input['VMD_OUTFLAG'][0]=1
    else: 
      self.input['VMD_OUTFLAG'][0]=4
    self.changed['VMD_OUTFLAG'] = self.input['VMD_OUTFLAG']
    return

  def macro(self, flag=2):
    self.input['MACRO_OUTFLAG']=flag
    self.changed['MACRO_OUTFLAG'] = self.input['MACRO_OUTFLAG']
    return

  def vflag(self, flag, freq):
    self.input['VELOCITY_OUTFLAG'] = [flag, freq]
    self.changed['VELOCITY_OUTFLAG'] = self.input['VELOCITY_OUTFLAG']
    return
  
  def pflag(self, flag, freq):
    self.input['PRESSURE_OUTFLAG'] = [flag, freq]
    self.changed['PRESSURE_OUTFLAG'] = self.input['PRESSURE_OUTFLAG']
    return

  def eteflag(self, flag, freq):
    self.input['ETEVTCF_OUTFLAG'] = [flag, freq]
    self.changed['ETEVTCF_OUTFLAG'] = self.input['ETEVTCF_OUTFLAG']
    return
  ############### /END alternative input setting methods #####################




  
  ############### Calculate certain functional properties #####################
  # Estabilish location of results directory 
  def resultsDir(self):
    if not self.runDir:
      direc = self.sourceDir + '/results' 
      direc = os.path.abspath(direc)
      return direc
    else: 
      direc = self.runDir + '/results' 
      direc = os.path.abspath(direc)
      return direc

  # Calculate number of particles
  def calcNparticles(self):
    np = 4
    for i in self.input['INITIALNUNITS']:
      np *= i
    return np

  # Calculate volume of domain
  def calcDomainVolume(self):
    return self.calcNparticles()/self.input['DENSITY']

  # Calculate length of domain (as array of [Lx, Ly, Lz])
  def calcDomainLength(self):
    nx=self.input['INITIALNUNITS'][0]
    ny=self.input['INITIALNUNITS'][1]
    nz=self.input['INITIALNUNITS'][2]
    n=float(nx*ny*nz)
    V=float(self.calcDomainVolume())
    nV=math.pow(V/n, (1.0/3.0))
    xdomain = nV*nx
    ydomain = nV*ny
    zdomain = nV*nz
    return [xdomain, ydomain, zdomain]

  # Calculate number of processors
  def calcNproc(self):
    np=1
    for i in self.input['PROCESSORS']:
      np *= i
    return np

  # Check if simulation is to be run in serial
  def isSerial(self):
    if self.calcNproc() > 1:
      return False
    else:
      return True

  def exeFile(self):
    if self.isSerial():
      return self.serial_exe
    else:
      return self.parallel_exe

  # Check if Lees-Edwards bc's are to be used
  def le_state(self):
    for i in range(0, 3):
      if self.input['PERIODIC'][i] == 2:
        return True
    return False


  # Return direction of LE 
  def le_dir(self):
    if self.le_state():
      for i in range(0, 3):
        if self.input['PERIODIC'][i] == 2:
          return i
    else:
      print 'Error: LE pbc not enabled'
      return 

  # Work out the Lees-Edwards shear rate
  def le_sr(self):
    if self.le_state():
      sr = self.input['DEFINE_SHEAR'][3] 
      # Check for velocity/rate flag
      if self.input['DEFINE_SHEAR'][2] == 1:
        return sr
      else:
        gamma_dot = sr/self.calcDomainLength()[self.le_dir()]
        return gamma_dot
    else:
      return 0.0
 
  # Work out the Lees-Edwards shear velocity 
  def le_sv(self):
    if self.le_state():
      sr = self.input['DEFINE_SHEAR'][3] 
      # Check for velocity/rate flag
      if self.input['DEFINE_SHEAR'][2] == 0:
        return sr
      else:
        gamma_dot = sr*self.calcDomainLength()[self.le_dir()]
        return gamma_dot
    else:
      return 0.0
  ############### /END functional properties #####################



  ############### SETUP INPUT FILES AND RUN SIMULATION #####################
  # Change line in input file.
  def changeLine(self, variableName, variables):
    os.chdir(self.sourceDir)
    shutil.copy(self.inputFile, self.inputFile+'.bak')
    fid = open(self.inputFile)
    lineno = 1
    status = 0
    while 1:
      tline = fid.readline()
      if not tline:
        status = 1
      if variableName in tline:
        if '#' not in tline:
          lineno += 1
          break
      lineno += 1
      if status != 0:
        print 'Error:', variableName, 'not found in MD.in'
        return
    if type(variables) is list:
      for i in range(0,len(variables)):
        sedstr = "sed '" + str(lineno) + "s/.*/" + str(variables[i]) + "/' <MD.in >new && mv new MD.in" 
        os.system(sedstr)
        lineno += 1
    else:
      sedstr = "sed '" + str(lineno) + "s/.*/" + str(variables) + "/' <MD.in >new && mv new MD.in" 
      os.system(sedstr)
      lineno += 1
    return


  # Create new input file
  def createInput(self, infile=''):
    if self.post:
      return

    if self.runDir:
      os.chdir(self.runDir) # ?
    else:
      os.chdir(self.sourceDir)

    # If no input file is specified then modify existing one (MD.in)
    # i.e. to retain existing comments in .in file.
    if not infile:
      for k, v in self.input.iteritems():
        self.changeLine(k, v)
    # Otherwise create a new .in file with just key-value pairs
    # i.e. no comments explaining what each key represents
    else:
      self.inputFile = infile
      # Remove existing file
      if os.path.exists(self.inputFile):
        os.remove(self.inputFile)
      f = open(self.inputFile, 'w')
      for k, v in self.input.iteritems():
        print>>f, k
        if type(v) is list:
          for i in range(0,len(v)):
            print>>f, v[i]
        else:
          print>>f, v
        print>>f, ''
      f.close()
    return


  # Compile exe if needed
  def compile(self):
    print 'compile'
    if not os.path.exists(self.exeFile()):
      os.system('make clean')
      print 'path exists'
      if self.debug:
        exstring = 'make debug_'
      else:
        exstring = 'make PLATFORM=' + self.platform + ' '
      
      if self.isSerial():
        exstring += 's'
      else:
        exstring += 'p'
      os.system(exstring)
    return


  # Sets up run directory - creates directory if it doesn't exist and
  # copies over input file and source executable.
  def setupRunDir(self):
    self.runDir = os.path.abspath(self.runDir)
    # Check if run and results directories exist and if not, create them
    if not os.path.exists(self.runDir):
      os.makedirs(self.runDir)

    if not os.path.exists(self.resultsDir()):
      os.makedirs(self.resultsDir())

    # Copy over serial and parallel .exe files and input file
    os.chdir(self.sourceDir)
    shutil.copy(self.sourceDir+self.exeFile(), self.runDir+self.exeFile())
    shutil.copy(self.sourceDir+'/MD.in', self.runDir+'/MD.in')
    os.chdir(self.runDir)

    return

     

  # Run simulation
  def run(self, restartFile=''):
    if self.post:
      return
   
    self.compile()

    if self.runDir:
      self.setupRunDir()
      os.chdir(self.runDir)
    else: 
      os.chdir(self.sourceDir)
    
    if restartFile:
      self.restartFile = restartFile
      self.restart = True
    
    outf = ' 2>&1 | tee ' + self.simOutFile
    #outf = ' > ' + self.simOutFile
    
 #   if self.runBG:
 #     outf += ' &'
    
    ctime = time.time()
    
    if self.isSerial():
      if not self.restart:
        #os.system('./md-s.exe' + outf)
        exstring = './' + self.exeFile() + outf
      else:
        # Keep track of total steps following restart. 
        self.finalStep += self.input['NSTEPS']
        #os.system('./md-s.exe -r ' + self.restartFile + outf)
        exstring = './' + self.exeFile() + ' -r ' + self.restartFile + outf
    else:
      if not self.restart:
        #os.system('mpiexec ./md-p.exe -n ' + self.nproc + outf)
        exstring = 'mpiexec -n ' + str(self.calcNproc()) + ' ./' + self.exeFile() + outf
      else:
        # Keep track of total steps following restart. 
        self.finalStep += self.input['NSTEPS']
        #os.system('mpiexec ./md-p.exe -n ' + self.nproc + ' -r ' + self.restartFile + outf)
        exstring = 'mpiexec -n ' + str(self.calcNproc()) + ' ./' + self.exeFile() + ' -r ' + self.restartFile + outf
   
    os.system(exstring)
    self.simTime = time.time() - ctime
    return


  # Check if simulation_progress=nsteps
  def checkProgress(self):
    os.chdir(self.resultsDir)
    f = open('simulation_progress')
    line = f.readline()
    prog = line.split()
    self.simProgress=int(prog[0])
    f.close()
    if self.simProgress == self.finalStep:
      return True
    else:
      return False 


  # Check if simulation has ended from sim output file (for 'Time taken' or 'Error')
  def checkOutfile(self):
    for line in reversed(open(self.simOutFile, 'r').readlines()):
      if 'Time taken' in line:
        return True
      elif 'Error' in line:
        print line.rstrip()
        return False
    return False 
  ### /END input files etc. #############################################   



  ############### POST PROCESSING #####################
  def copyInOut(self, dest):
    if self.runDir:
      os.chdir(self.runDir)
    else:
      os.chdir(self.sourceDir)
    shutil.copy(self.simOutFile, dest)
    shutil.copy(self.inputFile, dest)
    return


  # Copy all files from /results to dest
  def copyResults(self, dest):
    os.chdir(self.resultsDir())
    if not os.path.exists(dest):
      os.makedirs(dest)
    resultFiles = os.listdir('./')

    for file_name in resultFiles:
        full_file_name = os.path.join('.', file_name)
        if (os.path.isfile(full_file_name)):
            shutil.copy(full_file_name, dest)
    self.copyInOut(dest)
    return

  def clearResults(self):
    os.chdir(self.resultsDir())
    resultFiles = os.listdir('./')

    for file_name in resultFiles:
        full_file_name = os.path.join('.', file_name)
        if (os.path.isfile(full_file_name)):
            os.remove(full_file_name)
    return


  def copyRestart(self, dest):
    if self.runDir:
      os.chdir(self.runDir) # ?
    else:
      os.chdir(self.sourceDir)
    shutil.copy(self.restartFile, dest)
    return



  # Calculate viscosity and standard deviation from pressure file
  # Based on D Trevelyan's viscometrics.py
  # Only works for shear in x of y planes atm
  def viscosity(self):
    if self.input['PRESSURE_OUTFLAG'][0] == 1:
      if not self.le_state():
        print 'Error: le boundaries not enabled'
        return
      else:
        os.chdir(self.resultsDir())
        stressfile = 'pvirial'
        if os.path.isfile(stressfile):
          f = open(stressfile, 'rb')
          Ndoubles = os.path.getsize(stressfile)/8
          Nrecords = Ndoubles/9
          data = []
          for i in range(0, Ndoubles):
            data.append(struct.unpack('d', f.read(8)))
          f.close()
          self.Pxy = []
          rec = 0
          for i in range(0, Ndoubles, 9):
            self.Pxy.append(data[i+3])
            rec += 1
          self.eta = []
          for rec in range(Nrecords):
            etaval = -self.Pxy[rec][0]/self.le_sr()
            self.eta.append(etaval)
          self.visc_mu = average(self.eta)
          self.visc_sig = sigma(self.eta, self.visc_mu)
          self.changed['VISCOSITY'] = str(self.visc_mu) + ' +/- ' + str(self.visc_sig) 
          return [self.visc_mu, self.visc_sig]
        else: 
          print 'Error: cannot find pvirial'
          return

    else: 
      print 'Error: pressure outflag != 1'
      return


  # Returns output of macroscopic properties column defined by varName
  # as a list.
  #   e.g. 'iter', 'Temp', 'KE', 'PE', 'TE', 'Pressure'
  def macrolist(self, varName):
    os.chdir(self.resultsDir())
    f = open('./macroscopic_properties', 'r')
    titleline = f.readline()
    titleline = titleline.replace('\n', '')
    titleline = titleline.replace(' ', '')
    #titleline.remove('\n')
    titleline = titleline.split(';')
    elem = False
    for i in range(0, len(titleline)):
      if titleline[i] == varName:
        elem = i
        break
    #f.close()
    if not elem:
      print 'Error: variable name not found in macroscopic properties'
      return
    varlist = [] 
    for line in f.readlines():
      currentline = line
      currentline = currentline.replace('\n', '')
      currentline = currentline.split(';')
      #currentline.remove('\n')
      varlist.append(float(currentline[elem]))
    return varlist


  def read_header(self, varName):
    os.chdir(self.resultsDir())
    f = open('simulation_header','r')    # Create ascii file object
    headerdata = f.readlines()                      # Read all lines into list
    for i in range(len(headerdata)):                # Loop all lines 
        headerdata[i] = headerdata[i].split(';')    # Split list elements into
                                                    # sub-lists by semicolon 
        for j in range(len(headerdata[i])):         # Loop through sublist
            headerdata[i][j] = (headerdata[i][j].
                                strip().replace(
                                '(','').replace(
                                ')',''))            # Strip white space and ()
        vars()[headerdata[i][1]] = headerdata[i][2] # Create variables with name 
        if varName == headerdata[i][1]:
          f.close()
          return headerdata[i][2]                   # headerdata[i][1] and value
                                                    # headerdata[i][2]

    f = open('simulation_progress','r')  # Create ascii file object
    Nsteps = f.read().strip()                       # Read number of timesteps 
    f.close()


  def mslice(self, disp=0):
    from array import array
    if not np:
      print 'Error: numpy is not installed. mslice() will not work'
      return
    # Retrieve important variables from simulation_header with read_header.py
    Nsteps        = int(self.read_header('Nsteps'))
    tplot         = int(self.read_header('tplot'))
    Nmass_ave     = int(self.read_header('Nmass_ave'))
    initialstep   = int(self.read_header('initialstep'))
    m_outflag     = int(self.read_header('mass_outflag'))
    gnbins        =[int(self.read_header('gnbins1')),
                        int(self.read_header('gnbins2')),
                                        int(self.read_header('gnbins3'))]
    Nmass_records = int(math.floor((Nsteps-initialstep+disp)/(tplot*Nmass_ave)))
    nbins         = int(gnbins[m_outflag-1])

    os.chdir(self.resultsDir())
    f = open('mslice','rb')             # Create binary file object
    mslice = array('i')                            # Initialise array of integers
    mslice.fromfile(f,nbins*Nmass_records)         # Read binary file into array
    mslice = reshape(mslice,(Nmass_records,nbins)) # Reshape array into 3D
    f.close()                                      # Close file object
    return mslice


  def vslice_analytic(self, y, t, mu): 
    pi = math.pi
    mr = 2.061/self.input['DENSITY']
    u0 = self.le_sv()*0.5
    vdirn = self.input['VELOCITY_OUTFLAG'][0]-1
    Y = self.calcDomainLength()[vdirn]*0.5
    total = y*u0/Y
    for n in range(1, 10000):
      l = (n*pi/Y)
      l2 = l*l 
      un = 2*u0*pow(-1, n)/(n*pi)
      un *= (math.exp(-l*l*mr*t))
      un *= math.sin(l*y)
      total += un
    return total


  # Time average of vslice_analytic between t0 and tf, using nsteps
  def vslice_analytic_range(self, y, t0, tf, mu, nsteps):
    dt = float(tf-t0)/float(nsteps-1)
    t = t0
    total = 0.0
    for i in range(nsteps):
      total += self.vslice_analytic(y, t, mu)
      t += dt
    av = total/float(nsteps)
    return av


  def vslice(self, plot=False, analytic=False):
    from array import array

    # Check if numpy is installed
    if not np:
      print 'Error: numpy is not installed. vslice() will not work'
      return

    # Averages a list of values
    def mean(values):                                 
        return float(sum(values)/len(values))

    # Initialise simulation constants
    nd            = 3
    Nsteps        = int(self.read_header('Nsteps'))
    tplot         = int(self.read_header('tplot'))
    Nvel_ave      = int(self.read_header('Nvel_ave'))
    lesd          = int(self.read_header('le_sd')) - 1
    initialstep   = int(self.read_header('initialstep'))      # Initial timestep
    v_outflag     = int(self.read_header('velocity_outflag')) # Slice direction
    gnbins        =[int(self.read_header('gnbins1')),         # Number of bins in each
                    int(self.read_header('gnbins2')),         # direction
                    int(self.read_header('gnbins3'))]
    Nvel_records  = int(math.floor(Nsteps-initialstep)/(tplot*Nvel_ave))
    nbins         = gnbins[v_outflag-1]               # Python counts from 0

    vdirn = self.input['VELOCITY_OUTFLAG'][0]-1
    Y = self.calcDomainLength()[vdirn]

    # Generate mass slice
    os.chdir(self.resultsDir())
    mslice = self.mslice()
    
    f = open('../results/vslice','rb')                # File object f
    vslice = array('d')                               # Array to store data
    vslice.fromfile(f,Nvel_records*nd*nbins)          # Get vslice data
    vslice = reshape(vslice,(Nvel_records,nd,nbins))  # Reshape array
    f.close()

    nt        = 16                  # Split vslice into nt time regions
    vprofile  = [[0.0]*nbins]*nt    # Initialise array of average velocity profiles
    sp_coords = [0.0]*nbins         # Initialise shear plane coordinates
    plotrange = [0, 1, 3, 7, 15]    # Plot in powers of 2

    # Create array of arrays to store values for post-processing.
    results = []
    for i in range(len(plotrange)+1):
      results.append([])
    count = 1

    for t in plotrange:
      for cell in range(nbins):                     # Loop over nbins cells
        sp_coords[cell] = float(cell+0.5)/nbins   # Center of cell coordinates
        varray = vslice[t*(Nvel_records/nt)+1:(t+1)*(Nvel_records/nt),
                        lesd,cell]               # Array of velocities for
                                                  # this cell over Nvel_records/
                                                  # nt timesteps
        marray = mslice[t*(Nvel_records/nt)+1:(t+1)*(Nvel_records/nt),
                        cell]                     # Total mass over Nvel_records
                                                  # /nt timsteps
        vprofile[t][cell] = mean(varray/marray)   # Average velocity profile
        outstring = (str(sp_coords[cell]).rjust(16) + str(vprofile[t][cell]).rjust(32))
        results[count].append(vprofile[t][cell])
      count += 1

    # Set the first array of arrays to the spatial co-ordinates, y/Y, between -0.5 and 0.5
    for sp in sp_coords:
      y = sp-0.5
      results[0].append(y)

    # If the plot option is not specified, simply return the array of arrays.
    if not plot: 
      return results

    # Otherwise plot ... 
    else:
      label = ['u_x', 'y/Y']
      results_legend = [] 
      if analytic:
        analytic_legend = []

      t = 0.0

      # Generate appropriate data set titles e.g. '0 - 3', '3 - 6' etc.
      for i in range(len(plotrange)):
        t = plotrange[i]*self.input['NSTEPS']*self.input['DELTA_T']/float(nt)
        dt = self.input['NSTEPS']*self.input['DELTA_T']/float(nt)
        lstr = str(int(math.floor(t))) + ' - ' + str(int(math.floor(t+dt)))
        results_legend.append(lstr)
        if analytic:
          analytic_legend.append(lstr + str(' analytic'))

      # Add in the analytic solution if required.
      if analytic:
        analytic = []
        for i in range(len(plotrange)+1):
          analytic.append([])

        visc = 1.674

        # Populate spatial co-ordinates > 0 
        for yOverY in results[0]:
          if yOverY > 0:
            analytic[0].append(yOverY)

        for i in range(0, len(plotrange)):
          # Create new array for analytic solution at time t.
          t = plotrange[i]*self.input['NSTEPS']*self.input['DELTA_T']/float(nt)
          dt = self.input['NSTEPS']*self.input['DELTA_T']/float(nt)

          # Populate using spatial values (results[0]) at time midpoint
          for yOverY in analytic[0]:
            y = Y*yOverY
            analytic[i+1].append(self.vslice_analytic_range(y, t, t+dt, visc, 10))
          
      
      # Generate analyis object for plotting. 
      p = Analysis(results, label, results_legend, 'vslice.pdf', [], [-0.5, 0.5])
      p.array2 = analytic
      p.legends2 = analytic_legend

      # Plot (y,x) instead of (x,y)
      p.reverse = True
      p.reverse2 = True

      # Plot the second half of the results with lines instead of points

      p.plot()

      return


  def vslice_radial(self, disp, plotrange):
    from array import array

    # Check if numpy is installed
    if not np:
      print 'Error: numpy is not installed. vslice() will not work'
      return

    # Averages a list of values
    def mean(values):                                 
        return float(sum(values)/len(values))

    # Initialise simulation constants
    nd            = 3
    Nsteps        = int(self.read_header('Nsteps'))
    tplot         = int(self.read_header('tplot'))
    Nvel_ave      = int(self.read_header('Nvel_ave'))
    lesd          = int(self.read_header('le_sd')) - 1
    initialstep   = int(self.read_header('initialstep'))      # Initial timestep
    v_outflag     = int(self.read_header('velocity_outflag')) # Slice direction
    gnbins        =[int(self.read_header('gnbins1')),         # Number of bins in each
                    int(self.read_header('gnbins2')),         # direction
                    int(self.read_header('gnbins3'))]
    Nvel_records  = int(math.floor(Nsteps-initialstep+disp)/(tplot*Nvel_ave))
    nbins         = gnbins[v_outflag-1]               # Python counts from 0

    vdirn = self.input['VELOCITY_OUTFLAG'][0]-1
    Y = self.calcDomainLength()[vdirn]

    # Generate mass slice
    os.chdir(self.resultsDir())
    mslice = self.mslice(disp)
   
    f = open('../results/vslice','rb')                # File object f
    vslice = array('d')                               # Array to store data
    vslice.fromfile(f,Nvel_records*nd*nbins)          # Get vslice data
    vslice = reshape(vslice,(Nvel_records,nd,nbins))  # Reshape array
    f.close()

    nt        = Nvel_records
    vprofile  = [[0.0]*nbins]*nt    # Initialise array of average velocity profiles
    sp_coords = [0.0]*nbins         # Initialise shear plane coordinates

   # print 'lesd = ', lesd
   # print 'Nvel_records = ', Nvel_records
   # print 'Nsteps = ', Nsteps
   # print 'init step = ', initialstep
   # print 'tplot = ', tplot
   # print 'nvelav = ', Nvel_ave
   # print 'nd', nd

    # Create array of arrays to store values for post-processing.
    results = []
    for i in range(len(plotrange)+1):
      results.append([])
    count = 1

   # print 'nbins = ', nbins
   # print 'plotrange = ', plotrange
   # print 'mslice = ', mslice
   # print 'vslice = ', vslice
    for t in [0]:
      for cell in range(nbins):                     # Loop over nbins cells
        sp_coords[cell] = float(cell+0.5)/nbins   # Center of cell coordinates
        varray = vslice[t,lesd,cell]
        marray = mslice[t,cell]
        vprofile[t][cell] = varray/marray   # Average velocity profile
        outstring = (str(sp_coords[cell]).rjust(16) + str(vprofile[t][cell]).rjust(32))
        results[count].append(vprofile[t][cell])
      count += 1

    # Set the first array of arrays to the spatial co-ordinates from 0 to Y
    for sp in sp_coords:
      y = sp*Y
      results[0].append(y)

    return results

  def etev(self, justY=False, plot=False):
    import os
    import sys 
    import numpy as np
    import scipy.signal as sp

    # Get important information from simulation_header with read_header.py
    delta_t = float(self.read_header('delta_t'))
    tplot   = int(self.read_header('tplot'))
    Nsteps  = int(self.read_header('Nsteps'))
    nchains = int(self.read_header('nchains'))

    os.chdir(self.resultsDir())
    fileloc = 'etev'                   # File location

    etev   = np.fromfile(fileloc,dtype=float)     # Read from binary file
    Ndt    = len(etev)/(3*nchains)                # Find number of timesteps
    etev_x = np.empty((nchains,Ndt),float)    
    etev_y = np.empty((nchains,Ndt),float)
    etev_z = np.empty((nchains,Ndt),float)

    # Begin rearrange of etev
    for i in range(nchains):
      etev_x[i]  = etev[0+i::3*nchains]
      etev_y[i]  = etev[nchains+i::3*nchains]
      etev_z[i]  = etev[2*nchains+i::3*nchains]
    del(etev)
    # End rearrange

    # Calculate auto-correlation functions
    C = np.zeros(Ndt)
    for i in range(nchains):
      auto_cx = sp.fftconvolve(etev_x[i],etev_x[i][::-1]) # Correlation is 
                                                          # convolution in reverse
      auto_cy = sp.fftconvolve(etev_y[i],etev_y[i][::-1]) 
      auto_cz = sp.fftconvolve(etev_z[i],etev_z[i][::-1]) 
      l = len(auto_cx)                                    # Only positive half
      C += auto_cx[l/2:l]+auto_cy[l/2:l]+auto_cz[l/2:l]

    C = C/float(nchains)                                    # Ensemble average
    C = C/C[0]                                              # Normalise

    eteResults = [[], []]
    for i in range(len(C)):    
      eteResults[0].append(str(i*tplot*delta_t)) 
      eteResults[1].append(str(C[i]))

    if not plot:
      if not justY:
        return eteResults
      else:
        return eteResults[1]
    else:
      label = ['t', 'C(t)']
      legend = ['']
      p = Analysis(eteResults, label, legend, 'etev.pdf', [], [0, 1], 'x')
      p.plot()
      return

 



  def writeSummary(self, optstring=''):
    os.chdir(self.initDir)
    f = open('Simulations.txt', 'a')
    print>>f, '###########################################'
    print>>f, 'Date: \t', datetime.datetime.now()
    print>>f, 'Simulation duration: \t', self.simTime
    for k,v in self.changed.iteritems():
      outstr = str(k) + ': \t' + str(v)
      print>>f, outstr
    if optstring:
      print>>f, optstring
    print>>f, '\n'
    print>>f, '\n'
    f.close()

    return

# END OF SIMULATION CLASS
#######################################################################################






# ----------------------------------------------------------------------------
# Class Analysis 
#
# Encapsulates the process of writing arrays to file and plotting them using gnuplot.
# 
# - Allows the creation of analysis objects. On initialisation: 
#     - (default) supply array of arrays of data. Assumes first array contains x co-ords
#       and the other n >= 1 arrays contain the corresponding sets of y co-ords.
#       e.g.    >>  vslice_analysis = Analysis(vslice_results)
#     - (optional) specify array of x/y labels as strings
#       e.g.    >>  vslice_analysis = Analysis(vslice_results, ['xlabel', 'ylabel'])
#     - (optional) specify array of plot legends as strings. len(legends) must equal n.
#       e.g.    >>  legends = ['data 1', 'data 2', data 3']
#               >>  vslice_analysis = Analysis(vslice_results, labels, legends)
#     - (optional) specify output file  
#       e.g.    >>  vslice_analysis = Analysis(vslice_results, labels, legends, 'vslice.pdf')
#     - (optional) specify xrange and yrange
#       e.g.    >>  vslice_analysis = Analysis(vslice_results, labels, legends, 'vslice.pdf', 
#                                                 [0, 3.141], [-1,1])
#     - (optional) specify if x and or y axis should be plotted with a logscale
#       e.g.    >>  vslice_analysis = Analysis(vslice_results, labels, legends, 'vslice.pdf', 
#                                               rangeX, rangeY, 'x')
#               >>  vslice_analysis = Analysis(vslice_results, labels, legends, 'vslice.pdf', 
#                                               rangeX, rangeY, 'xy')
#
# - Creates array of plotstyles, assuming all are to be plotted as points ("p"). 
#
# - All of the above can also be modified after the creation of the object 
#   e.g.    >>  vslice_analysis = Analysis(vslice_results)
#           >>  vslice_analysis.arrays = new_vslice_results
#           >>  vslice_analysis.labels = new_labels
#           >>  vslice_analysis.legends = new_legends
#           >>  vslice_analysis.rangeX = [0, 2*3.141]
#           >>  vslice_analysis.rangeY = [0, 1] 
#           >>  vslice_analysis.log = 'y' 
#           >>  vslice_analysis.plotStyle[2] = "lp" 
#           >>  vslice_analysis. 
#
# - Can set keyword 'reverse' to plot (y,x) instead of (x,y)
#   e.g.    >>  vslice_analysis.reverse = True
#
# - There is a method to write the array of arrays to file in column format:
#   e.g.    >>  vslice_analysis.writeArr('data.txt')
#
# - There is a plot method which first uses writeArr to write to plotfile.data, then
#     creates a gnuplot script using any labels/legends/ranges/plotstyles and then executes
#     the gnuplot script.
#   e.g.    >>  vslice_analysis.plot()

class Analysis():
  def __init__(self, array, labels=[], legends=[], plotfile='analysis.pdf', rangeX=[], rangeY=[], log=''):
    # fn plot() to accept this.
    self.array = array
    self.labels = labels
    self.legends = legends
    self.plotfile = plotfile
    self.rangeX = rangeX 
    self.rangeY = rangeY
    self.log = log
    
    self.plotStyle = []
    for i in range(len(self.array)):
      self.plotStyle.append("p")

    # Plot order is reversed i.e. plot blah u 2:1, 3:1 ... instead of 
    # 1:2, 1:3 etc. Useful for e.g. vslice.  
    self.reverse = False

    # Specify a second set of results
    self.array2 = []
    self.legends2 = []
    self.plotStyle2 = []
    self.reverse2 = False

  # Write array of arrays to file in column format 
  # Assume the array to be written is self.array unless otherwise specified.
  def writeArr(self, outfile, arr=[]):
    if not arr:
      array = self.array
    else:
      array = arr

    g = open(outfile, 'w')
    num_arr = len(array)
    len_arr = len(array[0])  
    count = 0
    for arr in array:
      if len(arr) != len_arr:
        print 'Error: length of arrays dont match'
        print count, len(arr), len_arr
        return
      count += 1

    for i in range(0, len_arr):
      out_str = ''
      for j in range(0, num_arr):
        out_str += str(array[j][i])
        out_str += '\t'
      print>>g,out_str
    g.close()
    return


  # Creates plot from list of lists using gnuplot (must be installed)
  # Plots first list as x
  # Plots remaining lists as y_1, y_2 ... 
  # x/y labels determined by list 'labels' and 'legend'
  # x/y ranges can be specified by lists i.e. [0, 10] 
  # log axes in x and or y can also be enabled with a string e.g. 'x', 'xy' 
  def plot(self):
    datafile = self.plotfile + ".data"
    datafile2 = self.plotfile + ".data2"
    gpfile = self.plotfile + ".gnu"
    num_arr = len(self.array)
    num_arr2 = len(self.array2)

    # Set blank legend titles if none specified.
    if not self.legends:
      for i in range(num_arr - 1):
        self.legends.append('')


    if len(self.labels) != 2:
      print 'Error: 2 plot titles required i.e. for x and y'
      return
    
    if len(self.legends) != num_arr - 1:
      print 'Error: not enough legend titles'
      print 'num_arr', num_arr
      print 'num_leg', len(self.legends)
      return

    # Write arrays to temporary files for plotting
    self.writeArr(datafile)
    if self.array2:
      self.writeArr(datafile2, self.array2)
  
    gpout = []
    gpout.append("set term pdf enhanced")
    gpout.append("set output '" + self.plotfile + "'")
    
    if self.labels:
      gpout.append("set xlabel '" + self.labels[0] + "'") 
      gpout.append("set ylabel '" + self.labels[1]+ "'")
    
    if self.log:
      gpout.append("set log " + self.log)

    if self.rangeX:
      gpstring = "set xrange [ " + str(self.rangeX[0]) + " : " + str(self.rangeX[1]) + " ] " 
      gpout.append(gpstring)

    if self.rangeY:
      gpstring = "set yrange [ " + str(self.rangeY[0]) + " : " + str(self.rangeY[1]) + " ] " 
      gpout.append(gpstring)

    gpout.append("set xtics auto; set ytics auto")
    gpout.append("a='" + datafile + "'")

    # Define line styles
    for i in range(1, num_arr-1):
      n=str(i)
      gpstring = "set style line "+n+" lt "+n+" pt "+n
      gpout.append(gpstring)

    # Generate plot string (plot 'blah.txt' u 1:2, ... etc)
    plotstring = "plot "
    for i in range (0, num_arr-1):
      # DEPRECATED - IGNORE ::: 
      # Plot lines in same colours as points
      # ie if there are 6 sets of y data (+ the x axis), where the last three are 
      # to be plotted as lines, then match styles e.g. 1:4, 2:5, 3:6.  
      #if self.plotStyle[i] == "l":
      #  ylen = int((len(self.array)-1)/2)
      #else:
      ######### IGNORE ####
      ylen = 0
      lsn = str(i+1-ylen)

      if not self.reverse:
        plotstring += "a u 1:" + str(i+2) + " w " + self.plotStyle[i] + " ls " + lsn + " t '" + self.legends[i] + "'"
      else:
        plotstring += "a u " + str(i+2) + ":1 w " + self.plotStyle[i] + " ls " + lsn + " t '" + self.legends[i] + "'"
      if i < num_arr-2:
        plotstring += ", "

    # If there is a second array to be plotted, add that to the script too.
    if self.array2:
      # If no plot style is specified for this array, assume lines are to be used.
      if not self.plotStyle2:
        for i in range(len(self.array2)):
          self.plotStyle2.append("l")

      for i in range (0, num_arr2-1):
        ylen = 0
        lsn = str(i+1-ylen)
        gpout.append("b='" + datafile2 + "'")
        if not self.reverse2:
          plotstring += ", b u 1:" + str(i+2) + " w " + self.plotStyle2[i] + " ls " + lsn + " t '" + self.legends2[i] + "'"
        else:
          plotstring += ", b u " + str(i+2) + ":1 w " + self.plotStyle2[i] + " ls " + lsn + " t '" + self.legends2[i] + "'"
    gpout.append(plotstring)

    # Write gnuplot script from gpout array and execute
    f = open(gpfile, 'w')
    for k in gpout:
      print>>f,k
    f.close()
    plotcmd = "gnuplot " + gpfile
    os.system(plotcmd)

    return

# END OF ANALYSIS CLASS
#######################################################################################
