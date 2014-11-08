Program eampa
! University of Birmingham
! Ben Palmer
!
! This package has been designed with several aims in mind.
! 1. make batch of input files for PWscf
! 2. Analyse EAM potentials
! 3. Fit polynomial splines to EAM functions and output in tabulated form
! 4. Fit these functions to DFT generated data
! 5. Predicted bulk properties from EAM potentials
!
! Internal units:
! - energy, eV
! - length, Angstrom
! - forces, ev/Angstrom
!
!Setup Modules
  Use kinds          ! data kinds
  Use msubs           ! mpi module
  Use constants      ! physical constants module
  Use units          ! unit conversion and normalisation 
  Use general        ! string functions
  Use maths          ! maths functions
  Use globals        ! declare all globals
  Use initialise     ! initialise program
  Use loadData       ! load important data
  Use readinput      ! read input
  Use readEAM        ! read EAM potential file  
  Use prepEAM        ! read EAM potential file  
  Use readConfig     ! read config file  
  Use neighbourList  ! make neighbour list 
  Use prep           ! prepare before calculations etc 
  Use calcEAM
  Use calcEval
  Use optimise
  Use testEAM
  Use pwBatch        ! read config file  
  Use output
  Use clean
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'  
!Variables
  Integer(kind=StandardInteger) :: error    
!store start time
  Call cpu_time(programStartTime)
!-------------------------------------------------------------- 
!--- Load data and initialise
!Init MPI
  Call MPI_Init(error)   
! Run globals module:
! Initialise the default values for the global variables
  Call initGlobals()
! Run initialisation module:
! Make and store output/temp directories
! Create output files
! Init a file cleanup list
  Call runInitialise()  
! Run load data module:
! Loads isotope data into 4 arrays
! Any other data useful should be loaded here
  Call loadIsotopeData()
!-------------------------------------------------------------- 
!--- Read Input Files
! Read user input file:
  Call readUserInput()
! Optional read ins
  If(optionReadEAM.eq.1)Then
    Call readEAMFile()
    Call runPrepEAM()
  End If
  If(optionReadConf.eq.1)Then
    Call readConfigFile()
  End If
  If(optionNeighbourList.eq.1)Then
    Call makeNeighbourList()
  End If
! Prepare for calculations
  Call runPrep()
!-------------------------------------------------------------- 
!--- Eval/Calculate config energies
  If(optionCalcEnergies.eq.1)Then
    Call calcEnergies()
  End If  
  If(optionEval.eq.1)Then
    Call evaluate()
  End If
  If(optionEvalFull.eq.1)Then
    Call evaluate()
  End If
  
!-------------------------------------------------------------- 
!--- Optimise input EAM potential functions
  If(optionOptimise.eq.1)Then
    Call runOptimise()
  End If  
  
!-------------------------------------------------------------- 
!--- Test input EAM potential functions  
  If(optionTestEAM.eq.1)Then
    Call runTestAnalysis()
  End If  
  
!-------------------------------------------------------------- 
!--- PW Batch Files
  If(optionRunPWBatch.eq.1)Then
    Call runPWBatch()
  End If
  
  
  
!--------------------------------------------------------------   
!--- Clean up and finalise 
  Call runClean()
  
  
! Store end time
  Call cpu_time(programEndTime)
! Store Time    
  Call storeTime(100,programEndTime-programStartTime)    
! Output times
  Call outputCpuTimes()  
! Call end output to terminal 
  Call outputEndT()
! Finalise MPI
  Call MPI_Finalize(error)
    
  
End Program eampa