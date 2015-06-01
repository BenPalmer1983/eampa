PROGRAM eampa
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
! Setup Modules
  Use kinds          ! data kinds
  Use msubs          ! mpi module
  Use constants      ! physical constants module
  Use units          ! unit conversion and normalisation
  Use general        ! string functions
  Use maths          ! maths functions
  Use mMaths         ! mpi maths functions
  Use globals        ! declare all globals
  Use initialise     ! initialise program
  Use loadData       ! load important data
  Use readinput      ! read input
  Use readEAM        ! read eam file
  Use readConfig     ! read user/dft configurations
  Use neighbourList  ! build neighbour list
  Use preCalc
  Use calcEAM
  Use output
  Use finalise
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Variables
  Integer(kind=StandardInteger) :: error
! store start time
  Call cpu_time(programStartTime)
! --------------------------------------------------------------
! --- Load data and initialise
! Init MPI
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
! --------------------------------------------------------------
! --- Read Input Files
! Read user input file into memory
  Call readUserInput()
! Read the EAM functions/functionals into memory  
  Call readEAMFile()
! Read config file, dft output files and save these as one config file  
  Call readConfigFile()
! Save input files
  Call outputInputFiles()
! Make neighbour list
  Call makeNeighbourList()
! Pre calculation  
  Call runPreCalc()
  
  Call calcEnergies()
  
  Call runFinalise()

! Store end time
  Call cpu_time(programEndTime)
  
! Finalise MPI
  Call MPI_Finalize(error)

End Program eampa
