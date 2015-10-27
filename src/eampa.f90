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
! - stress, ev/ang^3
!
! Printed units:
! - energy, eV
! - length, Angstrom
! - forces, ev/Ang
! - stress, GPa
!
! Setup Modules
  Use libBP          ! data kinds
  Use types
  Use msubs          ! mpi module
  Use typesM
  Use globals        ! declare all globals
  Use initialise     ! initialise program
  Use loadData       ! load important data
  Use readinput      ! read input
  Use readEAM        ! read eam file
  Use makePotential
  Use readConfig     ! read user/dft configurations
  Use bpConfig
  Use neighbourList  ! build neighbour list
  Use neighbourListBP  ! build neighbour list for bulk property configs
  Use preCalc
  Use calcEAM
  Use bpCalcEAM
  Use eval
  Use opti
  Use evalBP
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
!  Set up MPI Derived Data Types
!  Call SetUpMPITypes()  
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
! -----------------------------------------------------------------------  
! --- Run Type: EVAL
! -----------------------------------------------------------------------
  If(eampaRunType.eq."EVAL")Then
! Read the EAM functions/functionals into memory
    Call readEAMFile()             ! readEAM.f90 
! Read config file, dft output files and save these as one config file
    Call readConfigFile()          ! readConfig.f90 
! Read bulk property configs
    Call readBpConfigFile()        ! readConfig.f90
! Prepare a 256 FCC config and 128 BCC config for bulk property testing
    Call prepareBPConfig()         ! bpConfig.f90
! Save input files
    Call outputInputFiles()        ! output.f90
! Make neighbour list
    Call makeNeighbourList()       ! neighbourList.f90
! Make bulk property neighbour list    
    Call makeNeighbourListBP()     ! neighbourListBP.f90
! Pre calculation
    Call runPreCalc()              ! preCalc.f90
! Pre calc summary
    Call outputPreCalcSummaryT()   ! output.f90
! Calculate stress/energy/force of configuration/s
    quietOverride = .false.
    Call evalEAM()                 ! eval.f90
! Bulk Properties 
    Call evalBulkProperties(.false.,.true.)      ! evalBP.f90
! Output summary on input EAM potential    
    Call outputEAMSummaryT()
! Finalise
    Call runFinaliseEval()         ! finalise.f90  prints out run time data
  End If
! -----------------------------------------------------------------------  
! --- Run Type: OPTI
! -----------------------------------------------------------------------
  If(eampaRunType.eq."OPTI")Then
  print *,"opti"
! Read the EAM functions/functionals into memory
    Call readEAMFile()            ! readEAM.f90
! Read config file, dft output files and save these as one config file
    Call readConfigFile()         ! readConfig.f90
! Read bulk property configs
    Call readBpConfigFile()       ! readConfig.f90
! Prepare a 256 FCC config and 128 BCC config for bulk property testing
    Call prepareBPConfig()        ! bpConfig.f90
! Save input files
    Call outputInputFiles()
! Make neighbour list
    Call makeNeighbourList()
! Make bulk property neighbour list    
    Call makeNeighbourListBP()
! Pre calculation
    Call runPreCalc()
! Pre calc summary
    Call outputPreCalcSummaryT()
! Run optimise
    Call optiEAM()
! Finalise
    Call runFinaliseEval()         ! finalise.f90  prints out run time data
  End If
! -----------------------------------------------------------------------  
! --- END
! -----------------------------------------------------------------------
! Store end time
  Call cpu_time(programEndTime)
! Finalise MPI
  Call MPI_Finalize(error)

End Program eampa
