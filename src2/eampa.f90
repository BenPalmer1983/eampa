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
!  Use globals        ! declare all globals
!  Use initialise     ! initialise program

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Variables
  Integer(kind=StandardInteger) :: error
! --------------------------------------------------------------
! --- Load data and initialise
! Init MPI
  Call MPI_Init(error)
! Run globals module:
! Initialise the default values for the global variables
!  Call initGlobals()

! Finalise MPI
  Call MPI_Finalize(error)

End Program eampa
