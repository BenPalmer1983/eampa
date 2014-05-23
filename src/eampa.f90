Program eampa

! University of Birmingham
! Ben Palmer
!

!Setup Modules
  Use kinds				!data kinds
  Use constants			!physical constants module
  Use units				!unit conversion and normalisation 
  Use strings		    !string functions
  Use maths				!maths functions
  Use initialise		! initialise program
  Use input				! input
  Use prep     			! prep module
  Use calc			    ! calc
  Use prepeam			! calc
  Use run			    ! calc
  Use pwbatch			! calc
  Use output			! output
  
!Include MPI header
  Include 'mpif.h'
  
!Variables
  Integer(kind=StandardInteger) :: error  
  



!run initialisation module
  Call runInitialise()

!read input files, potential and configurations
  Call runInput()  

!prepare data, make neighbour list
  If(optionRunPrep.eq.1)Then
    Call runPrep()
  End If

!prepare eam potential
  If(optionRunPrepEAM.eq.1)Then 
    Call runPrepeam()
  End If

!start calculations
  If(optionRunProcesses.eq.1)Then
    Call runProcesses()
  End If
  
 !Run pwscf batch file
  If(optionRunPWBatch.eq.1)Then
    Call runPWBatch()
  End If
  
!Finalise MPI
  Call MPI_Finalize(error)

  
  
  
  
  
  
  
End