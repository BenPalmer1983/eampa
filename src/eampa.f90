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
  Use prepdl     		! prep dlpoly files
  Use preppotfit     	! prep potfit files
  Use run			    ! calc
  Use pwbatch			! calc
  Use output			! output
  Use clean
  
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

!prepare dlpoly files
  If(optionRunDLPrep.eq.1)Then 
    Call runPrepdl()
  End If
  
!prepare potfit files 
  If(prepPotfitOption(1:1).eq."Y".or.optionRunPrepPotfit.eq.1)Then
    Call runPrepPotfit()
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

  
!Clean up  
  Call runClean()
  
  
  
  
End