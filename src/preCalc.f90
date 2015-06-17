Module preCalc

! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: preCalc
! Updated: 18th May 2015
! --------------------------------------------------------------!
! Description:
! Sits between reading input files and calculation
! --------------------------------------------------------------!


! Setup Modules
  Use kinds
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use globals        ! declare all globals
  Use initialise     ! initialise program
  Use output
  
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Subroutines
  Public :: runPreCalc        !Subroutine

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

  contains

! Run all the input subroutines
  Subroutine runPreCalc()
! Internal subroutine variables
    Integer(kind=StandardInteger) :: configID

   
! Assign config to process map
    Do configID=1,configCount    
      processMap(configID) = mod(configID-1,mpiProcessCount)      
    End Do
    Call outputProcessMapT()
    
    
! Synch MPI processes
    Call M_synchProcesses()    
  End Subroutine
  
  
  
  
  
  
  
  
  
  
End Module preCalc