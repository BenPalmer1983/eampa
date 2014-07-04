Module datagen

!--------------------------------------------------------------!
! Prep EAM File Subroutines                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 1st May 2014
!----------------------------------------

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use output
  Use prep

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!Privacy of functions/subroutines/variables
  Private    
!Public subroutines
  Public :: runDatagen
  
  
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Subroutines
!-------------------------------------------------------------------------  
! 

  Subroutine runDatagen() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
    If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Run datagen module"
	End If
	
	
  End Subroutine runDatagen
  
 
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 

  
  
End Module datagen 
