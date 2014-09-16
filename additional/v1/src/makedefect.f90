Module makedefect

!--------------------------------------------------------------!
! Make Defect Configs Subroutines                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 6th August 2014
!----------------------------------------

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!Privacy of functions/subroutines/variables
  Private    
!Public subroutines
  Public :: makeDefectConfigs
  
  
  
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

  Subroutine makeDefectConfigs() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
    If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  !print *,ProgramTime(),"Run datagen module"
	End If
	
	
  End Subroutine makeDefectConfigs
  
 
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 

  
  
End Module makedefect 
