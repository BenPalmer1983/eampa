Module prepeam

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
  Public :: runPrepeam
  
  
  
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


  Subroutine runPrepeam() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
	Character(len=128) :: outPath, potPath, potPathDlpoly
!Options	
	!PRP1 
	
	outPath = trim(outputDirectory)
	potPath = trim(outPath)//"/"//trim(eamPreparedFile)
	potPathDlpoly = trim(potPath)//".dlpoly"
	
    If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Prepare potential file"
	End If
	print *,trim(outputDirectory)
	Call eamForceZBLCore(eamKey,eamData) 
	Call setPotentialDerivatives(eamKey,eamData) 
	Call outputEAMPrepData(eamKey, eamData, trim(potPath)) 
	Call outputEAMPrepDataDlpoly(eamKey, eamData, trim(potPathDlpoly))
	Call outputPrepareEAM()	        !Print output to file
	!Call outputPrepareDlpolyEAM()	
	!Call outputPrepareLammpsEAM()	
	  
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Saved to ",trim(eamPreparedFile)
	End If
	
	
  End Subroutine runPrepeam
  
 
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 

  
  
End Module prepeam  