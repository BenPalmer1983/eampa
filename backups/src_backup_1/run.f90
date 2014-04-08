Module run

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use prep
  Use calc
  
!force declaration of all variables
  Implicit None
  
!declare global variables 

!Privacy of functions/subroutines/variables
  Private
  
!Variables
  
!Subroutines
  Public :: runProcesses		  
  
!Functions
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!------------------------------------------------------------------------!
! runProcesses
!------------------------------------------------------------------------!
  
  Subroutine runProcesses()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k
!run
    Call runProcessesAction()
	
  End Subroutine runProcesses    
  
  
  
  
  
  
  Subroutine runProcessesAction()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")		
!write to output file
	write(999,"(F8.4,A2,A33)") ProgramTime(),"  ",&
	"Start calculations (runProcesses)"


!run subroutines	
	
!choose type of calculation: eval energy forces bulkmodulus elasticconstants 
	If(calcRunType.eq.1)Then	!Energy
	  Call calcConfigEnergies(.true.,.true.)
	  Call calcDifference(.true.)	  
	End If
	
	If(calcRunType.eq.5)Then	!Evaluate	  
	  !Call makeReducedEAMSet()
	  !Call makeTrialEAMSet()
	  !Call calcEval()  
	End If
	
	If(calcRunType.eq.6)Then	!Optimise
	  !Call calcConfigEnergies(.true.,.true.)
	  !Call calcDifference(.true.)	  
	End If
	
	
	
!close the output file
    close(999) 	
	
  End Subroutine runProcessesAction 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

End Module run