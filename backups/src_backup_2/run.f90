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
	Real(kind=DoubleReal) :: x,a,b,xV
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
	  eamDataSet = 3
	  Call makeReducedEAMSet()
	  Call makeTrialEAMSet()
	  Call calcEval()  	  
	  Call makeReducedEAMSet()
	  Call makeTrialEAMSet()
	  Call calcEval()  
	End If
	
	If(calcRunType.eq.6)Then	!Optimise
	  Call runOptimise()
	End If
	
	
!close the output file
    close(999) 	
	
  End Subroutine runProcessesAction 
  
  
  
  
  
  
  
  Subroutine runOptimise()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k
  
 !Set to use trial data points in calculations
    eamDataSet = 3
!Make reduced set of points
	Call makeReducedEAMSet()
!Store this set of points as the optimum
	Call storeReducedToOpt()
!Make set of points to use in calculations
	Call makeTrialEAMSet()
!Evaluate
	Call calcEval()
	print *,configurationEnergy(1),trialResidualSquareSum
	  
	Do i=1,10
!Load optimum reduced set of EAM pot data points	
	  Call loadOptToReduced()
!Vary reduced data points
      eamDataReduced = VaryPoints(eamDataReduced,0.1D0,0.1D0,1)
!Make set of points to use in calculations
	  Call makeTrialEAMSet()
!Evaluate
	  Call calcEval()	
	  !print *,configurationEnergy(1),trialResidualSquareSum
	End Do
	
	!Call calcConfigEnergies(.true.,.true.)
	!Call calcDifference(.true.)	  
  
  
  
  End Subroutine runOptimise 
  
  
  
  
  
  

End Module run