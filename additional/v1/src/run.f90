Module run

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
  Use calc
  Use optimise
  
!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
  
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
	Integer(kind=StandardInteger) :: pointToVary
	Real(kind=DoubleReal) :: x,a,b,xV,differenceRSS,tempDouble
!Temp testing vars
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: distributionPoints
	Integer(kind=StandardInteger) :: distributionDataPoints
	Real(kind=DoubleReal) :: distX,distXInterval 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
!Write to output file
	If(mpiProcessID.eq.0)Then
!open output file	
	  open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
	  status="old",position="append",action="write")	
!write to output file
	  write(999,"(F8.4,A2,A33)") ProgramTime(),"  ",&
	  "Start calculations (runProcesses)"
!close the output file
      close(999) 	
    End If
	
!run subroutines	

!------------------------------------------------
! Calculation types
!------------------------------------------------
!ENER     1
!BMOD     3
!ECON     4
!EVAL     5
!OPTI     6
!PRP1    15
!PRP2    16
!PRP3    17
!TEST    99


!ENER    1     Calculate energy, forces, stresses of all input configurations
	If(calcRunType.eq.1)Then
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"Energy, force, stress calculation"
	  End If
	  eamDataSet = 1 	 
	  Call calcConfigEnergies()	!Calculate energies	  
	  Call calcOutput()	        !Print output to file
	  Call outputForces()
	End If
!BMOD    3     Calculate bulk modulus of all input configurations	
!ECON    4     Calculate elastic constants of all input configurations	
!EVAL    5     Evaluate energy, stress, forces
	If(calcRunType.eq.5)Then	!Evaluate	
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"Evaluate configurations (energy, stresses, forces, bulk properties)"
	  End If
	  eamDataSet = 1 	 	  
	  Call eamForceZBLCore(eamKey,eamData) 
	  Call setPotentialDerivatives(eamKey,eamData) 	  
	  Call calcEval()  
	  Call calcOutput()	        !Print output to file
	  Call outputForces()
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"RSS: ",trialResidualSquareSum
	  End If
	End If
!OPTI    6     Optimise the input potential	
	If(calcRunType.eq.6)Then	!Optimise
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"Optimise potential"
	  End If
	  Call storeEAMToFile(eamKey, eamData, trim(outputDirectory)//"/"//"vanilla.pot") 
	  Call eamForceZBLCore(eamKey,eamData) 
	  Call setPotentialDerivatives(eamKey,eamData) 
	  Call storeEAMToFile(eamKey, eamData, trim(outputDirectory)//"/"//"inputPrepared.pot") 
	  Call optimisePotential()	  
	End If
!OPTI    7     Optimise the input potential	full evaluation
	If(calcRunType.eq.7)Then	!Optimise
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"Optimise potential"
	  End If
	  Call storeEAMToFile(eamKey, eamData, trim(outputDirectory)//"/"//"vanilla.pot") 
	  Call eamForceZBLCore(eamKey,eamData) 
	  Call setPotentialDerivatives(eamKey,eamData) 
	  Call storeEAMToFile(eamKey, eamData, trim(outputDirectory)//"/"//"inputPrepared.pot") 
	  Call optimisePotential()	  
	End If
!EVAF    8     Evaluate bulk properties of all input configurations	
	If(calcRunType.eq.8)Then	!Evaluate	
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"Evaluate configurations (energy, stresses, forces, bulk properties)"
	  End If
	  eamDataSet = 1 	 	  
	  Call eamForceZBLCore(eamKey,eamData) 
	  Call setPotentialDerivatives(eamKey,eamData) 	  
	  Call calcEvalFull()  
	  Call calcOutput()	        !Print output to file
	  Call outputForces()
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		print *,ProgramTime(),"RSS: ",trialResidualSquareSum
	  End If
	End If
!TEST
	If(calcRunType.eq.99)Then    !Test
	  	  
	End If
	
	

	
  End Subroutine runProcessesAction 
  
  
  
  
!------------------------------------------------------------------------!
! Run Subroutines
!------------------------------------------------------------------------!  
 

End Module run