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
	Real(kind=DoubleReal) :: x,a,b,xV,differenceRSS
!Temp testing vars
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: distributionPoints
	Integer(kind=StandardInteger) :: distributionDataPoints
	Real(kind=DoubleReal) :: distX,distXInterval 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
!Write to output file
	If(mpiProcessID.eq.0)Then
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")		
!write to output file
	  write(999,"(F8.4,A2,A33)") ProgramTime(),"  ",&
	  "Start calculations (runProcesses)"
!close the output file
      close(999) 	
    End If
	
!run subroutines	
	
!choose type of calculation: eval energy forces bulkmodulus elasticconstants 
	If(calcRunType.eq.1)Then	!Calculates energy of configurations only (and forces/stresses) 
	  eamDataSet = 1 	 
	  Call calcConfigEnergies()	!Calculate energies	  
	  Call calcOutput()	        !Print output to file
	  Call calcOutputForces()
	End If
	
	If(calcRunType.eq.5)Then	!Evaluate	  
	  eamDataSet = 1 	  
	  Call calcEvalFull()  
	  Call calcOutput()	        !Print output to file
	  Call calcOutputForces()
	End If
	
	If(calcRunType.eq.6)Then	!Optimise
	  Call runOptimise()
	  
!Testing start
	  !Set data points	
	  !distributionDataPoints = 21
	  !Allocate(distributionPoints(1:distributionDataPoints,1:2))
	  !distX = 0.0D0
	  !distXInterval = 1.0D0/(distributionDataPoints-1)
!Make array
      !Do i=1,distributionDataPoints
	  !  distributionPoints(i,1) = distX
	  !  distributionPoints(i,2) = exp(-1*((distX-0.0D0)**2)/(2*0.28D0**2))
	 !   distX = distX + distXInterval
	 !   print *,i,distributionPoints(i,1),distributionPoints(i,2)
	 ! End Do	
	  !x = 0.8500000001D0
	 ! yArray = PointInterpolationArr(distributionPoints,x,5,1,21,"Y")
	  !print *,x,yArray(1),yArray(2)
!Testing end
	  
	  
	  
	End If
	
	
	
  End Subroutine runProcessesAction 
  
  
  
  
  
  
  
  Subroutine runOptimise()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, potPoints
	Integer(kind=StandardInteger) :: reducedPointCount, bestPotCounter
	Real(kind=DoubleReal) :: startingRSS, bestRSS
	Real(kind=DoubleReal) :: varyAmount, unperturbed, perturbAmount, x, y, dy, ddy
	Real(kind=DoubleReal) :: potX, potY, potYb, wobbliness
	Integer(kind=StandardInteger) :: potStart, potEnd
	Integer(kind=StandardInteger) :: potStartNew, potLengthNew
	Integer(kind=StandardInteger) :: varyType
	Integer(kind=StandardInteger), Dimension(:), Allocatable :: potType
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: interpPoints, makeJacobianMatrix
	Real(kind=DoubleReal), Dimension(:), Allocatable :: yArray
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: jacobianMatrix, jacobianMatrixT
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: hessianMatrix
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: potentialResponseResults
	Character(len=3) :: bestPotCounterChar
  
 !Set to use trial data points in calculations
    eamDataSet = 3
	varyType = 4 !1 = all points 2 = individual points
!Make reduced set of points
	Call makeReducedEAMSet()
!Store this set of points as the optimum
	Call storeReducedToOpt()
!Make set of points to use in calculations
	Call makeTrialEAMSet()
!Evaluate
	Call calcEval()
	startingRSS = trialResidualSquareSum
	bestRSS = trialResidualSquareSum
	
	If(varyType.eq.1)Then
!Vary all points at once
	  Do i=1,1000
!Load optimum reduced set of EAM pot data points	
	    Call loadOptToReduced()
!Vary reduced data points
        eamDataReduced = VaryPoints(eamDataReduced,0.01D0)
!Make set of points to use in calculations
	    Call makeTrialEAMSet()
!Evaluate
	    Call calcEval()
		print *,i,point,trialResidualSquareSum
		If(trialResidualSquareSum.gt.10)Then
		  Exit
		End If
		
!If better, store these points
        If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
	      bestRSS = trialResidualSquareSum
		  Call storeReducedToOpt()
	    End If	
	  End Do
	End If
	
	
	
    If(varyType.eq.2)Then
!Make pot type index
	  Allocate(potType(1:size(eamDataReduced,1)))
	  k = 0
	  Do i=1,size(eamKeyReduced,1)
	    Do j=eamKeyReduced(i,4),(eamKeyReduced(i,4)+eamKeyReduced(i,5)-1)
		  k = k + 1		  
		  potType(k) = eamKeyReduced(i,3)
		End Do
	  End Do
!Vary points one at a time
	  Do i=1,13000
!Load optimum reduced set of EAM pot data points	
	    Call loadOptToReduced()
!Point to vary
        point = mod(i,size(eamDataReduced,1))
		If(potType(point).eq.1)Then
		  varyAmount = 0.1D0
		End If
		If(potType(point).eq.2.or.potType(point).eq.4)Then
		  varyAmount = 0.1D0
		End If
		If(potType(point).eq.3.or.potType(point).eq.4)Then
		  varyAmount = 0.001D0
		End If
!Vary reduced data points
        eamDataReduced = VaryPoints(eamDataReduced,varyAmount,point)
!Make set of points to use in calculations
	    Call makeTrialEAMSet()
!Evaluate
	    Call calcEval()
		If(mpiProcessID.eq.0)Then
		  print *,i,point,varyAmount,trialResidualSquareSum
		End If
!If better, store these points
        If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
	      bestRSS = trialResidualSquareSum
		  Call storeReducedToOpt()
	    End If	
	  End Do
	End If
	
	
	
    If(varyType.eq.3)Then
!Make pot type index
	  Allocate(potType(1:size(eamDataReduced,1)))
	  k = 0
	  Do i=1,size(eamKeyReduced,1)
	    Do j=eamKeyReduced(i,4),(eamKeyReduced(i,4)+eamKeyReduced(i,5)-1)
		  k = k + 1		  
		  potType(k) = eamKeyReduced(i,3)
		End Do
	  End Do
!Calculate the unperturbed response function value
      Call loadOptToReduced()
	  Call calcEval()
      unperturbed = trialResidualSquareSum
      !Allocate(makeJacobianMatrix(1:size(eamDataReduced,1),1:4))  hessianMatrix
	  Allocate(interpPoints(1:4,1:2))
	  Allocate(jacobianMatrix(1:1,1:5))
	  Allocate(jacobianMatrixT(1:5,1:1))
!loop through each parameter and perturb each
	  !Do i=1,size(eamDataReduced,1)
	  k = 0
	  Do i=114,118
	    k = k+1
        If(potType(i).eq.1)Then
		  perturbAmount = 1.0D-04
		ElseIf(potType(i).eq.2.or.potType(i).eq.4)Then
		  perturbAmount = 1.0D-04
		ElseIf(potType(i).eq.3.or.potType(i).eq.5)Then
		  perturbAmount = 1.0D-04
		Else
		  perturbAmount = 1.0D-05
		End If      
!Perturb EAM potential
!First pertubation/first point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set
		eamDataReduced(i,2) = eamDataReduced(i,2)-perturbAmount
		Call calcEval()
		interpPoints(1,1) = eamDataReduced(i,2)	
		interpPoints(1,2) = 1.0D0*trialResidualSquareSum
!Store unperturbed amount/second point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set		
		interpPoints(2,1) = eamDataReduced(i,2)	
		interpPoints(2,2) = 1.0D0*unperturbed
!Second pertubation/third point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set
		eamDataReduced(i,2) = eamDataReduced(i,2)+perturbAmount
		Call calcEval()
		interpPoints(3,1) = eamDataReduced(i,2)	
		interpPoints(3,2) = 1.0D0*trialResidualSquareSum
!Second pertubation/third point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set
		eamDataReduced(i,2) = eamDataReduced(i,2)+2.0D0*perturbAmount
		Call calcEval()
		interpPoints(4,1) = eamDataReduced(i,2)	
		interpPoints(4,2) = 1.0D0*trialResidualSquareSum
!calculate first and second derivatives
		If(interpPoints(1,2).eq.interpPoints(2,2).and.&
		  interpPoints(1,2).eq.interpPoints(3,2).and.&
		  interpPoints(1,2).eq.interpPoints(4,2))Then
		  y = 0.0D0
		  dy = 0.0D0
		Else
		  yArray = PointInterpolationArr(interpPoints,interpPoints(2,1),4)
		  y = yArray(1)
		  dy = yArray(2)
		End If		
!make Jacobian matrix and transpose
	    jacobianMatrix(1,k) = dy
	    jacobianMatrixT(k,1) = dy
		hessianMatrix = matmul(jacobianMatrixT,jacobianMatrix)	
		
		print *,i
		print *,interpPoints(1,1),interpPoints(1,2)
		print *,interpPoints(2,1),interpPoints(2,2)
		print *,interpPoints(3,1),interpPoints(3,2)
		print *,interpPoints(4,1),interpPoints(4,2)
		print *,i,unperturbed,y,dy
        print *,""
		
	  End Do
	  
	  print *,""
	  Do i=1,5
	    print *,i,hessianMatrix(i,1),hessianMatrix(i,2),hessianMatrix(i,3),&
	    hessianMatrix(i,4),hessianMatrix(i,5)
	  End Do
	  print *,""
	  
	End If
	
	
	
	!Call makeTrialEAMSet()
	
    If(varyType.eq.4)Then	
!Plan
! Perturb each point
! Find which points affect potential
! Reduce eam set further
!
!Store the input eamReducedPoints
	  Call storeReducedToOpt()
! Choose set of eam data points
      eamDataSet = 3
!Set variables
      potPoints = size(eamDataReduced,1)
!Allocate array	
	  Allocate(potentialResponseResults(1:size(eamDataReduced,1),1:6))
!Make pot type index
	  Allocate(potType(1:potPoints))
	  k = 0
	  Do i=1,size(eamKeyReduced,1)
	    Do j=eamKeyReduced(i,4),(eamKeyReduced(i,4)+eamKeyReduced(i,5)-1)
		  k = k + 1		  
		  potType(k) = eamKeyReduced(i,3)
		End Do
	  End Do
!Calculate the unperturbed response function value
      Call loadOptToReduced()
	  Call makeTrialEAMSet()
	  Call calcEval()
	  wobbliness = eamTotalWobbliness(eamKeyTrial, eamDataTrial, 10000)
      unperturbed = trialResidualSquareSum
	  If(mpiProcessID.eq.0)Then
	    print *,"Unperturbed, starting RSS:",unperturbed
	    print *,"Total potential wobbliness:",wobbliness
	  End If
	  Allocate(interpPoints(1:4,1:2))
!loop through each parameter and perturb each
	  k = 0
	  Do i=1,size(eamDataReduced,1)
	    k = k+1
        If(potType(i).eq.1)Then
		  perturbAmount = 1.0D-04
		ElseIf(potType(i).eq.2.or.potType(i).eq.4)Then
		  perturbAmount = 1.0D-04
		ElseIf(potType(i).eq.3.or.potType(i).eq.5)Then
		  perturbAmount = 1.0D-04
		Else
		  perturbAmount = 1.0D-05
		End If      
!Perturb EAM potential
!First pertubation/first point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set
		eamDataReduced(i,2) = eamDataReduced(i,2)-perturbAmount
		Call makeTrialEAMSet()
		Call calcEval()
		interpPoints(1,1) = eamDataReduced(i,2)	
		interpPoints(1,2) = 1.0D0*trialResidualSquareSum
!Store unperturbed amount/second point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set		
		interpPoints(2,1) = eamDataReduced(i,2)	
		interpPoints(2,2) = 1.0D0*unperturbed
!Second pertubation/third point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set
		eamDataReduced(i,2) = eamDataReduced(i,2)+perturbAmount
		Call makeTrialEAMSet()
		Call calcEval()
		interpPoints(3,1) = eamDataReduced(i,2)	
		interpPoints(3,2) = 1.0D0*trialResidualSquareSum
!Second pertubation/third point
		Call loadOptToReduced()		!Load the optimised set of reduced points to the reduce eam set
		eamDataReduced(i,2) = eamDataReduced(i,2)+2.0D0*perturbAmount
		Call makeTrialEAMSet()
		Call calcEval()
		interpPoints(4,1) = eamDataReduced(i,2)	
		interpPoints(4,2) = 1.0D0*trialResidualSquareSum
!Pot points		
		potX = eamDataReduced(i,1)	
		potY = interpPoints(2,1)
!calculate first and second derivatives
		If(interpPoints(1,2).eq.interpPoints(2,2).and.&
		  interpPoints(1,2).eq.interpPoints(3,2).and.&
		  interpPoints(1,2).eq.interpPoints(4,2))Then
		  y = 0.0D0
		  dy = 0.0D0
		Else
		  yArray = PointInterpolationArr(interpPoints,interpPoints(2,1),4)
		  y = yArray(1)
		  dy = yArray(2)
		End If		
!Store potential response results
        potentialResponseResults(i,1) = i
        potentialResponseResults(i,2) = y
        potentialResponseResults(i,3) = dy
        potentialResponseResults(i,4) = 0.0D0
        potentialResponseResults(i,5) = potX
        potentialResponseResults(i,6) = potY		
	  End Do
!Store reduced set of potentials
	  k = 0
	  reducedPointCount = 0
	  Do i=1,size(eamKeyReduced,1)
	    potStart = eamKeyReduced(i,4)
		potEnd = (eamKeyReduced(i,4)+eamKeyReduced(i,5)-1)
		potLengthNew = 0
		potStartNew = reducedPointCount + 1
	    Do j=potStart,potEnd
		  k = k + 1		  
		  If(j.eq.potStart.or.j.eq.potEnd)Then
		    reducedPointCount = reducedPointCount + 1
			potLengthNew = potLengthNew + 1
			potentialResponseResults(k,4) = 1.0D0
		  Else
		    If(potentialResponseResults(k,3).ne.0)Then
		      reducedPointCount = reducedPointCount + 1
			  potLengthNew = potLengthNew + 1
			  potentialResponseResults(k,4) = 1.0D0
			End If
		  End If
		  If(mpiProcessID.eq.0.and.potentialResponseResults(k,4).eq.1.0D0)Then
		    !print *,i,k,potentialResponseResults(k,5),potentialResponseResults(k,6)
          End If	
		End Do
!Update eamKeyReduced
		eamKeyReduced(i,4) = potStartNew
		eamKeyReduced(i,5) = potLengthNew
		!If(mpiProcessID.eq.0)Then
        !  print *,i,potStartNew,potLengthNew   
        !End If		  
	  End Do
!Make new reduced set of points
	  Deallocate(eamDataReduced)
	  Allocate(eamDataReduced(1:reducedPointCount,1:2))
	  reducedPointCount = 0
      Do i=1,potPoints
	    If(potentialResponseResults(i,4).eq.1.0D0)Then
		  reducedPointCount = reducedPointCount + 1
	      eamDataReduced(reducedPointCount,1) = potentialResponseResults(i,5)
	      eamDataReduced(reducedPointCount,2) = potentialResponseResults(i,6)
		!If(mpiProcessID.eq.0)Then
		!  print *,i,eamDataReduced(reducedPointCount,1),eamDataReduced(reducedPointCount,2)
		!End If  
		End If
	  End Do
!Synch MPI processes
      Call synchMpiProcesses()
!Run evaluation
      Call makeTrialEAMSet()
	  Call calcEval()
	  Call calcOutput()	
	  Call storeReducedToOpt()
!Starting RSS
	  print *,"Starting RSS:",trialResidualSquareSum
!Find better points
!Make pot type index
      Deallocate(potType)
	  Allocate(potType(1:size(eamDataReduced,1)))
	  k = 0
	  Do i=1,size(eamKeyReduced,1)
	    Do j=eamKeyReduced(i,4),(eamKeyReduced(i,4)+eamKeyReduced(i,5)-1)
		  k = k + 1		  
		  potType(k) = eamKeyReduced(i,3)
		End Do
	  End Do
!Vary points one at a time
      bestPotCounter = 0
	  Do i=1,(30*size(eamDataReduced,1))
!Load optimum reduced set of EAM pot data points	
	    Call loadOptToReduced()
!Point to vary
        point = mod(i,size(eamDataReduced,1))
		If(potType(point).eq.1)Then
		  varyAmount = 0.05D0
		End If
		If(potType(point).eq.2.or.potType(point).eq.4)Then
		  varyAmount = 0.05D0
		End If
		If(potType(point).eq.3.or.potType(point).eq.4)Then
		  varyAmount = 0.001D0
		End If
!Vary reduced data points
        potYb = eamDataReduced(point,2)
        eamDataReduced = VaryPoints(eamDataReduced,varyAmount,point)
!Make set of points to use in calculations
	    Call makeTrialEAMSet()
!Get function wobbliness
        wobbliness = eamTotalWobbliness(eamKeyTrial, eamDataTrial, 10000)
!Evaluate
	    Call calcEval()
		If(mpiProcessID.eq.0)Then
		  potX = eamDataReduced(point,1)
		  potY = eamDataReduced(point,2)
		  !print *,i,point,potX,potYb,potY,varyAmount,trialResidualSquareSum
		  print *,i,point,potX,potYb,potY,varyAmount,trialResidualSquareSum,wobbliness
		End If
!If better, store these points
        If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
		  bestPotCounter = bestPotCounter + 1
	      bestRSS = trialResidualSquareSum
		  Call storeReducedToOpt()
		  Call storeEAMToFile(eamKeyReduced, eamDataReduced, "bestreduced.pot")
		  Call storeEAMToFile(eamKeyTrial, eamDataTrial, "besttrial.pot")
		  Call storeEAMToFile(eamKeyReduced, eamDataReduced,&
		  "bestreduced"//".pot",10000)
	    End If	
	  End Do
	  
	End If
	
	
	
	
	
!Run final test
    Call loadOptToReduced()
	Call makeTrialEAMSet()
	Call calcEval()
	
	Print *,"Best RSS:",bestRSS
	
  
  
  End Subroutine runOptimise 
  
  
  
  
  
  

End Module run