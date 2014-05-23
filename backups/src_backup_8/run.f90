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
	  If(mpiProcessID.eq.0)Then
	    print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	  End If 
	End If
	
	If(calcRunType.eq.6)Then	!Optimise
	  Call runOptimise()	  
	End If
	
	If(calcRunType.eq.7)Then    !Eval Trial
	  eamDataSet = 3 	 
      Call makeReducedEAMSet()
	  Call makeTrialEAMSet()	  
	  Call calcEvalFull()  
	  Call calcOutput()	        !Print output to file
	  Call calcOutputForces()
	  If(mpiProcessID.eq.0)Then
	    print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	  End If 
	End If
	
	
	
	
	If(calcRunType.eq.99)Then    !Test
	  eamDataSet = 3	!eamData
	  eamKeyTrial = eamKey
	  eamDataTrial = eamData
	  Call splinePotential(eamKeyTrial,eamDataTrial)	
	  Call synchMPIProcesses()
	  Call calcEval()
	  Call synchMPIProcesses()
	  If(mpiProcessID.eq.0)Then
	    print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	  End If 
	  
!Store starting potential as optimum
	  eamKeyOpt = eamKey
	  eamDataOpt = eamData  
	  
	  !Call storeEAMToFile(eamKey, eamData, "pot0.pot")
	  Call splinePotential(eamKeyTrial,eamDataTrial,110,5.0D0)
	  Call storeEAMToFile(eamKeyTrial,eamDataTrial,"pot1.pot")
	  Call synchMPIProcesses()
	  Call calcEval()
	  If(mpiProcessID.eq.0)Then
	    print *,"Input Potential RSS3:",trialResidualSquareSum  !Starting RSS
	  End If
	  

	  
	  Do pointToVary=1,130	 
	    eamKeyTrial = eamKeyOpt
	    eamDataTrial = eamDataOpt	
        tempDouble = eamDataTrial(pointToVary,2)
	    Call splinePotential(eamKey,eamDataTrial,pointToVary,100.0D0)
	    Call synchMPIProcesses()
	    eamDataSet = 3	!eamData
		Call calcEval()
	    Call synchMPIProcesses()
	    If(mpiProcessID.eq.0)Then
	      print *,pointToVary,trialResidualSquareSum,eamDataTrial(pointToVary,1),&
		  tempDouble,eamDataTrial(pointToVary,2)  !Starting RSS
	    End If 
	  End Do
	  
	  
	  
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
	Real(kind=DoubleReal) :: potX, potY, potYb, wobbliness, curvatureValue, curveLengthValue
	Integer(kind=StandardInteger) :: potStart, potEnd
	Integer(kind=StandardInteger) :: potStartNew, potLengthNew
	Integer(kind=StandardInteger) :: varyType
	Integer(kind=StandardInteger), Dimension(:), Allocatable :: potType
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: interpPoints, makeJacobianMatrix
	Real(kind=DoubleReal), Dimension(:), Allocatable :: yArray, curvature, curveLength
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: jacobianMatrix, jacobianMatrixT
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: hessianMatrix
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: potentialResponseResults
	Character(len=3) :: bestPotCounterChar
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOptimised	!reduced set of points optimised
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOptimised          !reduced set of points optimised
!--------------
! Process
! - Evaluate input potential
! - Make reduced potential set of points, spline to make trial potential, evaluate
! - Remove insignificant potential data points
! - 
! - 
! - 
! - 
!--------------
! Choose set of eam data points for eam calculations
    eamDataSet = 1	!eamData
	Call calcEval()
	If(mpiProcessID.eq.0)Then
	  print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	End If 
	!Call printEAM(eamKey,eamData)
    eamDataSet = 3	!eamDataTrial
!Store the input eamReducedPoints and output to files
	Call makeReducedEAMSet()
	Call storeReducedToOpt()
	Call makeTrialEAMSet()
	!Call printEAM(eamKeyTrial,eamDataTrial)
!Calculate RSS
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Reduced splined potential RSS:",trialResidualSquareSum  !Starting RSS
	End If  
!Remove insignificant points
	!Call runRemoveInsignificantPoints(potentialResponseResults)
	!Call storeReducedToOpt()
	!Call makeTrialEAMSet()
	!Call storeEAMToFile(eamKeyTrial, eamDataTrial, "preopt_trial.pot")
	!Call storeEAMToFile(eamKeyReduced, eamDataReduced, "preopt_reduced.pot")
!Calculate RSS
	!Call calcEval()
	!If(mpiProcessID.eq.0)Then
	!  print *,"Pre-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	!End If  	
	!Do i=1,size(potentialResponseResults,1)
	!  If(mpiProcessID.eq.0)Then
	!    print *,"prr",potentialResponseResults(i,4),potentialResponseResults(i,5),&
	!    potentialResponseResults(i,6)
	!  End If	
	!End Do
!Run SA optimisation 1
	Call runOptimisePointVary(eamKeyOptimised,eamDataOptimised)
	!Call storeToReduced(eamKeyOptimised,eamDataOptimised)   
!Set database set
	eamDataSet = 1	!eamDataTrial
!Test optimised potential
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Post-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	End If 
	Call storeEAMToFile(eamKey, eamData, "optimised.pot")	
	
  
  End Subroutine runOptimise 
  
  
  
  Subroutine runOptimisePointVary(eamKeyOptimised,eamDataOptimised)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, potPoints, fileCounter, totVarPoints
	Integer(kind=StandardInteger) :: reducedPointCount, bestPotCounter, tempCount, tempCounter
	Real(kind=DoubleReal) :: startingRSS, bestRSS
	Real(kind=DoubleReal) :: varyAmount, unperturbed, perturbAmount, x, y, dy, ddy
	Real(kind=DoubleReal) :: potX, potY, potYb, wobbliness, curveLengthValue
	Real(kind=DoubleReal) :: temp, tempInc
	Integer(kind=StandardInteger) :: potStart, potEnd
	Integer(kind=StandardInteger) :: potStartNew, potLengthNew
	Integer(kind=StandardInteger) :: varyType
	Integer(kind=StandardInteger), Dimension(:), Allocatable :: potType
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: interpPoints, makeJacobianMatrix
	Real(kind=DoubleReal), Dimension(:), Allocatable :: yArray, curveLength
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: jacobianMatrix, jacobianMatrixT
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: hessianMatrix
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: potentialResponseResults
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOptimised
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOptimised
	Character(len=64) :: fileName
!Synch processes	
	Call synchMpiProcesses()
! Choose set of eam data points for eam calculations
    eamDataSet = 3	!eamDataTrial
!Set best rss
	Call loadOptToReduced()
	Call makeTrialEAMSet(1)
	Call calcEval()
	bestRSS = trialResidualSquareSum
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
    bestPotCounter = 0
	fileCounter = 0
	totVarPoints = (saCycles*size(eamDataReduced,1))  
	tempCount = ceil(1.0D0*totVarPoints/saTIncs)
	tempInc = 1.0D0*(1.0D0*saTStart-1.0D0*saTEnd)/(1.0D0*saTIncs-1)
    tempCounter = 0     
    temp = saTStart
	Do i=1,totVarPoints
!Synch processes	
	  Call synchMpiProcesses()
!increment temperature counter
	  tempCounter = tempCounter + 1
!Force temp ge than 1	  
	  If(temp.lt.1.0D0)Then
	    temp = 1.0D0
	  End If
!Load optimum reduced set of EAM pot data points	
	  Call loadOptToReduced()
!Point to vary  
      point = mod(i,size(eamDataReduced,1))+1
	  If(potType(point).eq.1)Then
	    varyAmount = 0.01D0*temp
	  End If
	  If(potType(point).eq.2.or.potType(point).eq.4)Then
	    varyAmount = 0.01D0*temp
	  End If
	  If(potType(point).eq.3.or.potType(point).eq.4)Then
	    varyAmount = 0.0001D0*temp
	  End If
!Vary reduced data points
      potYb = eamDataReduced(point,2)
      eamDataReduced = VaryPoints(eamDataReduced,varyAmount,point)
!Set database set
	  eamDataSet = 3	!eamDataTrial
!Make set of points to use in calculations
      Call clearTrialEAM()
	  Call makeTrialEAMSet(1)
	  If(mpiProcessID.eq.0.and.i.eq.380)Then
	    Call storeEAMToFile(eamKeyTrial, eamDataTrial, "opt380.pot")
      End If
!Evaluate
	  Call calcEval()
	  If(mpiProcessID.eq.0)Then
		potX = eamDataReduced(point,1)
		potY = eamDataReduced(point,2)
		!print *,i,point,potX,potYb,potY,varyAmount,trialResidualSquareSum
		If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
		  !print *,temp,tempCounter,tempCount
		  !print *,i,totVarPoints,point,potX,potYb,potY,varyAmount,trialResidualSquareSum,"*"
		  print "(F8.4,I8,A4,I8,I8,F16.8,A2)",&
		  temp,i," of ",totVarPoints,point,trialResidualSquareSum," *"
		Else
		  !print *,temp,tempCounter,tempCount
          !print *,i,totVarPoints,point,potX,potYb,potY,varyAmount,trialResidualSquareSum
		  print "(F8.4,I8,A4,I8,I8,F16.8)",&
		  temp,i," of ",totVarPoints,point,trialResidualSquareSum
        End If		  
	  End If
!If better, store these points
      If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
		bestPotCounter = bestPotCounter + 1
	    bestRSS = trialResidualSquareSum
		Call storeReducedToOpt()
	    If(mpiProcessID.eq.0)Then
		  Call storeEAMToFile(eamKeyReduced, eamDataReduced, "bestreduced.pot")
		  Call storeEAMToFile(eamKeyTrial, eamDataTrial, "besttrial.pot")
		End If		
	  End If	!
      If(tempCounter.ge.tempCount.and.point.eq.size(eamDataReduced,1))Then
	    tempCounter = 0
		temp = temp - tempInc
	  End If
	End Do
!Synch MPI processes
	Call synchMpiProcesses()
	If(mpiProcessID.eq.0)Then
	  print *,"End opt"
	  Call storeEAMToFile(eamKeyOpt, eamDataOpt, "opt.pot")
    End If
!store optimised reduced set of data points
    If(Allocated(eamKeyOptimised))Then
	  Deallocate(eamKeyOptimised)
	End If
    If(Allocated(eamDataOptimised))Then
	  Deallocate(eamDataOptimised)
	End If
	eamKeyOptimised = eamKeyOpt
	eamDataOptimised = eamDataOpt
	Call loadOptToReduced() 
	Call makeTrialEAMSet(1)
!Also input the main input eam potential to new optimised
    !If(Allocated(eamKey))Then
	!  Deallocate(eamKey)
	!End If
    !If(Allocated(eamData))Then
	!  Deallocate(eamData)
	!End If
!Make full "trial" set of points from opt
    !eamKeyReduced = eamKeyOpt
    !eamDataReduced = eamDataOpt
    !Call makeTrialEAMSet(1)
	!eamKey = eamKeyTrial
	!eamData = eamDataTrial
!Synch MPI processes
	Call synchMpiProcesses()
!Set database set
	eamDataSet = 3	!eamDataTrial	
!Test optimised potential
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Optimised potential RSS: ",trialResidualSquareSum  !Starting RSS
	End If
	eamKey = eamKeyTrial
	eamData = eamDataTrial
	eamDataSet = 3	!eamDataTrial
    Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Optimised potential RSS (3): ",trialResidualSquareSum  !Starting RSS
	End If	
	eamDataSet = 1	!eamDataTrial
    Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Optimised potential RSS (1): ",trialResidualSquareSum  !Starting RSS
	End If
	
  End Subroutine runOptimisePointVary
  
  
  
  
  
  Subroutine runRemoveInsignificantPoints(potentialResponseResults)
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
!Set variables
    potPoints = size(eamDataReduced,1)
!Print out	
	If(mpiProcessID.eq.0)Then
	  print *,"Remove insignificant points in eam potentials"
	End If
!Calculate response of perturbing each point in the reduced eam set of points
    Call runPerturbationResponse(potentialResponseResults)	
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
	  End Do
!Update eamKeyReduced
	  eamKeyReduced(i,4) = potStartNew
	  eamKeyReduced(i,5) = potLengthNew	  
	  If(mpiProcessID.eq.0)Then
	    print *,i,potStartNew,potLengthNew
	  End If
	End Do
!Make new reduced set of points
	Deallocate(eamDataReduced)
	Allocate(eamDataReduced(1:reducedPointCount,1:2))
	reducedPointCount = 0
    Do i=1,potPoints	  
	  !If(mpiProcessID.eq.0)Then
	    !print *,i,potentialResponseResults(i,4),potentialResponseResults(i,5),&
	    !potentialResponseResults(i,6)
      !End If	 
	  If(potentialResponseResults(i,4).eq.1.0D0)Then
	    reducedPointCount = reducedPointCount + 1
	      eamDataReduced(reducedPointCount,1) = potentialResponseResults(i,5)
	      eamDataReduced(reducedPointCount,2) = potentialResponseResults(i,6) 
	  End If
	End Do
!Synch MPI processes
    Call synchMpiProcesses()  
  End Subroutine runRemoveInsignificantPoints 
  
  
  
  
  
  
  Subroutine runPerturbationResponse(potentialResponseResults)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, potPoints
	Integer(kind=StandardInteger) :: reducedPointCount, bestPotCounter
	Real(kind=DoubleReal) :: startingRSS, bestRSS
	Real(kind=DoubleReal) :: varyAmount, unperturbed, perturbAmount, x, y, yb, dyb, dy, ddy, dr
	Real(kind=DoubleReal) :: potX, potY, potYb, wobbliness
	Integer(kind=StandardInteger) :: potStart, potEnd
	Integer(kind=StandardInteger) :: potStartNew, potLengthNew
	Integer(kind=StandardInteger) :: varyType
	Integer(kind=StandardInteger), Dimension(:), Allocatable :: potType
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: interpPoints, makeJacobianMatrix
	Real(kind=DoubleReal), Dimension(:), Allocatable :: yArray, yArrayB
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: potentialResponseResults
! Choose set of eam data points for eam calculations
    eamDataSet = 3	!eamDataTrial
!Set variables
    potPoints = size(eamDataReduced,1)
!Allocate array	
    If(Allocated(potentialResponseResults))Then
	  Deallocate(potentialResponseResults)
	End If
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
    unperturbed = trialResidualSquareSum
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
		x = 1.0D0 * interpPoints(2,1)
		y = 0.0D0
		dy = 0.0D0
		ddy = 0.0D0
		dr = 0.0D0		
		yb = 0.0D0
		dyb = 0.0D0
	  Else   
	    x = 1.0D0 * interpPoints(2,1)
	    interpPoints = 1.0D0*interpPoints
		yArrayB = PointInterpolationArr(interpPoints,x,4)	
		yb = yArrayB(1)
		dyb = yArrayB(2)
		yArray = MatrixInterpolation(interpPoints,x)
		y = yArray(2)
		dy = yArray(3)
		ddy = yArray(4)
		dr = 1.0D0*abs(1.0D0*interpPoints(2,2)-1.0D0*interpPoints(1,2))+&
		     1.0D0*abs(1.0D0*interpPoints(3,2)-1.0D0*interpPoints(2,2))+&
		     1.0D0*abs(1.0D0*interpPoints(4,2)-1.0D0*interpPoints(3,2))
	  End If		
!Store potential response results
      potentialResponseResults(i,1) = i
      potentialResponseResults(i,2) = y
      potentialResponseResults(i,3) = dy
      potentialResponseResults(i,4) = ddy
      potentialResponseResults(i,5) = potX
      potentialResponseResults(i,6) = potY		
!Print out      
	  If(mpiProcessID.eq.0)Then
	    print *,i,potX,potY,yb,dyb,y,dy,ddy,dr
	  End If
	End Do
  End Subroutine runPerturbationResponse 
  
  
  

End Module run