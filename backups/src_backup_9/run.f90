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
	
		
	If(calcRunType.eq.8)Then    !Eval Spline
	  eamDataSet = 3 	 
      Call makeReducedEAMSet()
	  Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)	  
	  Call calcEvalFull()  
	  Call calcOutput()	        !Print output to file
	  Call calcOutputForces()
	  If(mpiProcessID.eq.0)Then
	    print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	  End If 
	End If
	
	
	
	
	If(calcRunType.eq.99)Then    !Test
	  	  
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
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOriginal	!reduced set of points optimised
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOriginal          !reduced set of points optimised
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
!
!
!----------------------------------------------------------------------------
! Step 1 - run evaluation for input potential
!----------------------------------------------------------------------------
! Choose set of eam data points for eam calculations
    eamDataSet = 1	!eamData
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	  !Call storeEAMToFile(eamKey, eamData, "pot1_input.pot")
	End If 
	eamKeyOriginal = eamKey
	eamDataOriginal = eamData
!----------------------------------------------------------------------------
! Step 2 - make reduced set and evaluate
!----------------------------------------------------------------------------
    eamDataSet = 3	!eamDataTrial
!Store the input eamReducedPoints and output to files
	Call makeReducedEAMSet()
	Call storeReducedToOpt()
	Call makeExpandedEAMSet(eamKeyReduced,eamDataReduced,eamKeyTrial,eamDataTrial)
	!Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)
	
	
!Calculate RSS
    eamDataSet = 3	!eamDataTrial
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Reduced splined potential RSS:",trialResidualSquareSum  !Starting RSS
	  !Call storeEAMToFile(eamKeyReduced, eamDataReduced, "pot2_reduced.pot")
	  !Call storeEAMToFile(eamKeyTrial, eamDataTrial, "pot3_spline.pot")
	End If  
!----------------------------------------------------------------------------
! Step 3 - make further reduced set and evaluate
!----------------------------------------------------------------------------	
!Remove insignificant points
	!Call runRemoveInsignificantPoints(potentialResponseResults)
	!Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)
!Calculate RSS
	Call calcEval()
	If(mpiProcessID.eq.0)Then
	  print *,"Pre-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	  !Call storeEAMToFile(eamKeyReduced, eamDataReduced, "pot4_reduced.pot")	
	  !Call storeEAMToFile(eamKeyTrial, eamDataTrial, "pot5_spline.pot")
	End If 
!----------------------------------------------------------------------------
! Step 3 - Run SA optimisation
!----------------------------------------------------------------------------
    !eamKeyOptimised = eamKeyTrial
	!eamDataOptimised = eamDataTrial
    !Call storeEAMToFile(eamKeyOptimised, eamDataOptimised, "optimised.pot")
	

	!If(mpiProcessID.eq.0)Then
	!  print *,eamDataReduced(115,2)
	!End If  
		
	!Do i=101,130	
	!  eamDataReduced(i,2) = eamDataReduced(i,2)*3.2D1
	!End Do  

	!If(mpiProcessID.eq.0)Then
	!  print *,eamDataReduced(115,2)
	!End If  
	
	
	
	!Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)
	If(mpiProcessID.eq.0)Then
	  print *,"Pre-optimise (SA1) potential RSS2: ",trialResidualSquareSum  !Starting RSS
	  Call storeEAMToFile(eamKeyReduced, eamDataReduced, "optimisedR.pot")
	  Call storeEAMToFile(eamKeyTrial, eamDataTrial, "optimised.pot")
	End If 
	
	eamDataReduced(110,1) = eamDataReduced(110,1)/1.01D0
	Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)
	If(mpiProcessID.eq.0)Then
	  print *,"Pre-optimise (SA1) potential RSS3: ",trialResidualSquareSum  !Starting RSS
	  !Call storeEAMToFile(eamKeyReduced, eamDataReduced, "optimised.pot")
	End If 
	
	
	!Call runOptimisePointVary(eamKeyReduced,eamDataReduced,eamKeyOptimised,eamDataOptimised)
	!eamKeyOptimised = eamKeyTrial
	!eamDataOptimised = eamDataTrial
	!Call synchMpiProcesses()
	!Call storeToReduced(eamKeyOptimised,eamDataOptimised) 
!Set database set
	!eamDataSet = 1	!eamDataTrial
	!eamKey = eamKeyOptimised
	!eamData = eamDataOptimised
!Test optimised potential
	!Call calcEval()
	!Call calcOutput()	        !Print output to file
	!If(mpiProcessID.eq.0)Then
	!  print *,size(eamKeyOptimised,1),size(eamDataOptimised,1)
	!  print *,"Post-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	!  Call storeEAMToFile(eamKeyOptimised, eamDataOptimised, "optimised.pot")
	!End If 
	!Call storeEAMToFile(eamKeyOptimised, eamDataOptimised, "optimised.pot")		
	!Call synchMpiProcesses()
	

  
  End Subroutine runOptimise 
  
  
  
  Subroutine runOptimisePointVary(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
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
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyInput
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataInput
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyBest
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataBest
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOptimised
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOptimised
	Character(len=64) :: fileName
!----------------------------------------------------------------------------
! Step 1 - calculate best rss
!----------------------------------------------------------------------------	
!Synch processes	
	Call synchMpiProcesses()
! Choose set of eam data points for eam calculations
    eamDataSet = 3	!eamDataTrial
!Set best rss
    eamKeyOptimised = eamKeyInput
    eamDataOptimised = eamDataInput
    eamKeyReduced = eamKeyInput
	Call makeTrialEAMSet(1,eamKeyOptimised,eamDataOptimised)
	Call calcEval()
	bestRSS = trialResidualSquareSum
	If(mpiProcessID.eq.0)Then
	  print *,"Best (starting) RSS: ",bestRSS  !Best RSS
	End If 
!----------------------------------------------------------------------------
! Step 2 - set potType key and set "temperature" details
!----------------------------------------------------------------------------		
!Synch processes	
	Call synchMpiProcesses()
!Make pot type index
	Allocate(potType(1:size(eamDataOptimised,1)))
	k = 0
	Do i=1,size(eamKeyOptimised,1)
	  Do j=eamKeyOptimised(i,4),(eamKeyOptimised(i,4)+eamKeyOptimised(i,5)-1)
	    k = k + 1		  
		potType(k) = eamKeyOptimised(i,3)
	  End Do
	End Do 
!Vary points one at a time
    bestPotCounter = 0
	fileCounter = 0
	totVarPoints = (saCycles*size(eamDataOptimised,1))  
	tempCount = ceil(1.0D0*totVarPoints/saTIncs)
	tempInc = 1.0D0*(1.0D0*saTStart-1.0D0*saTEnd)/(1.0D0*saTIncs-1)
    tempCounter = 0     
    temp = saTStart
!----------------------------------------------------------------------------
! Step 3 - vary points
!----------------------------------------------------------------------------	
	Do i=1,totVarPoints
	!Do i=81,90
!Synch processes	
	  Call synchMpiProcesses()
!increment temperature counter
	  tempCounter = tempCounter + 1
!Force temp ge than 1	  
	  If(temp.lt.1.0D0)Then
	    temp = 1.0D0
	  End If
!Load optimum reduced set of EAM pot data points	
	  !Call loadOptToReduced()
!Point to vary  
      point = mod(i,size(eamDataOptimised,1))+1
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
      eamDataReduced = VaryPoints(eamDataOptimised,varyAmount,point)
!Make set of points to use in calculations
      Call clearTrialEAM()
	  Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)
!Evaluate
	  eamDataSet = 3	!eamDataTrial
	  Call calcEval()
!If better, store these points
      If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
	    eamDataOptimised = eamDataReduced	!Store working set as the optimised set
		bestPotCounter = bestPotCounter + 1
	    bestRSS = trialResidualSquareSum
	    If(mpiProcessID.eq.0)Then
		  Call storeEAMToFile(eamKeyOptimised, eamDataOptimised, "bestreduced.pot")
		  Call storeEAMToFile(eamKeyTrial, eamDataTrial, "besttrial.pot")
		  print "(F8.4,I8,A4,I8,I8,F16.8,A2)",&
		  temp,i," of ",totVarPoints,point,trialResidualSquareSum," *"
		End If	
		eamKeyBest = eamKeyTrial
		eamDataBest = eamDataTrial
	  Else
	    If(mpiProcessID.eq.0)Then
		  print "(F8.4,I8,A4,I8,I8,F16.8)",&
		  temp,i," of ",totVarPoints,point,trialResidualSquareSum
		End If	
	  End If	!
      If(tempCounter.ge.tempCount.and.point.eq.size(eamDataOptimised,1))Then
	    tempCounter = 0
		temp = temp - tempInc
	  End If
	End Do
!Synch MPI processes
	Call synchMpiProcesses()
	
!temp	
    
    Call clearTrialEAM()
	Call makeTrialEAMSet(1,eamKeyOptimised,eamDataOptimised)
	
	
	eamDataSet = 3	!eamDataTrial
	Call calcEval()
	If(mpiProcessID.eq.0)Then	
	  print *,"RSS3: ",trialResidualSquareSum  !Starting RSS
	End If
	
	eamKey = eamKeyTrial
	eamData = eamDataTrial
!Evaluate
	eamDataSet = 1	!eamDataTrial
	Call calcEval()
	If(mpiProcessID.eq.0)Then	
	  print *,"RSS1: ",trialResidualSquareSum  !Starting RSS
	End If
	
	If(mpiProcessID.eq.0)Then
	  Call storeEAMToFile(eamKey, eamData, "optimised.pot")
	End If
	
	eamKey = eamKeyBest
	eamData = eamDataBest
!Evaluate
	eamDataSet = 1	!eamDataTrial
	Call calcEval()
	If(mpiProcessID.eq.0)Then	
	  print *,"RSS1: ",trialResidualSquareSum  !Starting RSS
	End If
	
	If(mpiProcessID.eq.0)Then
	  Call storeEAMToFile(eamKey, eamData, "optimised.pot")
	End If
	
!Stroe the full potential
    !Deallocate(eamKeyOptimised)
	!Deallocate(eamDataOptimised)
	eamKeyOptimised = eamKeyBest
    eamDataOptimised = eamDataBest
	
  End Subroutine runOptimisePointVary
  
  
  
  
  
  Subroutine runRemoveInsignificantPoints(potentialResponseResults)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, potPoints
	Integer(kind=StandardInteger) :: reducedPointCount, bestPotCounter
	Real(kind=DoubleReal) :: startingRSS, bestRSS
	Real(kind=DoubleReal) :: varyAmount, unperturbed, perturbAmount, x, y, dy, ddy
	Integer(kind=StandardInteger) :: potStart, potEnd
	Integer(kind=StandardInteger) :: potStartNew, potLengthNew
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