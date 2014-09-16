Module run

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use prep
  Use calc
  
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
      !Call makeReducedEAMSet()
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
      !Call makeReducedEAMSet()
	  !Call makeTrialEAMSet(1,eamKeyReduced,eamDataReduced)	  
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
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyR	!reduced set of points optimised
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataR         !reduced set of points optimised
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyW	!reduced set of points optimised
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataW         !reduced set of points optimised
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
	Call reReadEamPot(potentialFilePathTemp)	
    eamDataSet = 1	!eamData
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	End If 
!----------------------------------------------------------------------------
! Step 2 - make reduced set and evaluate
!----------------------------------------------------------------------------
!Store the input eamReducedPoints and output to files
	Call makeReducedEAMSet(eamKey,eamData,eamKeyR,eamDataR)		!Make reduced set
	Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)	!Make expanded set
	Call storeEAMToFile(eamKey, eamData, "temp.pot")			!Store to file	
	Call reReadEamPot("temp.pot")								!Re-read potential
	eamDataSet = 1	!eamDataTrial
	Call calcEval()
	Call calcOutput()	        !Print output to file
	If(mpiProcessID.eq.0)Then
	  print *,"Reduced splined potential RSS:",trialResidualSquareSum  !Starting RSS
	End If  
!----------------------------------------------------------------------------
! Step 3 - make further reduced set and evaluate
!----------------------------------------------------------------------------	
!Remove insignificant points
	!Call runRemoveInsignificantPoints(eamKeyR,eamDataR,potentialResponseResults)	
	Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)	!Make expanded set
	Call storeEAMToFile(eamKey, eamData, "temp.pot")			!Store to file	
	Call reReadEamPot("temp.pot")								!Re-read potential
!Calculate RSS
    eamDataSet = 1	!eamDataTrial
	Call calcEval()
	If(mpiProcessID.eq.0)Then
	  print *,"Pre-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	End If 
!----------------------------------------------------------------------------
! Step 3 - Run SA optimisation
!----------------------------------------------------------------------------
    Call makeReducedEAMSet(eamKey,eamData,eamKeyR,eamDataR)		!Make reduced set
    Call runOptimisePointVary(eamKeyR,eamDataR,eamKeyOptimised,eamDataOptimised)
	
	
	!Call makeExpandedEAMSet(eamKeyOptimised,eamDataOptimised,eamKey,eamData)	!Make expanded set
	!Call storeEAMToFile(eamKey, eamData, "opt.pot")			!Store to file	
	!Call reReadEamPot("opt.pot")								!Re-read potential
	
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
!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
  Subroutine runRemoveInsignificantPoints(eamKeyW,eamDataW,potentialResponseResults)
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
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyW	!eam key
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataW          !eam data
!Set variables
    potPoints = size(eamDataW,1)
!Print out	
	If(mpiProcessID.eq.0)Then
	  print *,"Remove insignificant points in eam potentials"
	End If
!Calculate response of perturbing each point in the reduced eam set of points
    Call runPerturbationResponse(eamKeyW,eamDataW,potentialResponseResults)	
!Store reduced set of potentials
	k = 0
	reducedPointCount = 0
	Do i=1,size(eamKeyW,1)
	  potStart = eamKeyW(i,4)
	  potEnd = (eamKeyW(i,4)+eamKeyW(i,5)-1)
	  potLengthNew = 0
	  potStartNew = reducedPointCount + 1
	  Do j=potStart,potEnd
	    k = k + 1		  
		If(j.eq.potStart.or.j.eq.potEnd)Then
		  reducedPointCount = reducedPointCount + 1
		  potLengthNew = potLengthNew + 1
		  potentialResponseResults(k,4) = 1.0D0
		Else
		  If(potentialResponseResults(k,3).ne.0.0D0)Then
		    reducedPointCount = reducedPointCount + 1
		    potLengthNew = potLengthNew + 1
			potentialResponseResults(k,4) = 1.0D0
		  End If
		End If	
	  End Do
!Update eamKeyW
	  eamKeyW(i,4) = potStartNew
	  eamKeyW(i,5) = potLengthNew	  
	  !If(mpiProcessID.eq.0)Then
	  !  print *,i,potStartNew,potLengthNew
	  !End If
	End Do
!Make new reduced set of points
	Deallocate(eamDataW)
	Allocate(eamDataW(1:reducedPointCount,1:2))
	reducedPointCount = 0
    Do i=1,potPoints	 	 
	  If(potentialResponseResults(i,4).eq.1.0D0)Then
	    reducedPointCount = reducedPointCount + 1
	      eamDataW(reducedPointCount,1) = potentialResponseResults(i,5)
	      eamDataW(reducedPointCount,2) = potentialResponseResults(i,6) 
	  End If
	End Do
!Synch MPI processes
    Call synchMpiProcesses()  
  End Subroutine runRemoveInsignificantPoints   
!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
  Subroutine runPerturbationResponse(eamKeyW,eamDataW,potentialResponseResults)
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
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyW	!eam key
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataW          !eam data
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyP	!eam key
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataP          !eam data
!Set variables
    potPoints = size(eamDataW,1)
!Allocate array	
    If(Allocated(potentialResponseResults))Then
	  Deallocate(potentialResponseResults)
	End If
	Allocate(potentialResponseResults(1:size(eamDataW,1),1:6))
!Make pot type index
	Allocate(potType(1:potPoints))
	k = 0
	Do i=1,size(eamKeyW,1)
	  Do j=eamKeyW(i,4),(eamKeyW(i,4)+eamKeyW(i,5)-1)
	    k = k + 1		  
	    potType(k) = eamKeyW(i,3)
	  End Do
	End Do
!Calculate the unperturbed response function value
    Call makeExpandedEAMSet(eamKeyW,eamDataW,eamKeyTrial,eamDataTrial)
    eamDataSet = 3	!eamDataTrial
	Call calcEval()
    unperturbed = trialResidualSquareSum
    Allocate(interpPoints(1:4,1:2))
!loop through each parameter and perturb each
	k = 0
	Do i=1,size(eamDataW,1)
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
	  eamKeyP = eamKeyW
	  eamDataP = eamDataW	  
	  eamDataP(i,2) = eamDataP(i,2)-perturbAmount
      Call makeExpandedEAMSet(eamKeyP,eamDataP,eamKeyTrial,eamDataTrial)
      eamDataSet = 3	!eamDataTrial
	  Call calcEval()
	  interpPoints(1,1) = eamDataP(i,2)	
	  interpPoints(1,2) = 1.0D0*trialResidualSquareSum
!Store unperturbed amount/second point
	  interpPoints(2,1) = eamDataW(i,2)	
	  interpPoints(2,2) = 1.0D0*unperturbed
!Second pertubation/third point
	  eamKeyP = eamKeyW
	  eamDataP = eamDataW	  
	  eamDataP(i,2) = eamDataP(i,2)+perturbAmount
      Call makeExpandedEAMSet(eamKeyP,eamDataP,eamKeyTrial,eamDataTrial)
	  Call calcEval()
	  interpPoints(3,1) = eamDataP(i,2)	
	  interpPoints(3,2) = 1.0D0*trialResidualSquareSum
!Third pertubation/fourth point
	  eamKeyP = eamKeyW
	  eamDataP = eamDataW	  
	  eamDataP(i,2) = eamDataP(i,2)+2.0D0*perturbAmount
      Call makeExpandedEAMSet(eamKeyP,eamDataP,eamKeyTrial,eamDataTrial)
	  Call calcEval()
	  interpPoints(4,1) = eamDataP(i,2)	
	  interpPoints(4,2) = 1.0D0*trialResidualSquareSum
!Pot points		
	  potX = eamDataW(i,1)	
	  potY = interpPoints(2,1)
!calculate first and second derivatives
	  If(interpPoints(1,2).eq.interpPoints(2,2).and.&
	    interpPoints(1,2).eq.interpPoints(3,2).and.&
	    interpPoints(1,2).eq.interpPoints(4,2))Then
		x = 1.0D0 * interpPoints(2,1)
		y = 0.0D0
		dy = 0.0D0
		ddy = 0.0D0
	  Else   
	    x = 1.0D0 * interpPoints(2,1)
	    interpPoints = 1.0D0*interpPoints
		yArray = MatrixInterpolation(interpPoints,x)
	    y = yArray(2)
		dy = yArray(3)
		ddy = yArray(4)
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
	    print *,i,potX,potY,y,dy,ddy
	  End If
	End Do
  End Subroutine runPerturbationResponse 
!---------------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------------
  Subroutine runOptimisePointVary(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, fileCount
	Real(kind=DoubleReal) :: bestRSS
	Real(kind=DoubleReal) :: x, xV
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: &
	eamKeyInput,eamKeyOptimised,eamKeyR
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: &
	eamDataInput,eamDataOptimised,eamDataR 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: randomSet
	Character(len=64) :: fileName
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,processCount,bufferSize
    Integer(kind=StandardInteger) :: processID
	Real(kind=DoubleReal) :: send, receive, buffer	
!MPI
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!----------------------------------------------------------------------------
! Step 1 - store optimum potential
!----------------------------------------------------------------------------		
	Call makeExpandedEAMSet(eamKeyInput,eamDataInput,eamKey,eamData)
	Call storeEAMToFile(eamKey, eamData, "tempOpt.pot")			!Store to file
!----------------------------------------------------------------------------
! Step 2 - calculate best rss
!----------------------------------------------------------------------------	
!Set best rss   
	eamDataSet = 1	
	Call calcEval()
	bestRSS = trialResidualSquareSum
	If(mpiProcessID.eq.0)Then
	  print *,"Best (starting) RSS: ",trialResidualSquareSum  !Best RSS
	End If 
!----------------------------------------------------------------------------
! Step 3 - vary points
!----------------------------------------------------------------------------
    fileCount = 1000
	Do i=40,50	
	  Call synchMpiProcesses()	!Synchronise processes
!Load potential from file
	  Call reReadEamPot("tempOpt.pot")								!Re-read potential
	  Call makeReducedEAMSet(eamKey,eamData,eamKeyR,eamDataR)		!Make reduced set
!Vary point
	  point = mod(i,size(eamDataR,1))+1
	  x = eamDataR(point,2)
	  !xV = VaryPointRand(x,40.0D0,0.1D0)
	  !randNumber = randomSet(i+1)
	  !xV = x*(1.0D0+(randNumber/100.0D0))
	  !xV = 1.01D0 * x
	  If(processID.eq.0)Then
	    xV = VaryPointRand(x,40.0D0,0.3D0) 
		xV = 1.01D0 * x
		send = xV
        Do j=1,(processCount-1)
          tag = 3140 + j	
          Call MPI_send(send,1,MPI_double_precision,j,tag,MPI_comm_world,status,error)
        End Do
      Else
        tag = 3140 + processID
	    Call MPI_recv(receive,1,MPI_double_precision,0,tag,MPI_comm_world,status,status,error)  
		xV = receive
      End If  
	  eamDataR(point,2) = xV
!expand, store and re-read
	  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)  !Make expanded set
!Run calculation
	  eamDataSet = 1	
	  Call calcEval()
	  Call calcOutput()	        !Print output to file
!If better, store
	  If(trialResidualSquareSum.lt.bestRSS)Then
	    bestRSS = trialResidualSquareSum		
	    Call storeEAMToFile(eamKey, eamData, "tempOpt.pot")			!Store to file
		eamKeyOptimised = eamKey
		eamDataOptimised = eamData
	    If(mpiProcessID.eq.0)Then
	      print *,"RSS: ",i,point,trialResidualSquareSum,x,xV,"*"
	    End If
	  Else
	    If(mpiProcessID.eq.0)Then
	      print *,"RSS: ",i,point,trialResidualSquareSum,x,xV
	    End If
	  End If	  	  
	End Do
	
	
	
	eamKey = eamKeyOptimised
	eamData = eamDataOptimised
	eamDataSet = 1	
	Call calcEval()
	
	If(mpiProcessID.eq.0)Then
	  Call storeEAMToFile(eamKey, eamData, "optimised.pot")
	  print *,"Optimised RSS: ",trialResidualSquareSum  !Best RSS
	End If
	
	
  End Subroutine runOptimisePointVary
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    Subroutine runOptimisePointVaryB(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
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
	
  End Subroutine runOptimisePointVaryB

  
    Subroutine runOptimisePointVaryC(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
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
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyR
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataR
	Character(len=64) :: fileName
!----------------------------------------------------------------------------
! Step 1 - calculate best rss
!----------------------------------------------------------------------------	
!Synch processes	
	Call synchMpiProcesses()
!Store optimal potential
    eamKeyOptimised = eamKeyInput
    eamDataOptimised = eamDataInput
	eamKeyR = eamKeyInput
	eamDataR = eamDataInput
!Set best rss
    Call makeExpandedEAMSet(eamKeyOptimised,eamDataOptimised,eamKeyTrial,eamDataTrial)
	eamDataSet = 3	
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
      eamDataR = VaryPoints(eamDataOptimised,varyAmount,point)
!Make set of points to use in calculations
	  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTrial,eamDataTrial)
!Evaluate
	  eamDataSet = 3	!eamDataTrial
	  Call calcEval()
!If better, store these points
      If((1.0D0*trialResidualSquareSum).lt.(1.0D0*bestRSS))Then
	    eamDataOptimised = eamDataR	!Store working set as the optimised set
		bestPotCounter = bestPotCounter + 1
	    bestRSS = trialResidualSquareSum
	    If(mpiProcessID.eq.0)Then
		  Call storeEAMToFile(eamKeyOptimised, eamDataOptimised, "bestreduced.pot")
		  Call storeEAMToFile(eamKeyTrial, eamDataTrial, "besttrial.pot")
		  print "(F8.4,I8,A4,I8,I8,F16.8,A2)",&
		  temp,i," of ",totVarPoints,point,trialResidualSquareSum," *"
		End If	
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
!Store optimised
	Call makeExpandedEAMSet(eamKeyOptimised,eamDataOptimised,eamKeyTrial,eamDataTrial)	
	eamDataSet = 3	
	Call calcEval()
	If(mpiProcessID.eq.0)Then
	  Call storeEAMToFile(eamKeyTrial, eamDataTrial, "optimised.pot")
	  print *,"Optimised RSS: ",trialResidualSquareSum  !Best RSS
	End If
	
	
  End Subroutine runOptimisePointVaryC
  
  
  

End Module run