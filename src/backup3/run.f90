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
	  eamDataSet = 1 	 
	  Call calcConfigEnergies()	!Calculate energies	  
	  Call calcOutput()	        !Print output to file
	  Call calcOutputForces()
	End If
!BMOD    3     Calculate bulk modulus of all input configurations	
!ECON    4     Calculate elastic constants of all input configurations	
!EVAL    5     Evaluate bulk properties of all input configurations	
	If(calcRunType.eq.5)Then	!Evaluate	
	  eamDataSet = 1 	  
	  Call calcEvalFull()  
	  Call calcOutput()	        !Print output to file
	  Call calcOutputForces()
	End If
!OPTI    6     Optimise the input potential	
	If(calcRunType.eq.6)Then	!Optimise
	  Call runOptimise()	  
	End If
!PRP1   15     Prepare a potential file	
	If(calcRunType.eq.15)Then	!output potential in formatted file
	  Call setPotentialDerivatives(eamKey,eamData) 
	  Call storeEAMToFile(eamKey, eamData, trim(eamPreparedFile)) 
	End If

	
	
	
	!
	
	
	If(calcRunType.eq.99)Then    !Test
	  	  
	End If
	
	
	
	
	!EVAL    6     Optimise the input potential	
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
	
	
  End Subroutine runProcessesAction 
  
  
  
  
!------------------------------------------------------------------------!
! Run Subroutines
!------------------------------------------------------------------------!  
  
!------------------------------------------------------------------------!
! Run Optimise Subroutines
!------------------------------------------------------------------------!  
  
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
	Integer(kind=StandardInteger) :: i, j, k, point, fileCount, startI, endI
	Real(kind=DoubleReal) :: bestRSS
	Real(kind=DoubleReal) :: x, xV
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: &
	eamKeyInput,eamKeyOptimised,eamKeyR,eamKeyBest,eamKeyTemp
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: &
	eamDataInput,eamDataOptimised,eamDataR,eamDataBest,eamDataTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: randomSet
	Character(len=64) :: fileName
	Logical :: writeFile
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,processCount,bufferSize
    Integer(kind=StandardInteger) :: processID
	Real(kind=DoubleReal) :: send, receive, buffer	
	Call calcEval()
!----------------------------------------------------------------------------
! Step 1 - store starting potential
!----------------------------------------------------------------------------		
    eamKeyOptimised = eamKeyInput	
    eamDataOptimised = eamDataInput
!----------------------------------------------------------------------------
! Step 2 - vary points
!----------------------------------------------------------------------------		
    startI = 0
	endI = 1000
	Do i=startI,endI
!Synchronise processes		
	  Call synchMpiProcesses()	
!Load optimum reduced points
      eamKeyR = eamKeyOptimised	
      eamDataR = eamDataOptimised	 
!Set best rss
	  If(i.eq.startI)Then
	    !make expanded set of points
	    Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTemp,eamDataTemp)
!Store temp potential	  
	    Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"tempfull.pot")	  
!Synchronise processes		  
	    Call synchMpiProcesses()	
!Load potential from file on to all processes	  
	    Call reReadEamPot("tempfull.pot")	
!Calculate rss
	    eamDataSet = 1		
        Call calcEval()
		bestRSS = trialResidualSquareSum
		If(mpiProcessID.eq.0)Then
	      print *,"Start RSS: ",trialResidualSquareSum
		End If  
	  End If
!Vary the point
	  point = mod(i,size(eamDataR,1))+1
	  x = eamDataR(point,2)
!Vary on master, share to workers
	  If(mpiProcessID.eq.0)Then
	    xV = VaryPointRand(x,40.0D0,0.05D0) 
		send = xV
        Do j=1,(mpiProcessCount-1)
          tag = 3140 + j	
          Call MPI_send(send,1,MPI_double_precision,j,tag,MPI_comm_world,status,error)
        End Do
      Else
        tag = 3140 + mpiProcessID
	    Call MPI_recv(receive,1,MPI_double_precision,0,tag,MPI_comm_world,status,status,error)  
		xV = receive
      End If  
	  eamDataR(point,2) = xV
!make expanded set of points
	  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTemp,eamDataTemp)
!Store temp potential	  
	  Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"tempfull.pot")	  
!Synchronise processes		  
	  Call synchMpiProcesses()	
!Load potential from file on to all processes	  
	  Call reReadEamPot("tempfull.pot")	
!Calculate rss
	  eamDataSet = 1		
      Call calcEval()	
!Check if optimum
	  If(trialResidualSquareSum.lt.bestRSS)Then	
	    bestRSS = trialResidualSquareSum
		eamKeyOptimised = eamKeyR
		eamDataOptimised = eamDataR
		Call storeEAMToFile(eamKey, eamData, "opt.pot")
		eamKeyBest = eamKey
		eamDataBest = eamData
		If(mpiProcessID.eq.0)Then
	      print *,"RSS: ",i,point,trialResidualSquareSum,x,xV,"*"
		End If  
	  Else	    
		If(mpiProcessID.eq.0)Then
	      print *,"RSS: ",i,point,trialResidualSquareSum,x,xV
		End If  
	  End If	
	End Do
!Test optimum
	Call reReadEamPot("opt.pot")	
!Calculate rss
	eamDataSet = 1		
    Call calcEval()	
	If(mpiProcessID.eq.0)Then
	  print *,"Optimum RSS: ",trialResidualSquareSum
!Store temp potential	  
	  Call storeEAMToFileMaster(eamKeyBest,eamDataBest,"optimum.pot")
	End If  
  End Subroutine runOptimisePointVary


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  Subroutine runOptimisePointVaryWorking(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, fileCount, startI
	Real(kind=DoubleReal) :: bestRSS
	Real(kind=DoubleReal) :: x, xV
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: &
	eamKeyInput,eamKeyOptimised,eamKeyR,eamKeyBest,eamKeyTemp
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: &
	eamDataInput,eamDataOptimised,eamDataR,eamDataBest,eamDataTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: randomSet
	Character(len=64) :: fileName
	Logical :: writeFile
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,processCount,bufferSize
    Integer(kind=StandardInteger) :: processID
	Real(kind=DoubleReal) :: send, receive, buffer	
	Call calcEval()
!----------------------------------------------------------------------------
! Step 1 - store starting potential
!----------------------------------------------------------------------------		
    eamKeyOptimised = eamKeyInput	
    eamDataOptimised = eamDataInput
!----------------------------------------------------------------------------
! Step 2 - vary points
!----------------------------------------------------------------------------		
    startI = 30
	Do i=startI,40
!Synchronise processes		
	  Call synchMpiProcesses()	
	  If(mpiProcessID.eq.0)Then
!Load optimum reduced points
        eamKeyR = eamKeyOptimised	
        eamDataR = eamDataOptimised	 
!Vary the point
	    point = mod(i,size(eamDataR,1))+1
	    x = eamDataR(point,2)
	    xV = VaryPointRand(x,40.0D0,0.3D0)
	    eamDataR(point,2) = xV
!make expanded set of points
	    Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTemp,eamDataTemp)
!Store temp potential	  
	    Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"tempfull.pot")	  
	  End If	
!Synchronise processes		  
	  Call synchMpiProcesses()	
!Load potential from file on to all processes	  
	  Call reReadEamPot("tempfull.pot")	
!Calculate rss
	  eamDataSet = 1		
      Call calcEval()	
!Check if optimum
	  If(i.eq.startI.or.trialResidualSquareSum.lt.bestRSS)Then	
	    bestRSS = trialResidualSquareSum
		eamKeyOptimised = eamKeyR
		eamDataOptimised = eamDataR
		Call storeEAMToFile(eamKey, eamData, "opt.pot")
		eamKeyBest = eamKey
		eamDataBest = eamData
		If(mpiProcessID.eq.0)Then
	      print *,"RSS: ",i,point,trialResidualSquareSum,x,xV,"*"
		End If  
	  Else	    
		If(mpiProcessID.eq.0)Then
	      print *,"RSS: ",i,point,trialResidualSquareSum,x,xV
		End If  
	  End If	
	End Do
!Test optimum
	Call reReadEamPot("opt.pot")	
!Calculate rss
	eamDataSet = 1		
    Call calcEval()	
	If(mpiProcessID.eq.0)Then
	  print *,"Optimum RSS: ",trialResidualSquareSum
!Store temp potential	  
	  Call storeEAMToFileMaster(eamKeyBest,eamDataBest,"optimum.pot")
	End If  
  End Subroutine runOptimisePointVaryWorking

End Module run