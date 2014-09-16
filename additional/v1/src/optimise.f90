Module optimise

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

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
  
!declare global variables  

!Privacy of functions/subroutines/variables
  Private
!Variables  
!Subroutines  
  Public :: optimisePotential   
!Functions
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

 
!------------------------------------------------------------------------!
! Run Optimise Subroutines
!------------------------------------------------------------------------!  
  
  Subroutine optimisePotential()
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
	!Call reReadEamPot(potentialFilePathTemp)	
	Call eamForceZBLCore(eamKey,eamData) 			!If set, force hard zbl core
	Call setPotentialDerivatives(eamKey,eamData) 	!Fill in derivatives
    eamDataSet = 1	!eamData
	saveFileForces = "Y"        !Save forces	
	If(calcRunType.eq.6)Then	!Optimise normal
	  Call calcEval()
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	    print *,ProgramTime(),"First Eval"
	  End If 
	Else If(calcRunType.eq.7)Then	!Optimise full
	  Call calcEvalFull() 
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	    print *,ProgramTime(),"First Eval (Full)"
	  End If  
	End If
	Call calcOutput()	        !Print output to file
	saveFileForces = "N"        !Stop saving forces
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Input Potential RSS:",trialResidualSquareSum  !Starting RSS
	End If 
!----------------------------------------------------------------------------
! Step 2 - make reduced set and evaluate
!----------------------------------------------------------------------------
!Store the input eamReducedPoints and output to files
	Call makeReducedEAMSet(eamKey,eamData,eamKeyR,eamDataR)		!Make reduced set
	Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)	!Make expanded set
	Call eamForceZBLCore(eamKey,eamData) 			            !If set, force hard zbl core
	Call setPotentialDerivatives(eamKey,eamData) 	            !Fill in derivatives
	Call storeEAMToFile(eamKey, eamData, trim(tempDirectory)//"/"//"reduced.pot")			!Store to file	
	Call reReadEamPot(trim(tempDirectory)//"/"//"reduced.pot")								!Re-read potential
	eamDataSet = 1	!eamDataTrial
	If(calcRunType.eq.6)Then	!Optimise normal
	  Call calcEval()
	Else If(calcRunType.eq.7)Then	!Optimise full
	  Call calcEvalFull()  
	End If
	Call calcOutput()	        !Print output to file
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Reduced splined potential RSS:",trialResidualSquareSum  !Starting RSS
	End If  
!----------------------------------------------------------------------------
! Step 3 - make further reduced set and evaluate
!----------------------------------------------------------------------------	
!Remove insignificant points
	!Call runRemoveInsignificantPoints(eamKeyR,eamDataR,potentialResponseResults)	
	Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)	!Make expanded set
	Call eamForceZBLCore(eamKey,eamData) 			            !If set, force hard zbl core
	Call setPotentialDerivatives(eamKey,eamData) 	            !Fill in derivatives
	Call storeEAMToFile(eamKey, eamData, trim(tempDirectory)//"/"//"preOpt.pot")		!Store to file		
	Call reReadEamPot(trim(tempDirectory)//"/"//"preOpt.pot")						!Re-read potential
!Calculate RSS
    eamDataSet = 1	!eamDataTrial
	If(calcRunType.eq.6)Then	!Optimise normal
	  Call calcEval()
	Else If(calcRunType.eq.7)Then	!Optimise full
	  Call calcEvalFull()  
	End If
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Pre-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	End If 
!----------------------------------------------------------------------------
! Step 3 - Run SA optimisation
!----------------------------------------------------------------------------
    Call makeReducedEAMSet(eamKey,eamData,eamKeyR,eamDataR)		!Make reduced set
    Call runOptimisePointVary(eamKeyR,eamDataR,size(eamDataR,1))
!----------------------------------------------------------------------------
! Step 4 - Run SA optimisation
!----------------------------------------------------------------------------	
	If(calcRunType.eq.6)Then	!Optimise normal
	  Call calcEval()
	Else If(calcRunType.eq.7)Then	!Optimise full
	  Call calcEvalFull()  
	End If
	Call calcOutput()	        !Print output to file
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Final RSS: ",trialResidualSquareSum  !Starting RSS
	End If 
  
  End Subroutine optimisePotential 
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
	If(calcRunType.eq.6)Then	!Optimise normal
	  Call calcEval()
	Else If(calcRunType.eq.7)Then	!Optimise full
	  Call calcEvalFull()  
	End If
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
	  If(calcRunType.eq.6)Then	!Optimise normal
	    Call calcEval()
	  Else If(calcRunType.eq.7)Then	!Optimise full
	    Call calcEvalFull()  
	  End If
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
	  If(calcRunType.eq.6)Then	!Optimise normal
	    Call calcEval()
	  Else If(calcRunType.eq.7)Then	!Optimise full
	    Call calcEvalFull()  
	  End If
	  interpPoints(3,1) = eamDataP(i,2)	
	  interpPoints(3,2) = 1.0D0*trialResidualSquareSum
!Third pertubation/fourth point
	  eamKeyP = eamKeyW
	  eamDataP = eamDataW	  
	  eamDataP(i,2) = eamDataP(i,2)+2.0D0*perturbAmount
      Call makeExpandedEAMSet(eamKeyP,eamDataP,eamKeyTrial,eamDataTrial)
	  If(calcRunType.eq.6)Then	!Optimise normal
	    Call calcEval()
	  Else If(calcRunType.eq.7)Then	!Optimise full
	    Call calcEvalFull()  
	  End If
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
  Subroutine runOptimisePointVary(eamKeyRInput,eamDataRInput,dataPoints)
!---------------------------------------------
! 
!
!---------------------------------------------
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, n, point, pointN, fileCount, startI, endI
	Integer(kind=StandardInteger) :: potStart, potEnd, potLength, dataPoints
	Integer(kind=StandardInteger) :: reducedKeyPoints, reducedDataPoints
	Integer(kind=StandardInteger) :: expandedKeyPoints, expandedDataPoints
	Real(kind=DoubleReal) :: bestRSS
	Real(kind=DoubleReal) :: x, xV
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyRInput
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataRInput
	Integer(kind=StandardInteger), Dimension(1:100,1:5) :: eamKeyROpt,eamKeyRVary,eamKeyOptimum
    Real(kind=DoubleReal), Dimension(1:100000,1:2) :: eamDataROpt,eamDataRVary,eamDataOptimum
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: potTypeKey
	Integer(kind=StandardInteger), Dimension(1:dataPoints) :: list
	Logical :: varyPoint
!mpi
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom	
!--------------------------------------------
! Step 1 - Prepare variables and data
!--------------------------------------------
    reducedKeyPoints = size(eamKeyRInput,1)
    reducedDataPoints = size(eamDataRInput,1)
	expandedKeyPoints = size(eamKey,1)
	expandedDataPoints = size(eamData,1)
	eamKeyROpt = 0
	eamKeyRVary = 0
	eamKeyOptimum = 0
	eamDataROpt = 0.0D0
	eamDataRVary = 0.0D0
	eamDataOptimum = 0.0D0
    Do i=1,reducedKeyPoints
	  Do j=1,5
	    eamKeyROpt(i,j) = eamKeyRInput(i,j)
	    eamKeyRVary(i,j) = eamKeyRInput(i,j)
	  End Do
    End Do	  
    Do i=1,reducedDataPoints
	  Do j=1,2
	    eamDataROpt(i,j) = eamDataRInput(i,j)
	    eamDataRVary(i,j) = eamDataRInput(i,j)
	  End Do
    End Do	   
!--------------------------------------------
! Step 2 - Loop through and vary
!--------------------------------------------	
!Make a list of shuffled numbers
    list = 0
    If(mpiProcessID.eq.0)Then
      list = NumberList(dataPoints,500)
	End If
!Send out list
    If(mpiProcessID.eq.0)Then
!SEND by master process	 
      Do i=1,(mpiProcessCount-1)
	    processTo = i
        tag = 2000 + i
		Call MPI_send(list,size(list,1),&
		MPI_integer,processTo,tag,MPI_comm_world,error)
	  End Do	
	Else
!RECV by worker processes
      processFrom = 0
	  tag = 2000 + mpiProcessID
      Call MPI_recv(list,size(list,1),&
	  MPI_integer,processFrom,tag,MPI_comm_world,status,error)
	End If
	startI = 0
	endI = saCycles*size(eamDataRInput,1)
	!endI = 200
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Start optimising, max steps: ",endI
	End If
	Do n=startI,endI
!Load optimum
	  Do i=1,reducedDataPoints
	    Do j=1,2
	      eamDataRVary(i,j) = eamDataROpt(i,j)
	    End Do
      End Do	
!Vary point
      varyPoint = .false. 
	  If(n.gt.startI)Then
	    pointN = mod(n-1,reducedDataPoints)+1
		point = list(pointN)
	    varyPoint = .true.
	  End If	  
	  If(varyPoint.eqv..true.)Then
!Vary on master, share to workers
	    If(mpiProcessID.eq.0)Then
		  x = eamDataRVary(point,2)
	      xV = VaryPointRand(x,saTemp,saSpreadFactor,1.2D0,0.01D0,1) 
	      eamDataRVary(point,2) = xV		  
        End If  
	  End If
	  Call MPI_sendData2DDP(eamDataRVary,size(eamDataRVary,1),size(eamDataRVary,2))
!Copy vary array back
      Do i=1,reducedDataPoints
	    Do j=1,2
	      eamDataRInput(i,j) = eamDataRVary(i,j)
	    End Do
      End Do	
	  Call synchMpiProcesses()	
      If(mpiProcessID.eq.0)Then
	    Call makeExpandedEAMSet(eamKeyRInput,eamDataRInput,eamKey,eamData)                !Expand 
	    Call eamForceZBLCore(eamKey, eamData) 
      End If
	  !Call MPI_sendData2DInt(eamKey,size(eamKey,1),size(eamKey,2))
	  Call MPI_sendData2DDP(eamData,size(eamData,1),size(eamData,2))
!Calculate rss
	  eamDataSet = 1		
	  If(calcRunType.eq.6)Then	!Optimise normal
	    Call calcEval()
	  Else If(calcRunType.eq.7)Then	!Optimise full
	    Call calcEvalFull()  
	  End If
!Store first/optimum 	  
	  If(n.eq.0)Then
        bestRSS = trialResidualSquareSum
		If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	      print *,ProgramTime(),"Start RSS: ",trialResidualSquareSum
		End If 
      Else
	    If(trialResidualSquareSum.lt.bestRSS)Then	
		  bestRSS = trialResidualSquareSum
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",n,pointN,point,x,xV,trialResidualSquareSum,"****"
	      End If 
!store optimum reduced
          If(mpiProcessID.eq.0)Then
	        Do i=1,reducedDataPoints
	          Do j=1,2
	            eamDataROpt(i,j) = eamDataRVary(i,j)
	          End Do
            End Do	
		  End If	
		  Call MPI_sendData2DDP(eamDataROpt,size(eamDataROpt,1),size(eamDataROpt,2))
!store optimum full	
          If(mpiProcessID.eq.0)Then
		    Do i=1,size(eamKey,1)
	          Do j=1,5
	            eamKeyOptimum(i,j) = eamKey(i,j)
	          End Do
            End Do	
	        Do i=1,size(eamData,1)
	          Do j=1,2
	            eamDataOptimum(i,j) = eamData(i,j)
	          End Do
            End Do	
		  End If		  
	      !Call MPI_sendData2DInt(eamKeyOptimum,size(eamKeyOptimum,1),size(eamKeyOptimum,2))
		  Call MPI_sendData2DDP(eamDataOptimum,size(eamDataOptimum,1),size(eamDataOptimum,2))
	    Else	   
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",n,pointN,point,x,xV,trialResidualSquareSum
	      End If  
	    End If	
	  End If		  
	End Do
!--------------------------------------------
! Step 3 - Store optimum
!--------------------------------------------	
	Do i=1,size(eamData,1)
	  Do j=1,2
	    eamData(i,j) = eamDataOptimum(i,j)
	  End Do
    End Do	
	Call storeEAMToFileMaster(eamKey,eamData,trim(outputDirectory)//"/optimum.pot")
!--------------------------------------------
! Step 4 - Check and output result
!--------------------------------------------	
	Call reReadEamPot(trim(outputDirectory)//"/optimum.pot")						!Re-read potential
!Calculate RSS
    eamDataSet = 1	!eamDataTrial
	If(calcRunType.eq.6)Then	!Optimise normal
	  Call calcEval()
	Else If(calcRunType.eq.7)Then	!Optimise full
	  Call calcEvalFull()  
	End If
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Optimised RSS: ",trialResidualSquareSum  !Starting RSS
	End If 
	
  End Subroutine runOptimisePointVary


  
  
  
  
  

End Module optimise