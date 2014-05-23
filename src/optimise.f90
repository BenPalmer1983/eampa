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
	Call reReadEamPot(potentialFilePathTemp)	
	Call eamForceZBLCore(eamKey,eamData) 			!If set, force hard zbl core
	Call setPotentialDerivatives(eamKey,eamData) 	!Fill in derivatives
    eamDataSet = 1	!eamData
	Call calcEval()
	Call calcOutput()	        !Print output to file
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
	Call storeEAMToFile(eamKey, eamData, "temp.pot")			!Store to file	
	Call reReadEamPot("temp.pot")								!Re-read potential
	eamDataSet = 1	!eamDataTrial
	Call calcEval()
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
	Call storeEAMToFile(eamKey, eamData, "temp/preOpt.pot")		!Store to file	
	Call reReadEamPot("temp/preOpt.pot")						!Re-read potential
!Calculate RSS
    eamDataSet = 1	!eamDataTrial
	Call calcEval()
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Pre-optimise (SA1) potential RSS: ",trialResidualSquareSum  !Starting RSS
	End If 
!----------------------------------------------------------------------------
! Step 3 - Run SA optimisation
!----------------------------------------------------------------------------
    Call makeReducedEAMSet(eamKey,eamData,eamKeyR,eamDataR)		!Make reduced set
    Call runOptimisePointVary(eamKeyR,eamDataR)
	
  
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
  Subroutine runOptimisePointVary(eamKeyRInput,eamDataRInput)
!---------------------------------------------
! 
!
!---------------------------------------------
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, n, point, fileCount, startI, endI
	Integer(kind=StandardInteger) :: potStart, potEnd, potLength
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
	Logical :: varyPoint
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
	    point = mod(n,reducedDataPoints)+1
	    varyPoint = .true.
	  End If	  
	  If(varyPoint.eqv..true.)Then
!Vary on master, share to workers
	    If(mpiProcessID.eq.0)Then
		  x = eamDataRVary(point,2)
	      xV = VaryPointRand(x,saTemp,saSpreadFactor) 
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
	  Call MPI_sendData2DInt(eamKey,size(eamKey,1),size(eamKey,2))
	  Call MPI_sendData2DDP(eamData,size(eamData,1),size(eamData,2))
!Calculate rss
	  eamDataSet = 1		
      Call calcEval()
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
	        print *,ProgramTime(),"RSS: ",n,point,x,xV,trialResidualSquareSum,"****"
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
	      Call MPI_sendData2DInt(eamKeyOptimum,size(eamKeyOptimum,1),size(eamKeyOptimum,2))
		  Call MPI_sendData2DDP(eamDataOptimum,size(eamDataOptimum,1),size(eamDataOptimum,2))
	    Else	   
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",n,point,x,xV,trialResidualSquareSum
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
	Call storeEAMToFileMaster(eamKey,eamData,"opt/optimum.pot")
	
	
	
	
	
	

  End Subroutine runOptimisePointVary


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    Subroutine runOptimisePointVary3(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
!---------------------------------------------
! 
!
!---------------------------------------------
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, fileCount, startI, endI
	Integer(kind=StandardInteger) :: potStart, potEnd, potLength
	Real(kind=DoubleReal) :: bestRSS
	Real(kind=DoubleReal) :: x, xV
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: &
	eamKeyInput,eamKeyOptimised,eamKeyR,eamKeyBest,eamKeyTemp
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: &
	eamDataInput,eamDataOptimised,eamDataR,eamDataBest,eamDataTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: potTypeKey
	Character(len=64) :: fileName, headerText
	Logical :: writeFile, varyPoint
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,processCount,bufferSize
    Integer(kind=StandardInteger) :: processID
	Real(kind=DoubleReal) :: send, receive, buffer	
!----------------------------------------------------------------------------
! Step 1 - store starting potential
!----------------------------------------------------------------------------	
    eamKeyR = eamKeyInput
	eamDataR = eamDataInput
	Call storeReducedEAM(eamKeyR,eamDataR,"temp/optReduced.tmp")
!----------------------------------------------------------------------------
! Step 2 - make potTypeKey
!----------------------------------------------------------------------------	
    Allocate(potTypeKey(1:size(eamDataInput,1)))
	Do i=1,size(eamKeyInput,1)
      potStart = eamKeyInput(i,4)
      potLength = eamKeyInput(i,5)
	  potEnd = potStart + potLength - 1
!pot type
	  potType = eamKeyInput(i,3) !1 PAIR, 2 DENS/SDEN, 3 EMBE/SEMB, 4 DDEN, 5 DEMB
	  Do j=potStart,potEnd
	    potTypeKey(j) = potType
	  End Do
	End Do	
!----------------------------------------------------------------------------
! Step 3 - vary points
!----------------------------------------------------------------------------		
    startI = 0
	!endI = saCycles*size(eamDataInput,1)
	endI = 100
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Start optimising, max steps: ",endI
	End If
	Do i=startI,endI
!Synchronise processes		
	  Call synchMpiProcesses()	
!Load optimum reduced points
	  Call loadReducedEAM(eamKeyR,eamDataR,"temp/optReduced.tmp")  
	  
!check load ok
!Prepare potential
	  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)                !Expand 
	  Call storeEAMToFileMaster(eamKey,eamData,"temp/tempfull.pot")	          !Save
	  Call synchMpiProcesses()	                                              !Synch 
	  Call reReadEamPot("temp/tempfull.pot")	                              !Reload
	  Call synchMpiProcesses()	                                              !Synch 
!Calculate rss
	  eamDataSet = 1		
      Call calcEval()	
	  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	      print *,ProgramTime(),"Test RSS: ",trialResidualSquareSum
	  End If 
	  
!Vary point
      varyPoint = .false. 
	  If(i.gt.startI)Then
	    point = mod(i,size(eamDataR,1))+1
	    varyPoint = .true.
	  End If	  
	  If(varyPoint.eqv..true.)Then
		x = eamDataR(point,2)
!Vary on master, share to workers
	    If(mpiProcessID.eq.0)Then
	      xV = VaryPointRand(x,40.0D0,0.0001D0) 
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
	  End If
!Prepare potential
	  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)                !Expand 
	  Call storeEAMToFileMaster(eamKey,eamData,"temp/tempfull.pot")	          !Save
	  Call synchMpiProcesses()	                                              !Synch 
	  Call reReadEamPot("temp/tempfull.pot")	                              !Reload
	  Call synchMpiProcesses()	                                              !Synch 
!Calculate rss
	  eamDataSet = 1		
      Call calcEval()		  
!If first iteration, store best
	  If(i.eq.0)Then
        bestRSS = trialResidualSquareSum
		If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	      print *,ProgramTime(),"Start RSS: ",trialResidualSquareSum
		End If 
      Else
	    If(trialResidualSquareSum.lt.bestRSS)Then	
	      Call synchMpiProcesses()	
		  Call storeReducedEAM(eamKeyR,eamDataR,"temp/optReduced.tmp")
		  bestRSS = trialResidualSquareSum
	      Call synchMpiProcesses()	
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",i,point,x,xV,trialResidualSquareSum,"****"
	      End If 
		  
!Check 1
	      eamDataSet = 1		
          Call calcEval()
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"Check1 RSS: ",trialResidualSquareSum,"****"
	      End If

!Check 2                        
	      Call synchMpiProcesses()	 
	      Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)                !Expand 
	      Call storeEAMToFileMaster(eamKey,eamData,"temp/tempfull.pot")	          !Save
	      Call synchMpiProcesses()	                                              !Synch 
	      Call reReadEamPot("temp/tempfull.pot")	                              !Reload
	      Call synchMpiProcesses()	 
	      eamDataSet = 1		
          Call calcEval()	
	      If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"Check2 RSS: ",trialResidualSquareSum,"****"
	      End If 
		  
		  
!Check 3              
	      Call synchMpiProcesses()	 
	      Call loadReducedEAM(eamKeyR,eamDataR,"temp/optReduced.tmp")  		  !Synch 
		  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)                !Expand 
	      Call storeEAMToFileMaster(eamKey,eamData,"temp/tempfull.pot")	          !Save
	      Call synchMpiProcesses()	                                              !Synch 
	      Call reReadEamPot("temp/tempfull.pot")	                              !Reload
	      Call synchMpiProcesses()	 
	      eamDataSet = 1		
          Call calcEval()	
	      If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"Check3 RSS: ",trialResidualSquareSum,"****"
	      End If 
		  
		  
	    Else	   
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",i,point,x,xV,trialResidualSquareSum
	      End If  
	    End If	
	  End If		  
	End Do
!----------------------------------------------------------------------------
! Step 4 - test optimum
!----------------------------------------------------------------------------		
!Synchronise processes		
	Call synchMpiProcesses()	
!Load optimum reduced points
	Call loadReducedEAM(eamKeyR,eamDataR,"temp/optReduced.tmp")  
!Prepare potential
	Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKey,eamData)                !Expand 
	Call storeEAMToFileMaster(eamKey,eamData,"opt/optimum.pot")	            !Save
	Call synchMpiProcesses()	                                            !Synch 
	Call reReadEamPot("opt/optimum.pot")	                                !Reload
	Call synchMpiProcesses()	                                            !Synch 
!Calculate rss
	eamDataSet = 1		
    Call calcEval()	
	If(mpiProcessID.eq.0)Then
	  print *,"Optimum RSS Check: ",trialResidualSquareSum
	End If 
	
	
!Test optimum
	!Call makeExpandedEAMSet(eamKeyOptimised,eamDataOptimised,eamKey,eamData)
!Calculate rss
	!eamDataSet = 1		
    !Call calcEval()	
	!If(mpiProcessID.eq.0)Then
	!  print *,"Optimum RSS: ",trialResidualSquareSum
	!End If  
!Store temp potential	  
	!Call storeEAMToFileMaster(eamKey,eamData,"opt/optimum.pot")
	!Call reReadEamPot("opt/optimum.pot")	
	!eamDataSet = 1	
    !Call calcEval()	
	!If(mpiProcessID.eq.0)Then
	!  print *,"Optimum RSS Check: ",trialResidualSquareSum
	!End If 	
	
	
	!Call reReadEamPot("temp/opt.pot")	
	!Call eamForceZBLCore(eamKey,eamData) 			!If set, force hard zbl core
	!Call setPotentialDerivatives(eamKey,eamData) 	!Fill in derivatives
!Calculate rss
	!eamDataSet = 1		
    !Call calcEval()	
	!If(mpiProcessID.eq.0)Then
	!  print *,"Optimum RSS: ",trialResidualSquareSum
!Store temp potential	  
	!  Call storeEAMToFileMaster(eamKeyBest,eamDataBest,"opt/optimum.pot")
	!End If  
  End Subroutine runOptimisePointVary3

  
  
   Subroutine runOptimisePointVaryWorking2(eamKeyInput,eamDataInput,eamKeyOptimised,eamDataOptimised)
!---------------------------------------------
! 
!
!---------------------------------------------
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, point, fileCount, startI, endI
	Integer(kind=StandardInteger) :: potStart, potEnd, potLength
	Real(kind=DoubleReal) :: bestRSS
	Real(kind=DoubleReal) :: x, xV
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: &
	eamKeyInput,eamKeyOptimised,eamKeyR,eamKeyBest,eamKeyTemp
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: &
	eamDataInput,eamDataOptimised,eamDataR,eamDataBest,eamDataTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: potTypeKey
	Character(len=64) :: fileName
	Logical :: writeFile, varyPoint
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
!Store input potential	  
    Call storeEAMToFileMaster(eamKeyInput,eamDataInput,"temp/startingReduced.pot")	
	Call makeExpandedEAMSet(eamKeyInput,eamDataInput,eamKeyTemp,eamDataTemp)
    Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"temp/startingFull.pot")	
!----------------------------------------------------------------------------
! Step 2 - make potTypeKey
!----------------------------------------------------------------------------	
    Allocate(potTypeKey(1:size(eamDataInput,1)))
	Do i=1,size(eamKeyInput,1)
      potStart = eamKeyInput(i,4)
      potLength = eamKeyInput(i,5)
	  potEnd = potStart + potLength - 1
!pot type
	  potType = eamKeyInput(i,3) !1 PAIR, 2 DENS/SDEN, 3 EMBE/SEMB, 4 DDEN, 5 DEMB
	  Do j=potStart,potEnd
	    potTypeKey(j) = potType
	  End Do
	End Do	
!----------------------------------------------------------------------------
! Step 3 - vary points
!----------------------------------------------------------------------------		
    startI = 0
	endI = saCycles*size(eamDataInput,1)
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Start optimising, max steps: ",endI
	End If
	Do i=startI,endI
!Synchronise processes		
	  Call synchMpiProcesses()	
!Load optimum reduced points
      eamKeyR = eamKeyOptimised	
      eamDataR = eamDataOptimised	
      Call synchMpiProcesses()
!-----check optimum
	  Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTemp,eamDataTemp)
!Store temp potential	  
	  Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"temp/tempfullb.pot")	  
!Synchronise processes		  
	  Call synchMpiProcesses()	
!Load potential from file on to all processes	  
	  Call reReadEamPot("temp/tempfullb.pot")	
	  eamDataSet = 1		
      Call calcEval()	
	  If(mpiProcessID.eq.0)Then
		print *,ProgramTime(),"RSS check: ",i,point,trialResidualSquareSum
	  End If
!-----check optimum
	  
!Set best rss
	  If(i.eq.startI)Then
!make expanded set of points
	    Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTemp,eamDataTemp)
!Store temp potential	  
	    Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"temp/tempfull.pot")	
	    Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"opt/optimum.pot")	  
!Synchronise processes		  
	    Call synchMpiProcesses()	
!Load potential from file on to all processes	  
	    Call reReadEamPot("temp/tempfull.pot")	
!Calculate rss
	    eamDataSet = 1		
        Call calcEval()
		bestRSS = trialResidualSquareSum
		If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	      print *,ProgramTime(),"Start RSS: ",trialResidualSquareSum
		End If  
	  End If
	  
	  
	  
!Select the point
	  point = mod(i,size(eamDataR,1))+1
!Only vary the point if NOT in ZBL forced section of function	  
      varyPoint = .false.
	  If(potTypeKey(point).eq.1.and.eamDataR(point,1).ge.eamZBLPairUpper)Then
	    varyPoint = .true. 	!Check, if pair function
	  End If
	  If((potTypeKey(point).eq.2.or.potTypeKey(point).eq.4).and.&
	  eamDataR(point,1).ge.eamZBLDensCutoff)Then
	    varyPoint = .true. 	!Check, if dens function
	  End If
	  If((potTypeKey(point).eq.3.or.potTypeKey(point).eq.5).and.&
	  eamDataR(point,1).ge.eamZBLEmbeCutoff)Then
	    varyPoint = .true. 	!Check, if embe function
	  End If
	  If(varyPoint.eqv..true.)Then
	    x = eamDataR(point,2)
!Vary on master, share to workers
	    If(mpiProcessID.eq.0)Then
	      xV = VaryPointRand(x,40.0D0,0.001D0) 
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
	    !Call eamForceZBLCore(eamKey,eamData) 			!If set, force hard zbl core
	    !Call setPotentialDerivatives(eamKey,eamData) 	!Fill in derivatives
	    Call makeExpandedEAMSet(eamKeyR,eamDataR,eamKeyTemp,eamDataTemp)
!Store temp potential	  
	    Call storeEAMToFileMaster(eamKeyTemp,eamDataTemp,"temp/tempfull.pot")	  
!Synchronise processes		  
	    Call synchMpiProcesses()	
!Load potential from file on to all processes	  
	    Call reReadEamPot("temp/tempfull.pot")	
!Calculate rss
	    eamDataSet = 1		
        Call calcEval()	
!Check if optimum
	    If(trialResidualSquareSum.lt.bestRSS)Then	
	      bestRSS = trialResidualSquareSum
		  eamKeyOptimised = eamKeyR
		  eamDataOptimised = eamDataR
		  Call storeEAMToFile(eamKey, eamData, "temp/opt.pot")
		  eamKeyBest = eamKey
		  eamDataBest = eamData
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",i,point,x,xV,trialResidualSquareSum,"****"
	      End If 
	    Else	   
		  If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	        print *,ProgramTime(),"RSS: ",i,point,x,xV,trialResidualSquareSum
	      End If  
	    End If	
	  End If	
	End Do
!Test optimum
	Call reReadEamPot("temp/opt.pot")	
	!Call eamForceZBLCore(eamKey,eamData) 			!If set, force hard zbl core
	!Call setPotentialDerivatives(eamKey,eamData) 	!Fill in derivatives
!Calculate rss
	eamDataSet = 1		
    Call calcEval()	
	If(mpiProcessID.eq.0)Then
	  print *,"Optimum RSS: ",trialResidualSquareSum
!Store temp potential	  
	  Call storeEAMToFileMaster(eamKeyBest,eamDataBest,"opt/optimum.pot")
	End If  
  End Subroutine runOptimisePointVaryWorking2
  
  
  
  
  
  
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
	  Call storeEAMToFileMaster(eamKeyBest,eamDataBest,"opt/optimum.pot")
	End If  
  End Subroutine runOptimisePointVaryWorking



  
  

End Module optimise