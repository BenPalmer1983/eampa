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
	  Call calcConfigEnergies()	!Calculate energies	  
	  Call calcOutput()	        !Print output to file
	End If
	
	If(calcRunType.eq.5)Then	!Evaluate	  
	  eamDataSet = 3 	  
	  Call calcEval()  
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
	Integer(kind=StandardInteger) :: i, j, k, point
	Real(kind=DoubleReal) :: startingRSS, bestRSS
	Real(kind=DoubleReal) :: varyAmount, unperturbed, perturbAmount, y, dy, ddy
	Integer(kind=StandardInteger) :: varyType
	Integer(kind=StandardInteger), Dimension(:), Allocatable :: potType
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: interpPoints, makeJacobianMatrix
	Real(kind=DoubleReal), Dimension(:), Allocatable :: yArray
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: jacobianMatrix, jacobianMatrixT
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: hessianMatrix
  
 !Set to use trial data points in calculations
    eamDataSet = 3
	varyType = 3 !1 = all points 2 = individual points
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
	  Do i=1,5200
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
		print *,i,varyAmount,trialResidualSquareSum
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
	
	
	
	
!Run final test
    Call loadOptToReduced()
	Call makeTrialEAMSet()
	Call calcEval()
	
	Print *,bestRSS
	
  
  
  End Subroutine runOptimise 
  
  
  
  
  
  

End Module run