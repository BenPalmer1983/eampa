Module maths

!--------------------------------------------------------------!
! Maths functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 10th March 2014
!----------------------------------------


! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
!Make private
  Private  

! General Maths Functions
  Public :: Factorial  
  Public :: BinomialCoefficient
! Polynomial Related Functions
  Public :: SolvePolynomial  
  Public :: CalcPolynomial  
! Interpolation and Regression Functions 
  Public :: PolyFit  
  Public :: PolynomialInterpolation  
  Public :: QuadraticInterpolation  
  Public :: QuadraticInterpolationCalc  
  Public :: CubicInterpolationCalc  
  Public :: CalcResidualSquareSum  
  Public :: PointInterpolation
  Public :: PointInterpolationArr
  Public :: CalcSecondDerivativeMinimum
! Spline Functions   
  Public :: Spline
  Public :: MatrixInterpolation
! Matrix Functions  
  Public :: InvertSquareMatrix
! Random Number Related Functions  
  Public :: SetRandomSeedArray !Subroutine
  Public :: VaryPoint
  Public :: VaryPoints
! Laplace Transforms  
  Public :: GaverStehfestCoeffs
! Decay Equations    
  Public :: CalcIsotopeAmount  
! Rounding functions
  Public :: Ceil     
! Useful Miscellaneous Functions
! Misc Subroutines  

  
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 
! -
! 
! General Maths Functions
! - Factorial
!
! 
! Polynomial Related Functions
! - SolvePolynomial
! - CalcPolynomial
!  
!  
! Interpolation and Regression Functions 
! - PolyFit
! - PolynomialInterpolation 
! - PolynomialRegression
! - PolynomialRegressionVandermonde
! - QuadraticInterpolation
! - QuadraticInterpolationCalc
! - CubicInterpolationCalc
! - CalcResidualSquareSum
! - CalcSecondDerivativeMinimum
!  
!  
! Spline Functions 
! - FivePointInterpolation
! - FivePointInterpolationValue
!
!
! Matrix Functions
! - InvertSquareMatrix
! - SolveLinear
! 
! Random Number Related Functions
! - RandomSeed
! - RandomDataPoints
! - VaryPoint
!
! Laplace Transforms
! - InverseLaplaceTransform
! 
! Decay Equations  
! - Decay Functions
! 
! Useful Miscellaneous Functions
! - change ArraySize1D
! - change ArraySize2D

  
  
!------------------------------------------------------------------------!
! General Maths Functions
!------------------------------------------------------------------------! 

  Function Factorial(input) RESULT (output)
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: output
!calculate factorial
    output = 1
	Do i=1,input
	  output = i * output
	End Do  
  End Function Factorial
  
  Function FactorialQ(input) RESULT (output)
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,input
    Real(kind=QuadrupoleReal) :: output
!calculate factorial
    output = 1.0D0
	Do i=1,input
	  output = i * output
	End Do  
  End Function FactorialQ
  
  Function BinomialCoefficient(n,k) RESULT (c)
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: c,n,k
!calculate factorial
    c = Factorial(n)/(Factorial(n-k)*Factorial(k))
  End Function BinomialCoefficient
  
    
!------------------------------------------------------------------------!
! Polynomial Related Functions
!------------------------------------------------------------------------!     
  
  Function SolvePolynomial (coefficients, lower, upper) RESULT (output)    	
!force declaration of all variables
	Implicit None
!declare variables
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: upper, lower, output
	Real(kind=DoubleReal) :: x,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget, factor, difference
	Integer(kind=StandardInteger) :: i,j,k,maxLoops	
!Set values
	convergenceTarget = 0
	convergence = 1000
	convergenceThreshold = 0.00001
	maxLoops = 0
!set start value for x
	factor = 0.5
	difference = upper - lower
	x = lower + factor * difference	
	do while(convergence.gt.convergenceThreshold.and.maxLoops.le.10000)
	  maxLoops = maxLoops + 1
	  difference = factor * difference
	  y = 0
	  do i=0,size(coefficients)-1
	    y = y + x**(i) * coefficients(i)			
	  enddo
	  dydx = 0
	  do i=1,size(coefficients)-1
	    dydx = dydx + i * x**(i-1) * coefficients(i)			
	  enddo		  
	  convergence = abs(convergenceTarget - y)
	  if(convergence.gt.convergenceThreshold)then
	    if((dydx.lt.0.and.y.ge.0).or.(dydx.ge.0.and.y.lt.0))then
	      x = x + difference	
	    else
	      x = x - difference
	    endif
	  endif
	enddo	
	output = x
  End Function solvePolynomial   
  
   
  Function CalcPolynomial (polyCoefficients, x, derivative) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k, derivative
	Integer(kind=StandardInteger) :: coeffCount 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyWorking
    Real(kind=DoubleReal) :: x, y
!count coefficients
	coeffCount = size(polyCoefficients)	
!force no or positive derivative
	if(derivative.lt.0)then
	  derivative = 0
	endif
!run calc
	if(derivative.gt.coeffCount)then
!if order of derivative higher than coeffs/order polynomial, then answer is zero	
	  y = 0.0D0
	  Allocate(polyWorking(0:0))
!transfer to polyWorking
	  polyWorking(0) = 0.0D0
	  coeffCount = 1
	else
	  if(derivative.eq.0)then
!if no derivative, use input coeffs	
        Allocate(polyWorking(0:(coeffCount-1)))
!transfer to polyWorking
	    do i=0,(coeffCount-1)
		  polyWorking(i) = polyCoefficients(i)
		enddo
	  else
!loop through each derivative
	    do i=1,derivative
		  coeffCount = coeffCount - 1
		  do j=0,(coeffCount-1)
		    polyCoefficients(j) = (j+1)*polyCoefficients(j+1)
		  enddo
		enddo	
        Allocate(polyWorking(0:(coeffCount-1)))
!transfer to polyWorking
	    do i=0,(coeffCount-1)
		  polyWorking(i) = polyCoefficients(i)
		enddo
	  endif
	  y = 0.0D0
	  do i=0,(coeffCount-1)
	    y = y+1.0D0*polyWorking(i)*x**(1.0D0*i)
	  enddo 
	endif 
!returns y
  End Function CalcPolynomial
  
  
  
  
  
  Function PolynomialCoefficientDerivative (coefficientsIn) RESULT (coefficientsOut)
!Take array of polynomial coefficients and return derivative polynomial coefficients
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,order
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsIn, coefficientsOut
!Calculate
    order = size(coefficientsIn)-1
!Make coefficients
    If(order.lt.1)Then
!Allocate array
      Allocate(coefficientsOut(0:(order-1)))
	  coefficientsOut(0) = 0.0D0
	ElseIf(order.ge.1)Then
!Allocate array
      Allocate(coefficientsOut(0:(order-1)))
	  Do i=0,(order-1)
	    coefficientsOut(i) = (i+1)*coefficientsIn(i+1)
	  End Do
	End If  
  End Function PolynomialCoefficientDerivative
  
  
!------------------------------------------------------------------------!
! Interpolation and Regression Functions  
!------------------------------------------------------------------------! 
  
  Function PolyFit(inputPoints, order) RESULT (coefficients) 
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,order
	Real(kind=DoubleReal) :: rssThreshold, rssA, rssB, rssC, rssBest
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsA
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsB
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsTemp
    Integer(kind=StandardInteger), Dimension(1:1000) :: randomSeed
!set random seed
    !randomSeed = SetRandomSeedArray()
	Call SetRandomSeedArray()
!allocate arrays
    Allocate(coefficients(0:order))
    Allocate(coefficientsA(0:order))
    Allocate(coefficientsB(0:order))
    Allocate(coefficientsTemp(0:order))
!Use vandermonde method and calculate rss
	coefficientsA = PolynomialRegressionVandermonde(inputPoints, order)
	rssA = CalcResidualSquareSum(inputPoints,coefficientsA)
!calculate and compare rss for each
    rssThreshold = rssA
	coefficientsB = PolynomialRegression(inputPoints, order, rssThreshold)
	rssB = CalcResidualSquareSum(inputPoints,coefficientsB)
	if(rssA.lt.rssB)then
!use vandermonde result
      do i=0,order
	    coefficients(i) = coefficientsA(i) 
	  enddo
    else
!improve regression attempt if better than vandermonde method	  
	  rssBest = rssB
	  !print *,rssB
	  do i=1,4
	    rssThreshold = rssThreshold*0.2**i
		coefficientsTemp = PolynomialRegression(inputPoints, order, rssThreshold)
		rssB = CalcResidualSquareSum(inputPoints,coefficientsTemp)
		if(rssB.lt.rssBest)then
		  do j=0,order
		    coefficientsB(j) = coefficientsTemp(j)
			rssBest = rssB
		  enddo	
		else
          exit
        endif		  
	  enddo
	  do i=0,order
	    coefficients(i) = coefficientsB(i) 
	  enddo
	endif
  End Function PolyFit      
  
  
  
  Function PolynomialInterpolation(inputPoints) RESULT (coefficients)   
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i, col, row
	Integer(kind=StandardInteger) :: order, coefficientCount
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsTemp
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix
!set variables
    coefficientCount = size(inputPoints,1)
	order = size(inputPoints,1)-1
!Allocate arrays
    Allocate(xMatrix(1:coefficientCount,1:coefficientCount))
    Allocate(xMatrixInverse(1:coefficientCount,1:coefficientCount))
    Allocate(yMatrix(1:coefficientCount))
    Allocate(coefficientsTemp(1:coefficientCount))
    Allocate(coefficients(0:order))
!Fill arrays
	do row=1,coefficientCount
	  do col=1,coefficientCount
        xMatrix(row,col) = 1.0D0*inputPoints(row,1)**(1.0D0*(col-1))
      enddo
	  yMatrix(row) = 1.0D0*inputPoints(row,2)
    enddo
!invert xMatrix
    xMatrixInverse = InvertSquareMatrix(xMatrix)
!find coefficients
    coefficientsTemp = matMul(xMatrixInverse,yMatrix)
!move coefficients
    do i=1,coefficientCount
      coefficients(i-1) = coefficientsTemp(i)
	enddo
  End Function PolynomialInterpolation  
    
	
	
  Function PolynomialRegression(inputPoints, order, rssThresholdIn) RESULT &
  (coefficients)   
!force declaration of all variables
	Implicit None
!declare variables
    Real(kind=DoubleReal),optional :: rssThresholdIn
    Integer(kind=StandardInteger) :: i,j,k,col,row
	Integer(kind=StandardInteger) :: order, coefficientCount, polyCount, loops
	Real(kind=DoubleReal) :: randNumber, rss, bestRss, decayFactor, rssThreshold
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: samplePoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: tempCoefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: testCoefficients
!set variables	
	coefficientCount = order + 1
	if(present(rssThresholdIn))then
	  rssThreshold = rssThresholdIn
	else
	  rssThreshold = 1000
	endif
	if(rssThreshold.lt.1)then
	  !rssThreshold = 1
	endif
!allocate arrays
    Allocate(coefficients(0:order))
    Allocate(tempCoefficients(0:order))
    Allocate(testCoefficients(0:order))
    Allocate(samplePoints(1:coefficientCount,1:2))
!fill coefficients array with zeros
    do j=0,order
	  testCoefficients(j) = 0.0D0
	  coefficients(j) = 0.0D0
	  tempCoefficients(j) = 0.0D0
	enddo
!find best fit by interpolation between random data points
    polyCount = 0
	if(order.lt.3)then
	  loops = 25
	else  
	  loops = 50*((order-2)**2)
	endif
	do i=1,(48+order**3)
	  samplePoints = RandomDataPoints(inputPoints,coefficientCount,1)
	  tempCoefficients = PolynomialInterpolation(samplePoints)
	  rss = CalcResidualSquareSum(inputPoints,tempCoefficients)
	  if(rss.lt.rssThreshold)then 
	    if(polyCount.eq.0)then
	      do j=0,order
		    coefficients(j) = tempCoefficients(j)
		  enddo
		  bestRss = rss		
		  polyCount = polyCount + 1  
		else 
	      do j=0,order		
		    testCoefficients(j) = &
			   (polyCount * coefficients(j) + tempCoefficients(j)) /&
			   (polyCount + 1)
		  enddo
		  rss = CalcResidualSquareSum(inputPoints,testCoefficients)
		  if(rss.lt.bestRss)then
		    bestRss = rss	
		    polyCount = polyCount + 1
	        do j=0,order
		      coefficients(j) = 1.0D0*testCoefficients(j)
		    enddo
		  endif		  		  
		endif	  
	  endif	
	enddo
  End Function PolynomialRegression  
  

  
  Function PolynomialRegressionVandermonde(inputPoints, order)&
  RESULT (coefficients) 
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,k,col,row
	Integer(kind=StandardInteger) :: order, dataPoints, matrixSize, exponentValue
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
!set variables
    dataPoints = size(inputPoints,1)
	matrixSize = order+1
!allocate arrays
    Allocate(xMatrix(1:matrixSize,1:matrixSize))
    Allocate(yMatrix(1:matrixSize))
    Allocate(coefficients(0:order))
!Build Least Squares Fitting Vandermonde matrix
    do row=1,matrixSize
	  do col=1,matrixSize
	    exponentValue = 1.0D0*row+1.0D0*col-2.0D0
		xMatrix(row,col) = 0.0D0
		do k=1,dataPoints
	      xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*inputPoints(k,1)&
		  **exponentValue
		enddo
	  enddo
	enddo
	do row=1,matrixSize
	  exponentValue = 1.0D0*row-1.0D0
	  yMatrix(row) = 0.0D0
	  do k=1,dataPoints
	    yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*inputPoints(k,2)*&
		inputPoints(k,1)**exponentValue
	  enddo
	enddo
!invert xMatrix
	xMatrix = InvertSquareMatrix(xMatrix)
!multiply inverse by y to get coefficients
    yMatrix = matMul(xMatrix,yMatrix)
!save coefficients
    do i=0,order
	  coefficients(i) = yMatrix(i+1)
	enddo  
  End Function PolynomialRegressionVandermonde  
  
  
  
  Function QuadraticInterpolation(xA,yA,xB,yB,xC,yC) RESULT (coefficients)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC
	Real(kind=DoubleReal) :: aA,aB,aC
	Allocate(coefficients(0:2))  	
	aA = 1.0D0*(yA/((xA-xB)*(xA-xC)))
    aB = 1.0D0*(yB/((xB-xA)*(xB-xC)))
	aC = 1.0D0*(yC/((xC-xA)*(xC-xB)))
	coefficients(2) = 1.0D0*(aA+aB+aC)
	coefficients(1) = -1.0D0*(aA*(xB+xC)+aB*(xA+xC)+aC*(xA+xB))
	coefficients(0) = 1.0D0*(aA*xB*xC+aB*xA*xC+aC*xA*xB)
  End Function QuadraticInterpolation
  !------------------------------
  Function QuadraticInterpolationCalc(xA,yA,xB,yB,xC,yC,x) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC,x,y
	Allocate(coefficients(0:2)) 
	coefficients = QuadraticInterpolation(xA,yA,xB,yB,xC,yC)
	y = 1.0D0*coefficients(0)+1.0D0*coefficients(1)*x+1.0D0*coefficients(2)*x**2
  End Function QuadraticInterpolationCalc
  
  
  
  Function CubicInterpolationCalc(xA,yA,xB,yB,xC,yC,xD,yD,x) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC,xD,yD,x,y
	Real(kind=DoubleReal) :: aA,aB,aC,aD
	
	aA = 1.0D0*yA*(((x-xB)*(x-xC)*(x-xD))/((xA-xB)*(xA-xC)*(xA-xD)))
	aB = 1.0D0*yB*(((x-xA)*(x-xC)*(x-xD))/((xB-xA)*(xB-xC)*(xB-xD)))
	aC = 1.0D0*yC*(((x-xA)*(x-xB)*(x-xD))/((xC-xA)*(xC-xB)*(xC-xD)))
	aD = 1.0D0*yD*(((x-xA)*(x-xB)*(x-xC))/((xD-xA)*(xD-xB)*(xD-xC)))
		
	y = 1.0D0*(aA+aB+aC+aD)
  End Function CubicInterpolationCalc
  
  
  
  Function CalcResidualSquareSum(inputPoints,polyCoefficients) RESULT (rss) 
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: order
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
	Real(kind=DoubleReal) :: rss, x, y	
	order = size(polyCoefficients)-1	
	rss = 0.0D0
	do i=1,size(inputPoints,1)
	  x = 1.0D0*inputPoints(i,1)
	  y = calcPolynomial(polyCoefficients,x,0)
	  rss = rss + (y-inputPoints(i,2))**2  
	enddo  
  End Function CalcResidualSquareSum   
  
  
  
  Function PointInterpolation(inputPoints,x,subsetSizeIn,inputSetStartIn,&
    inputSetLengthIn,orderDataIn) RESULT (y) 
!Lagrange Formula Point Interpolation, returns y
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,n,k
	Integer(kind=StandardInteger) :: order, inputSetStart, inputSetEnd
    Real(kind=DoubleReal) :: x,y,dy,numerator,denominator,numeratorPart
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPointsSubset
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPointsInterp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
!Optional vars
	Integer(kind=StandardInteger), optional :: subsetSizeIn
	Integer(kind=StandardInteger), optional :: inputSetStartIn
	Integer(kind=StandardInteger), optional :: inputSetLengthIn
	Character(3), optional :: orderDataIn
!Sort vars
	Character(3) :: orderData
	Real(kind=DoubleReal) :: sortXA, sortXB, sortYA, sortYB
	Logical :: sorting
!Subset vars	
	Integer(kind=StandardInteger) :: subsetSize, dataSize, xPos, posOffset
	Real(kind=DoubleReal) :: xLower, xUpper
!Set variables	
	dataSize = size(inputPoints,1)
	If(present(inputSetStartIn))Then
	  inputSetStart = inputSetStartIn
	Else
	  inputSetStart = 1
	End If
	If(present(inputSetLengthIn))Then
	  inputSetEnd = inputSetStartIn+inputSetLengthIn-1
	  dataSize = inputSetLengthIn
	Else
	  inputSetEnd = dataSize
	End If	
!Order data - don't do by default as it takes time
    If(present(orderDataIn))Then
	  orderData = orderDataIn
	  If(StrToUpper(orderData(1:3)).eq."ASC")Then
	    sorting = .true.
		Do While(sorting.eqv..true.)
		  sorting = .false.
	      !Do i=1,(dataSize-1)
		  Do i=inputSetStart,(inputSetEnd-1)
		    If(inputPoints(i,1).gt.inputPoints(i+1,1))Then			  
			  sorting = .true.
			  sortXA = inputPoints(i,1)
			  sortXB = inputPoints(i+1,1)
			  sortYA = inputPoints(i,2)
			  sortYB = inputPoints(i+1,2)
			  inputPoints(i,1) = sortXB
			  inputPoints(i+1,1) = sortXA
			  inputPoints(i,2) = sortYB
			  inputPoints(i+1,2) = sortYA
			End If
		  End Do
		End Do
	  End If
	  If(StrToUpper(orderData(1:3)).eq."DES")Then
	    sorting = .true.
		Do While(sorting.eqv..true.)
		  sorting = .false.
	      !Do i=1,(dataSize-1)
		  Do i=inputSetStart,(inputSetEnd-1)
		    If(inputPoints(i,1).lt.inputPoints(i+1,1))Then			  
			  sorting = .true.
			  sortXA = inputPoints(i,1)
			  sortXB = inputPoints(i+1,1)
			  sortYA = inputPoints(i,2)
			  sortYB = inputPoints(i+1,2)
			  inputPoints(i,1) = sortXB
			  inputPoints(i+1,1) = sortXA
			  inputPoints(i,2) = sortYB
			  inputPoints(i+1,2) = sortYA
			End If
		  End Do
		End Do
	  End If
	End If
!If interpolating a subset of inputPoints, make the subset
	If(present(subsetSizeIn))Then
	  subsetSize = subsetSizeIn
!Check size
	  If(subsetSize.lt.2)Then
	    subsetSize = 2
	  End If
	  If(subsetSize.gt.dataSize)Then
	    subsetSize = dataSize
	  End If	  
!reduce set of data points
      xLower = inputPoints(inputSetStart,1)
	  xUpper = inputPoints(inputSetEnd,1)
	  If(x.lt.xLower)Then
	    xPos = 1
	  Elseif(x.gt.dataSize)Then
	    xPos = dataSize
	  Else
!Estimate position
		xPos = Floor(((x - xLower) / (xUpper - xLower)) * 1.0D0 * dataSize) + inputSetStart
		xLower = inputPoints(xPos,1)
		xUpper = inputPoints(xPos+1,1)
		If(x.ge.xLower.and.x.le.xUpper)Then
!Position found  
		Else
!Find position	
		  Do i=1,dataSize
		    xLower = inputPoints(xPos+i,1)
		    xUpper = inputPoints(xPos+i+1,1)
		    If((xPos+i+1).eq.inputSetEnd.or.(x.ge.xLower.and.x.le.xUpper))Then
			  xPos = xPos + i
			  Exit
			End If
		    xLower = inputPoints(xPos-i,1)
		    xUpper = inputPoints(xPos-i+1,1)
		    If((xPos-i).eq.inputSetStart.or.(x.ge.xLower.and.x.le.xUpper))Then
			  xPos = xPos-i
			  Exit
			End If		  
		  End Do
!Set offset
          posOffset = -1*floor(subsetSize/2.0D0)
		  If((xPos+posOffset).lt.inputSetStart)Then
		    posOffset = 0  
		  End If
		  If((xPos-posOffset).gt.inputSetEnd)Then
		    posOffset = -1*subsetSize  
		  End If
!Allocate subset array
          Allocate(inputPointsSubset(1:subsetSize,1:2))
		  Do i=1,subsetSize
		    inputPointsSubset(i,1) = inputPoints(xPos+posOffset+i-1,1)
		    inputPointsSubset(i,2) = inputPoints(xPos+posOffset+i-1,2)
		  End Do
          !Deallocate(inputPoints)
          Allocate(inputPointsInterp(1:subsetSize,1:2))
		  Do i=1,subsetSize
		    inputPointsInterp(i,1) = inputPointsSubset(i,1)
		    inputPointsInterp(i,2) = inputPointsSubset(i,2)
		  End Do
		  dataSize = subsetSize		  
		End If
	  End If
	End If
!set variables
    order = dataSize-1
!Allocate coefficients array	
	Allocate(coefficients(0:order))
!make coefficients for Lagrange's formula	
	Do n=0,order
	  numerator = 1.0D0
	  denominator = 1.0D0 
	  Do k=0,order
	    If(k.ne.n)Then
	      numerator=numerator*(x-inputPointsInterp(k+1,1))
		  denominator=denominator*(inputPointsInterp(n+1,1)-inputPointsInterp(k+1,1))
		End If  
	  End Do
	  coefficients(n)=1.0D0*(numerator/denominator)
    End Do  
!Calculate y
    y = 0.0D0
	Do n=0,order
	  y=y+inputPointsInterp(n+1,2)*coefficients(n)
	End Do
!make coefficients for Lagrange's formula for y'(x)
	Do n=0,order
	  numerator = 0.0D0
	  denominator = 1.0D0 
	  Do k=0,order
	    If(k.ne.n)Then
		  denominator=denominator*(inputPointsInterp(n+1,1)-inputPointsInterp(k+1,1))
		  numeratorPart = 1.0D0
		  Do i=0,order
		    If(i.ne.n.and.i.ne.k)Then
			  numeratorPart=numeratorPart*(x-inputPointsInterp(i+1,1)) 
			End If
		  End Do
		  numerator=numerator+numeratorPart
		End If  
	  End Do
	  coefficients(n)=1.0D0*(numerator/denominator)
    End Do  
!Calculate dy
    dy = 0.0D0
  	Do n=0,order
	  dy=dy+inputPointsInterp(n+1,2)*coefficients(n)
	End Do
	
  End Function PointInterpolation  
  
  
  
  
  
  
  Function PointInterpolationArr(inputPoints,x,subsetSizeIn,inputSetStartIn,&
    inputSetLengthIn,verboseIn,orderDataIn) RESULT (yArray) 
!Lagrange Formula Point Interpolation, returns y
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,n,k
	Integer(kind=StandardInteger) :: order, inputSetStart, inputSetEnd
	Integer(kind=StandardInteger) :: verbose
    Real(kind=DoubleReal) :: x,y,dy,numerator,denominator,numeratorPart
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPointsSubset
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPointsInterp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
!Optional vars
	Integer(kind=StandardInteger), optional :: subsetSizeIn
	Integer(kind=StandardInteger), optional :: inputSetStartIn
	Integer(kind=StandardInteger), optional :: inputSetLengthIn
	Character(3), optional :: orderDataIn
	Character(1), optional :: verboseIn
!Sort vars
	Character(3) :: orderData
	Real(kind=DoubleReal) :: sortXA, sortXB, sortYA, sortYB
	Logical :: sorting
!Subset vars	
	Integer(kind=StandardInteger) :: subsetSize, dataSize, xPos, posOffset
	Real(kind=DoubleReal) :: xLower, xUpper
!Set variables	
	dataSize = size(inputPoints,1)
	If(present(inputSetStartIn))Then
	  inputSetStart = inputSetStartIn
	Else
	  inputSetStart = 1
	End If
	If(present(inputSetLengthIn))Then
	  inputSetEnd = inputSetStartIn+inputSetLengthIn-1
	  dataSize = inputSetLengthIn
	Else
	  inputSetEnd = dataSize
	End If	
	verbose = 0
	If(present(verboseIn))Then
	  If(verboseIn(1:1).eq."y".or.verboseIn(1:1).eq."Y")Then
	    verbose = 1
	  End If
	End If
!Order data - don't do by default as it takes time
    If(verbose.eq.1)Then
      print *,"Interp start"
	  print *,"Input x:"
	  print *,"   ",x
	  print *,"Input points:"
	  Do i=inputSetStart,inputSetEnd
	    print *,"   ",inputPoints(i,1),inputPoints(i,2)
	  End Do
	End If
    If(present(orderDataIn))Then
	  orderData = orderDataIn
	  If(StrToUpper(orderData(1:3)).eq."ASC")Then
	    sorting = .true.
		Do While(sorting.eqv..true.)
		  sorting = .false.
	      !Do i=1,(dataSize-1)
		  Do i=inputSetStart,(inputSetEnd-1)
		    If(inputPoints(i,1).gt.inputPoints(i+1,1))Then			  
			  sorting = .true.
			  sortXA = inputPoints(i,1)
			  sortXB = inputPoints(i+1,1)
			  sortYA = inputPoints(i,2)
			  sortYB = inputPoints(i+1,2)
			  inputPoints(i,1) = sortXB
			  inputPoints(i+1,1) = sortXA
			  inputPoints(i,2) = sortYB
			  inputPoints(i+1,2) = sortYA
			End If
		  End Do
		End Do
	  End If
	  If(StrToUpper(orderData(1:3)).eq."DES")Then
	    sorting = .true.
		Do While(sorting.eqv..true.)
		  sorting = .false.
	      !Do i=1,(dataSize-1)
		  Do i=inputSetStart,(inputSetEnd-1)
		    If(inputPoints(i,1).lt.inputPoints(i+1,1))Then			  
			  sorting = .true.
			  sortXA = inputPoints(i,1)
			  sortXB = inputPoints(i+1,1)
			  sortYA = inputPoints(i,2)
			  sortYB = inputPoints(i+1,2)
			  inputPoints(i,1) = sortXB
			  inputPoints(i+1,1) = sortXA
			  inputPoints(i,2) = sortYB
			  inputPoints(i+1,2) = sortYA
			End If
		  End Do
		End Do
	  End If
	End If
!If interpolating a subset of inputPoints, make the subset
	If(present(subsetSizeIn))Then
	  subsetSize = subsetSizeIn
!Check size
	  If(subsetSize.lt.2)Then
	    subsetSize = 2
	  End If
	  If(subsetSize.gt.dataSize)Then
	    subsetSize = dataSize
	  End If
	Else
      subsetSize = dataSize	
    End If
	If(subsetSize.lt.dataSize)Then
!reduce set of data points
      xLower = inputPoints(inputSetStart,1)
	  xUpper = inputPoints(inputSetEnd,1)
	  posOffset = 0
	  If(x.lt.xLower)Then
	    xPos = 1
		posOffset = 0
	  Elseif(x.gt.dataSize)Then
	    xPos = dataSize
		posOffset = subsetSize
	  Else
!Estimate position
		xPos = Floor(((x - xLower) / (xUpper - xLower)) * 1.0D0 * dataSize) + inputSetStart
		xLower = inputPoints(xPos,1)
		xUpper = inputPoints(xPos+1,1)
		If(x.ge.xLower.and.x.le.xUpper)Then
!Position found  
	      If(verbose.eq.1)Then
            print *,"Position found"
		    print *,"   ",xPos,dataSize,posOffset,subsetSize,inputSetEnd
		  End If	
		Else
!Find position	
		  Do i=1,dataSize
		    xLower = inputPoints(xPos+i,1)
		    xUpper = inputPoints(xPos+i+1,1)
		    If((xPos+i+1).eq.inputSetEnd.or.(x.ge.xLower.and.x.le.xUpper))Then
			  xPos = xPos + i
			  Exit
			End If
		    xLower = inputPoints(xPos-i,1)
		    xUpper = inputPoints(xPos-i+1,1)
		    If((xPos-i).eq.inputSetStart.or.(x.ge.xLower.and.x.le.xUpper))Then
			  xPos = xPos-i
			  Exit
			End If		  
		  End Do
!Set offset
	      If(verbose.eq.1)Then
		    print *,"Set Offset"
		  End If
          posOffset = -1*floor(subsetSize/2.0D0)
		End If
	  End If
!adjust position with offset
	  If(verbose.eq.1)Then
        print *,"Pre-Adjust offset"
		print *,"   ",xPos,dataSize,posOffset,subsetSize,inputSetEnd
	  End If
	  If((xPos+posOffset).lt.inputSetStart)Then
		posOffset = 0  
	  End If
	  If((xPos-posOffset).gt.inputSetEnd)Then
	    posOffset = -1*subsetSize  
	  End If
	  If((xPos+subsetSize).gt.inputSetEnd)Then
	    posOffset = -1*subsetSize  
	  End If
	  If(verbose.eq.1)Then
        print *,"Adjust offset"
		print *,"   ",xPos,dataSize,posOffset,subsetSize,inputSetEnd
	  End If
!Allocate subset array
      Allocate(inputPointsSubset(1:subsetSize,1:2))	  
	  Do i=1,subsetSize
		inputPointsSubset(i,1) = inputPoints(xPos+posOffset+i-1,1)
		inputPointsSubset(i,2) = inputPoints(xPos+posOffset+i-1,2)
	    If(verbose.eq.1)Then
		  print *,inputPoints(xPos+posOffset+i-1,1),inputPoints(xPos+posOffset+i-1,2)
		End If
	  End Do
      Allocate(inputPointsInterp(1:subsetSize,1:2))
	  Do i=1,subsetSize
		inputPointsInterp(i,1) = inputPointsSubset(i,1)
		inputPointsInterp(i,2) = inputPointsSubset(i,2)
	  End Do
	  If(verbose.eq.1)Then
		print *,"Subset of points:"
		Do i=1,subsetSize
	      print *,"   ",inputPointsInterp(i,1),inputPointsInterp(i,2)
		End Do  
	  End If
	  dataSize = subsetSize	
	Else  
	  Allocate(inputPointsInterp(1:subsetSize,1:2))
	  Do i=1,subsetSize
		inputPointsInterp(i,1) = inputPoints(inputSetStart+i-1,1)
		inputPointsInterp(i,2) = inputPoints(inputSetStart+i-1,2)
	  End Do 
	End If
!set variables
    order = dataSize-1
!Allocate coefficients array	
	Allocate(coefficients(0:order))
!make coefficients for Lagrange's formula	
	Do n=0,order
	  numerator = 1.0D0
	  denominator = 1.0D0 
	  Do k=0,order
	    If(k.ne.n)Then
	      numerator=numerator*(x-inputPointsInterp(k+1,1))
		  denominator=denominator*(inputPointsInterp(n+1,1)-inputPointsInterp(k+1,1))
		End If  
	  End Do
	  coefficients(n)=1.0D0*(numerator/denominator)
    End Do  
!Calculate y
    y = 0.0D0
	Do n=0,order
	  y=y+inputPointsInterp(n+1,2)*coefficients(n)
	End Do
!make coefficients for Lagrange's formula for y'(x)
	Do n=0,order
	  numerator = 0.0D0
	  denominator = 1.0D0 
	  Do k=0,order
	    If(k.ne.n)Then
		  denominator=denominator*(inputPointsInterp(n+1,1)-inputPointsInterp(k+1,1))
		  numeratorPart = 1.0D0
		  Do i=0,order
		    If(i.ne.n.and.i.ne.k)Then
			  numeratorPart=numeratorPart*(x-inputPointsInterp(i+1,1)) 
			End If
		  End Do
		  numerator=numerator+numeratorPart
		End If  
	  End Do
	  coefficients(n)=1.0D0*(numerator/denominator)
    End Do  
!Calculate dy
    dy = 0.0D0
  	Do n=0,order
	  dy=dy+inputPointsInterp(n+1,2)*coefficients(n)
	End Do
!Store results in array
	Allocate(yArray(1:2))
	yArray(1) = y
	yArray(2) = dy
	If(verbose.eq.1)Then
	  print *,"Interp end"
	End If
  End Function PointInterpolationArr
  
  
  Function PointInterpolationFull(inputPoints,x) RESULT (output) 
!Lagrange Formula Point Interpolation, returns array of y, y'
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: n,k,i
	Integer(kind=StandardInteger) :: order
    Real(kind=DoubleReal) :: x,y,dy,numerator,denominator,numeratorPart
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: output
!set variables
    order = size(inputPoints,1)-1
!Allocate coefficients array	
	Allocate(coefficients(0:order))
	Allocate(output(1:2))
!make coefficients for Lagrange's formula for y(x)
	Do n=0,order
	  numerator = 1.0D0
	  denominator = 1.0D0 
	  Do k=0,order
	    If(k.ne.n)Then
	      numerator=numerator*(x-inputPoints(k+1,1))
		  denominator=denominator*(inputPoints(n+1,1)-inputPoints(k+1,1))
		End If  
	  End Do
	  coefficients(n)=1.0D0*(numerator/denominator)
    End Do  
!Calculate y
    y = 0.0D0
	Do n=0,order
	  y=y+inputPoints(n+1,2)*coefficients(n)
	End Do
!make coefficients for Lagrange's formula for y'(x)
	Do n=0,order
	  numerator = 0.0D0
	  denominator = 1.0D0 
	  Do k=0,order
	    If(k.ne.n)Then
		  denominator=denominator*(inputPoints(n+1,1)-inputPoints(k+1,1))
		  numeratorPart = 1.0D0
		  Do i=0,order
		    If(i.ne.n.and.i.ne.k)Then
			  numeratorPart=numeratorPart*(x-inputPoints(i+1,1)) 
			End If
		  End Do
		  numerator=numerator+numeratorPart
		End If  
	  End Do
	  coefficients(n)=1.0D0*(numerator/denominator)
    End Do  
!Calculate dy
    dy = 0.0D0
  	Do n=0,order
	  dy=dy+inputPoints(n+1,2)*coefficients(n)
	End Do
!return data as array
	output(1) = y
	output(2) = dy
  End Function PointInterpolationFull 
  
  
  
  
  
  
  Function CalcSecondDerivativeMinimum(inputPoints, orderIn) RESULT (secondDerivativeMinimum)  
!Takes input points, fits polynomial, finds minimum f'(x) = 0, then gives f''(x) at that point  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,order
!Optional
    Integer(kind=StandardInteger), optional :: orderIn
	Real(kind=DoubleReal) :: upperBound, lowerBound, boundDifference, boundMean
	Real(kind=DoubleReal) :: minimumValue, secondDerivative
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficientsDerivative
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients2ndDerivative
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension(1:2) :: secondDerivativeMinimum
!Set optional values
    If(present(orderIn))Then
	  order = orderIn
	Else
	  order = 5
	End If
!Allocate arrays
    Allocate(coefficients(0:order))
    Allocate(coefficientsDerivative(0:(order-1)))
    Allocate(coefficients2ndDerivative(0:(order-2)))
    coefficients = PolyFit(inputPoints,order)
!Get derivative
    coefficientsDerivative = PolynomialCoefficientDerivative(coefficients)
	!print *,coefficients(0),coefficients(1),&
	!coefficients(2),coefficients(3),coefficients(4),coefficients(5)
	!print *,coefficientsDerivative(0),coefficientsDerivative(1),&
	!coefficientsDerivative(2),coefficientsDerivative(3),coefficientsDerivative(4)
!Solve for f'(x) = 0
    Do i=1,size(inputPoints,1)
	  If(i.eq.1)Then
	    lowerBound = inputPoints(i,1)
	    upperBound = inputPoints(i,1)
	  Else	
		If(inputPoints(i,1).lt.lowerBound)Then
		  lowerBound = inputPoints(i,1)
		End If
		If(inputPoints(i,1).gt.upperBound)Then
		  upperBound = inputPoints(i,1)
		End If
	  End If
    End Do
	boundDifference = 1.0D0*upperBound - 1.0D0*lowerBound
	boundMean = 0.5D0 * boundDifference + 1.0D0*lowerBound
	!lowerBound = 1.0D0*boundMean - 2.0D0*boundDifference
	!upperBound = 1.0D0*boundMean + 2.0D0*boundDifference
	minimumValue=SolvePolynomial(coefficientsDerivative,lowerBound,upperBound)
	coefficients2ndDerivative = PolynomialCoefficientDerivative(coefficientsDerivative)
	secondDerivative = 0.0D0
	Do i=0,(order-2)
	  secondDerivative = secondDerivative + coefficients2ndDerivative(i)*minimumValue**i
	End Do
	secondDerivativeMinimum(1)=secondDerivative
	secondDerivativeMinimum(2)=minimumValue
  End Function CalcSecondDerivativeMinimum  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
! Spline Functions
!------------------------------------------------------------------------!   

  Function Spline(inputPoints,nIn,subsetStartIn,subsetLengthIn,orderDataIn) RESULT (functionXY) 
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,m
	Real(kind=DoubleReal) :: x, y, dx, xA, xB
	Integer(kind=StandardInteger) :: posOffset, interpPointsCount
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splinePoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: interpPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: functionXY
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray, coefficients, cIn
    Real(kind=DoubleReal) :: sortXA, sortXB, sortYA, sortYB
	Logical :: sorting
!optional variables
    Integer(kind=StandardInteger), optional :: subsetStartIn, subsetLengthIn, orderDataIn, nIn
    Integer(kind=StandardInteger) :: subsetStart, subsetLength, subsetEnd, orderData, n
!sort out optional
    If(Present(nIn))Then
      n = nIn
	Else
	  n = 1001
	End If
    If(Present(subsetStartIn))Then
      subsetStart = subsetStartIn
	Else
	  subsetStart = 1
	End If
    If(Present(subsetLengthIn))Then
      subsetLength = subsetLengthIn
	Else
	  subsetLength = size(inputPoints,1)
	End If
	subsetEnd = subsetStart+subsetLength-1
	If(Present(orderDataIn))Then
	  orderData = orderDataIn
	Else  
	  orderData = 0 
	End If
!n point interpolation, minimum interpPoints points needed
    interpPointsCount = 3
	If(subsetLength.ge.interpPointsCount)Then
!Allocate arrays
      Allocate(splinePoints(1:subsetLength,1:5))
	  Allocate(interpPoints(1:interpPointsCount,1:2))
      Allocate(coefficients(0:5))
      Allocate(functionXY(1:n,1:4))
!Order the subset of data
      If(orderData.eq.1)Then
	    sorting = .true.
	    Do While(sorting.eqv..true.)
	      sorting = .false.
		  Do i=subsetStart,(subsetEnd-1)
		    If(inputPoints(i,1).gt.inputPoints(i+1,1))Then			  
			  sorting = .true.
			  Call swapMatrixRows(inputPoints,i,i+1)
		    End If
			
          End Do		  
	    End Do
	  End If	
!Find dy/dx and d2y/dx2 at each point
      j=0
      Do i=subsetStart,subsetEnd
!Interp points
        x = inputPoints(i,1)
		y = inputPoints(i,2) 
        posOffset = 0	    
		If((i+interpPointsCount).gt.subsetEnd)Then
		  posOffset = ((i+interpPointsCount)-subsetEnd)-1
		End If
		Do k=1,interpPointsCount
		  interpPoints(k,1) = inputPoints(i-posOffset+(k-1),1)
		  interpPoints(k,2) = inputPoints(i-posOffset+(k-1),2)
		End Do
		yArray = MatrixInterpolation(interpPoints,x)
!store data
	    j=j+1
		splinePoints(j,1) = inputPoints(i,1)	!x
		splinePoints(j,2) = inputPoints(i,2)    !y
		splinePoints(j,3) = yArray(3)           !dy 
		splinePoints(j,4) = yArray(4)           !d2y 
		splinePoints(j,5) = yArray(5)           !d3y 
		!print *,j,subsetEnd,splinePoints(j,1),splinePoints(j,2),splinePoints(j,3),splinePoints(j,4)
	  End Do
	  x = splinePoints(1,1)
	  dx = 1.0D0*((splinePoints(size(splinePoints,1),1)-splinePoints(1,1))/(n-1))
	  Do i=1,n
	    functionXY(i,1) = x+(i-1)*dx
	  End Do
!spline between points
      m=1
      Do i=1,(size(splinePoints,1)-1)
	  !Do i=1,5
	    xA = splinePoints(i,1)
	    xB = splinePoints(i+1,1)
        coefficients = SplineTwoPoints(splinePoints,i,i+1)
	    Do j=1,n		
		  x = functionXY(j,1)
		  If(x.ge.xA.and.x.le.xB)Then
		    cIn = coefficients
		    functionXY(j,2) = CalcPolynomial(cIn, x, 0)
		    cIn = coefficients
		    functionXY(j,3) = CalcPolynomial(cIn, x, 1)
		    cIn = coefficients
		    functionXY(j,4) = CalcPolynomial(cIn, x, 2)
          End If 
	    End Do
	  End Do
    Else !If fewer than five points
	!Skip
	End If
	

	

  End Function Spline  
  
  
  
  
  
  Function SplineTwoPoints(splinePoints,pointA,pointB) RESULT (coefficients) 
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: expon, row, column
	Integer(kind=StandardInteger) :: point, pointA, pointB
	Real(kind=DoubleReal) :: x
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splinePoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix, xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix, cMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: multCoefficients 
!sorting
    Integer(kind=StandardInteger) :: sortingExpon
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: sortingValues	
	Logical :: sorting
	
!Allocate arrays
    Allocate(coefficients(0:5))
    Allocate(multCoefficients(1:2,1:6))
	Allocate(xMatrix(1:6,1:6))
	Allocate(xMatrixInverse(1:6,1:6))
	Allocate(yMatrix(1:6))
	Allocate(cMatrix(1:6))
	Allocate(sortingValues(1:6))
!Set values
	!Row 1
	x = splinePoints(pointA,1)
	Do i=1,6
	  xMatrix(1,i) = 1.0D0*x**(i-1)
	End Do
	!Row 2
	x = splinePoints(pointB,1)
	Do i=1,6
	  xMatrix(2,i) = 1.0D0*x**(i-1)
	End Do
	!Row 3
	x = splinePoints(pointA,1)
	xMatrix(3,1) = 0.0D0
	Do i=2,6
	  xMatrix(3,i) = 1.0D0*(i-1)*x**(i-2)
	End Do
	!Row 4
	x = splinePoints(pointB,1)
	xMatrix(4,1) = 0.0D0
	Do i=2,6
	  xMatrix(4,i) = 1.0D0*(i-1)*x**(i-2)
	End Do
	!Row 5
	x = splinePoints(pointA,1)
	xMatrix(5,1) = 0.0D0
	xMatrix(5,2) = 0.0D0
	Do i=3,6
	  xMatrix(5,i) = 1.0D0*(i-2)*(i-1)*x**(i-3)
	End Do
	!Row 6
	x = splinePoints(pointB,1)
	xMatrix(6,1) = 0.0D0
	xMatrix(6,2) = 0.0D0
	Do i=3,6
	  xMatrix(6,i) = 1.0D0*(i-2)*(i-1)*x**(i-3)
	End Do
!Make yMatrix
    row = 0
	Do j=1,3
      Do i=1,2
	    If(i.eq.1)Then
	      point = pointA
	    Else
	      point = pointB
	    End If
	    row = row + 1
		yMatrix(row) = splinePoints(point,j+1)
	  End Do
	End Do
!order the rows	
	Do i=1,6
	  sortingValues(i) = 0
	  Do j=1,6
	    sortingExpon= 7-j
		If(xMatrix(i,j).ne.0)Then
	      sortingValues(i) = sortingValues(i) + 2**sortingExpon
		End If
	  End Do 
	End Do  
	sorting = .true.
	Do While(sorting.eqv..true.)
	  sorting = .false.
	  Do i=1,5
	    If(sortingValues(i).lt.sortingValues(i+1))Then
	      sorting = .true.
		  Call swapMatrixRows1D(sortingValues,i,i+1)
		  Call swapMatrixRows2D(xMatrix,i,i+1)
		  Call swapMatrixRows1D(yMatrix,i,i+1)
	    End If
	  End Do	
	End Do
!Solve equation XC=Y
    xMatrixInverse = InvertSquareMatrix(xMatrix)
	cMatrix = matmul(xMatrixInverse,yMatrix)
!Copy to coefficients matrix
    Do i=0,5	
	  coefficients(i) = cMatrix(i+1)
	End Do
  End Function SplineTwoPoints


  
   
  Function MatrixInterpolation(inputPoints,x) RESULT (yArray)  
!yArray(1) = x, yArray(2) = y, yArray(3) = dy, yArray(4) = d2y, yArray(5) = d3y
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: points
	Real(kind=DoubleReal) :: x
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix, xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix, cMatrix, yArray
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients, cIn
	Logical :: sorting
	Real(kind=DoubleReal) :: sortXA, sortXB, sortYA, sortYB
!Set values	
	points = size(inputPoints,1)
!If more than 2 sets of x,y points	
    If(points.ge.2)Then
	  inputPoints = 1.0D0*inputPoints
!Order the input points, so any "zero" points are at the bottom	
	  sorting = .true.
	  Do while(sorting.eqv..true.)
	    sorting = .false.
		Do i=1,(points-1)
		  If(abs(inputPoints(i,1)).lt.abs(inputPoints(i+1,1)))Then
		    sorting = .true.
			Call swapMatrixRows(inputPoints,i,i+1)
		  End If
		End Do
	  End Do
!Allocate matrices
	  Allocate(xMatrix(1:points,1:points))
	  Allocate(yMatrix(1:points))
	  Allocate(cMatrix(1:points))
	  Allocate(coefficients(0:(points-1)))
!Fill array
      Do i=1,points
        Do j=1,points
          xMatrix(i,j) = 1.0D0*inputPoints(i,1)**(j-1)
        End Do	
        yMatrix(i) = 1.0D0*inputPoints(i,2)
      End Do
!Solve equation XC=Y
      xMatrixInverse = InvertSquareMatrix(xMatrix)
	  cMatrix = matmul(xMatrixInverse,yMatrix)
!Copy to coefficients matrix
      Do i=0,(points-1)	
	    coefficients(i) = cMatrix(i+1)
	  End Do
	  Allocate(yArray(1:5))
	  yArray = 0.0D0
!Calculate values
	  yArray(1) = x
	  cIn = coefficients
	  yArray(2) = CalcPolynomial(cIn, x, 0)
	  cIn = coefficients
	  yArray(3) = CalcPolynomial(cIn, x, 1)
	  If(points.gt.2)Then	  
	    yArray(4) = CalcPolynomial(cIn, x, 2)
	  End If
	  If(points.gt.3)Then	  
	    yArray(5) = CalcPolynomial(cIn, x, 3)
	  End If
	Else  	  
	  Allocate(yArray(1:1))
	  yArray(1) = x
	  yArray(2) = inputPoints(1,2)
    End If
  End Function MatrixInterpolation  
  
  
  
!------------------------------------------------------------------------!
! Matrix Functions
!------------------------------------------------------------------------! 
  
  Function InvertSquareMatrix(xMatrix) RESULT (xMatrixInverse)  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,row,col,rowb
	Integer(kind=StandardInteger) :: matrixSize
	Real(kind=DoubleReal) :: xA, xB
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixWorking
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: xMatrixRow
!if a square matrix
	if(size(xMatrix,1).eq.size(xMatrix,2))then
	  matrixSize = size(xMatrix,1)
!Allocate arrays
	  Allocate(xMatrixInverse(1:matrixSize,1:matrixSize))
	  Allocate(xMatrixWorking(1:matrixSize,1:2*matrixSize))	 
	  Allocate(xMatrixRow(1:2*matrixSize)) 
!Initialise arrays	
      xMatrixWorking = 0.0D0
      xMatrixInverse = 0.0D0  
      xMatrixRow = 0.0D0  
!Fill working array
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
	    enddo
	  enddo
      do row=1,matrixSize
	    do col=1,matrixSize
		  if(row.eq.col)then
	        xMatrixWorking(row,col+matrixSize) = 1.0D0
		  endif
	    enddo
	  enddo	
!make lower triangle of zeros	  
	  Do row=1,matrixSize-1
	    Do rowb=row+1,matrixSize	 
		  If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
		    Do col=1,(2*matrixSize) !loop over all columns
		      xMatrixRow(col) = 1.0D0*&
		      ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
			  xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
		    End Do
!replace row values
		    do col=1,(2*matrixSize) !loop over all columns
		      xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
		    End Do
	      End If
		End Do
!force zeros in the lower triangle
        Do rowb=row+1,matrixSize
	      xMatrixWorking(rowb,row) = 0.0D0
	    End Do
	  End Do
!re-force zeros in the lower triangle
      Do row=1,matrixSize
	    Do col=1,matrixSize
		  if(row.gt.col)then
		    xMatrixWorking(row,col) = 0.0D0
		  endif
		End Do
	  End Do	
!make upper triangle of zeros	
	  do row=matrixSize,2,-1
	    do rowb=row-1,1,-1	  
		  If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
		    Do col=1,(2*matrixSize) !loop over all columns
		      xMatrixRow(col) = 1.0D0*&
		      ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
			  xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
		    End Do
!replace row values
		    Do col=1,(2*matrixSize) !loop over all columns
		      xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
		    End Do
		  End If	
	    End Do
!force zeros in the upper triangle
        do rowb=row-1,1,-1
	      xMatrixWorking(rowb,row) = 0.0D0
	    enddo
	  enddo
!Divide rhs by diagonal on lhs and store in inverse
	  do row=1,matrixSize
	    do col=1,matrixSize
		  xMatrixInverse(row,col) = 1.0D0*&
		  xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
		enddo
	  enddo
	endif
  End Function InvertSquareMatrix  
  
    
    
  Function InvertSquareMatrixOld(xMatrixInput) RESULT (xMatrixInverse)  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,row,col,rowb
	Integer(kind=StandardInteger) :: matrixSize,optimiseSum,optimiseExponent
	Real(kind=DoubleReal) :: xA, xB
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInput
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixWorking
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: xMatrixRow
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: optimiseArray	
	Logical :: optimising	
!if a square matrix
	if(size(xMatrixInput,1).eq.size(xMatrixInput,2))then
	  matrixSize = size(xMatrixInput,1)
!Allocate arrays
	  Allocate(xMatrix(1:matrixSize,1:matrixSize))
	  Allocate(xMatrixInverse(1:matrixSize,1:matrixSize))
	  Allocate(xMatrixWorking(1:matrixSize,1:2*matrixSize))	 
	  Allocate(xMatrixRow(1:2*matrixSize)) 
	  Allocate(optimiseArray(1:matrixSize)) 
!Initialise arrays
      xMatrix = 0.0D0	
      xMatrixWorking = 0.0D0
      xMatrixInverse = 0.0D0  
      xMatrixRow = 0.0D0  
	  optimiseArray = 0
!copy out input array	  
	  xMatrix = 1.0D0*xMatrixInput
!Fill working array
      do row=1,matrixSize
	    do col=1,(2*matrixSize)
	      xMatrixWorking(row,col) = 0.0D0
	    enddo
	  enddo
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
	    enddo
	  enddo
!Fill inverse array
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixInverse(row,col) = 0.0D0
	    enddo
	  enddo
!optimise matrix order
      do row=1,matrixSize		!row at a time
	    optimiseSum = 0
	    do col=1,matrixSize 
	      if(xMatrixWorking(row,col).ne.0.0D0)then
	        optimiseExponent = matrixSize - col
		    optimiseSum = optimiseSum + 2**optimiseExponent
		  endif
	    enddo
        optimiseArray(row) = optimiseSum 
	  enddo	
	  optimising = .true.
      do while(optimising)
        optimising = .false.
        do row=1,(matrixSize-1)        !loop through rows
          if(optimiseArray(row).lt.optimiseArray(row+1))then
		    optimising = .true.
!reorder optimising array
            i = optimiseArray(row)
            j = optimiseArray(row+1)
		    optimiseArray(row) = j
		    optimiseArray(row+1) = i
!reorder xMatrixWorking
            do col=1,(2*matrixSize)   
		      xA = 1.0D0*xMatrixWorking(row,col)
              xB = 1.0D0*xMatrixWorking(row+1,col)
		      xMatrixWorking(row,col) = 1.0D0*xB
		      xMatrixWorking(row+1,col) = 1.0D0*xA
		    enddo
		  endif
	    enddo
      enddo	  
!Make identity in rhs matrix
      do row=1,matrixSize
	    do col=1,matrixSize
		  if(row.eq.col)then
	        xMatrixWorking(row,col+matrixSize) = 1.0D0
		  endif
	    enddo
	  enddo
!make lower triangle of zeros	  
	  do row=1,matrixSize-1
	    do rowb=row+1,matrixSize	  
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixRow(col) = 1.0D0*&
		    ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
			xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
		  enddo
!replace row values
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
		  enddo
	    enddo
!force zeros in the lower triangle
        do rowb=row+1,matrixSize
	      xMatrixWorking(rowb,row) = 0.0D0
	    enddo
	  enddo
!re-force zeros in the lower triangle
      do row=1,matrixSize
	    do col=1,matrixSize
		  if(row.gt.col)then
		    xMatrixWorking(row,col) = 0.0D0
		  endif
		enddo
	  enddo
!make upper triangle of zeros	
	  do row=matrixSize,2,-1
	    do rowb=row-1,1,-1	  
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixRow(col) = 1.0D0*&
		    ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
			xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
		  enddo
!replace row values
		  do col=1,(2*matrixSize) !loop over all columns
		    xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
		  enddo
	    enddo
!force zeros in the upper triangle
        do rowb=row-1,1,-1
	      xMatrixWorking(rowb,row) = 0.0D0
	    enddo
	  enddo
!Divide rhs by diagonal on lhs and store in inverse
	  do row=1,matrixSize
	    do col=1,matrixSize
		  xMatrixInverse(row,col) = 1.0D0*&
		  xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
		enddo
	  enddo
	  do row=1,matrixSize
		  print *,xMatrixInverse(row,1),xMatrixInverse(row,2),xMatrixInverse(row,3),&
		  xMatrixInverse(row,4),xMatrixInverse(row,5)
	  end do
	  print *,""
	endif
  End Function InvertSquareMatrixOld  	
	
	
 
  Function InvertSquareMatrixNew(xMatrix) RESULT (xMatrixInverse)  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,row,col,rowb
	Integer(kind=StandardInteger) :: matrixSize,optimiseSum,optimiseExponent
	Real(kind=DoubleReal) :: xA, xB
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixWorking
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: optimiseArray	
	Logical :: optimising	
!if a square matrix
	if(size(xMatrix,1).eq.size(xMatrix,2))then
	  matrixSize = size(xMatrix,1)
!Allocate arrays
	  Allocate(xMatrixInverse(1:matrixSize,1:matrixSize))
	  Allocate(xMatrixWorking(1:matrixSize,1:2*matrixSize))	 
	  Allocate(optimiseArray(1:matrixSize)) 
!Initialise arrays
      xMatrix = 1.0D0 * xMatrix	
      xMatrixWorking = 0.0D0
      xMatrixInverse = 0.0D0  
	  optimiseArray = 0
!Fill working array
      do row=1,matrixSize
	    do col=1,(2*matrixSize)
	      xMatrixWorking(row,col) = 0.0D0
	    enddo
	  enddo
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
	    enddo
	  enddo
!Fill inverse array
	  do row=1,matrixSize
	    do col=1,matrixSize
	      xMatrixInverse(row,col) = 0.0D0
	    enddo
	  enddo
!optimise matrix order
      do row=1,matrixSize		!row at a time
	    optimiseSum = 0
	    do col=1,matrixSize 
	      if(xMatrixWorking(row,col).ne.0.0D0)then
	        optimiseExponent = matrixSize - col
		    optimiseSum = optimiseSum + 2**optimiseExponent
		  endif
	    enddo
        optimiseArray(row) = optimiseSum 
	  enddo	
	  optimising = .true.
      do while(optimising)
        optimising = .false.
        do row=1,(matrixSize-1)        !loop through rows
          if(optimiseArray(row).lt.optimiseArray(row+1))then
		    optimising = .true.
!reorder optimising array
            i = optimiseArray(row)
            j = optimiseArray(row+1)
		    optimiseArray(row) = j
		    optimiseArray(row+1) = i
!reorder xMatrixWorking
            do col=1,(2*matrixSize)   
		      xA = 1.0D0*xMatrixWorking(row,col)
              xB = 1.0D0*xMatrixWorking(row+1,col)
		      xMatrixWorking(row,col) = 1.0D0*xB
		      xMatrixWorking(row+1,col) = 1.0D0*xA
		    enddo
		  endif
	    enddo
      enddo	  
!Make identity in rhs matrix
      do row=1,matrixSize
	    do col=1,matrixSize
		  if(row.eq.col)then
	        xMatrixWorking(row,col+matrixSize) = 1.0D0
		  endif
	    enddo
	  enddo
!print
      Do row=1,matrixSize
	    print "(E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6)",&
		xMatrixWorking(row,1),xMatrixWorking(row,2),xMatrixWorking(row,3),&
		xMatrixWorking(row,4),xMatrixWorking(row,5),xMatrixWorking(row,6),&
		xMatrixWorking(row,7),xMatrixWorking(row,8),xMatrixWorking(row,9),&
		xMatrixWorking(row,10)
	  End Do
	  
	  
	  
	endif
  End Function InvertSquareMatrixNew  	
	
	
	
	
  Function SolveLinear(xMatrix,yMatrix) RESULT (cMatrix)  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: row,column
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix, cMatrix  
  
  
  End Function SolveLinear  
  
  
  
  Function Jacobian (inputMatix) RESULT (outputMatrix)  
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: matrixLength, interpPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputMatix
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputMatrix
!Currently 1D jacobian with perturbations stored in the 2D of the input matrix
!   r1(x),r1(x+e),r1(x+2e)
!   r2(x),r2(x+e),r2(x+2e)
!   r3(x),x3(x+e),r3(x+2e)
!   rn(x),xn(x+e),rn(x+2e)
  
  
  
  End Function Jacobian
  
  
  
!------------------------------------------------------------------------!
! Random Number Related Functions
!------------------------------------------------------------------------!    
  
  Subroutine SetRandomSeedArray()
    Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: clockReturn, gap, seedTemp, multiple
    Integer(kind=StandardInteger), Dimension(1:1000) :: randomSeed
	Integer(kind=StandardInteger), Dimension(0:9) :: randomIntegers
  !random number seed from cpu
    randomIntegers(0) = 25441
    randomIntegers(1) = 37261
    randomIntegers(2) = 1622261
    randomIntegers(3) = 162982	
    randomIntegers(4) = 72635	
    randomIntegers(5) = 9927151
    randomIntegers(6) = 91
    randomIntegers(7) = 6452
    randomIntegers(8) = 448327
    randomIntegers(9) = 9253411
	j = 0
	Call SYSTEM_CLOCK(clockReturn)
	gap = clockReturn - 10*floor(1.0D0*clockReturn/10)
    do i=1,1000	  
	  k = j + gap
	  if(k.gt.9)then
	    k = k - 10
	  endif
	  seedTemp = i * (randomIntegers(j)-randomIntegers(k))
	  seedTemp = abs(seedTemp)
	  multiple = floor(1.0D0 * seedTemp / 1.0D0 * clockReturn)
	  seedTemp = seedTemp - multiple * clockReturn
	  randomSeed(i) = abs(clockReturn - seedTemp)
	  if(j.eq.9)then
	    j = 0
	  endif
	  j = j + 1
	enddo	
	Call RANDOM_SEED(put=randomSeed)
  End Subroutine SetRandomSeedArray
  
  
  
  Function RandomDataPoints(inputPoints,noPoints, tether) RESULT (outputPoints)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,noPoints,point,noInputPoints,tether
	Real(kind=DoubleReal) :: randNumber
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: selectedPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputPoints
	Logical :: pointOk
!set variables	
	noInputPoints = size(inputPoints,1) 
!Allocate arrays
    Allocate(selectedPoints(1:noPoints))
    Allocate(outputPoints(1:noPoints,1:size(inputPoints,2)))
!fill selected points array with 0s
	do i=1,noPoints
	  selectedPoints(i) = 0  
	enddo
	if(tether.eq.1)then
	  !tether start and end points
	  selectedPoints(1) = 1
	  selectedPoints(noPoints) = noInputPoints
	endif
	if(tether.eq.2)then
	  !tether start and end points
	  Call RANDOM_NUMBER(randNumber)
	  if(randNumber.lt.0.5)then
	    selectedPoints(1) = 1
	  else
	    selectedPoints(1) = 2
	  endif
	  Call RANDOM_NUMBER(randNumber)
	  if(randNumber.lt.0.5)then
	    selectedPoints(noPoints) = noInputPoints
	  else
	    selectedPoints(noPoints) = noInputPoints-1
	  endif
	endif
!pick random points
    i = 1
	do while(i.le.noPoints)
	  Call RANDOM_NUMBER(randNumber)  
	  point = ceiling(randNumber*noInputPoints)
	  pointOk = .true.
	  do j=1,noPoints
	    if(selectedPoints(j).eq.point)then
		  pointOk = .false.
		endif		
	  enddo
	  if(pointOk)then
		selectedPoints(i) = point
	    i = i + 1
	  endif
	enddo
!Transfer data to output array
    do i=1,noPoints
	  do j=1,size(inputPoints,2)
        outputPoints(i,j) = inputPoints(selectedPoints(i),j)
	  enddo
	enddo
  End Function RandomDataPoints  
  
  
  Function VaryPoint(x,T) RESULT (xV)  
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i, distributionDataPoints
	Integer(kind=StandardInteger) :: factor
	Real(kind=DoubleReal) :: x,T,xV
	Real(kind=DoubleReal) :: distX, distY, distXInterval
	Real(kind=DoubleReal) :: sigma, mu, piNum, normalisationFactor
	Real(kind=DoubleReal) :: randNumber, variation
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: distributionPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Logical :: addVariation
!Make distribution set
    distributionDataPoints = 21
	piNum = 3.1415926535898D0
	sigma = 0.28D0
	mu = 0.0D0	
!Allocate distribution data array
    Allocate(distributionPoints(1:distributionDataPoints,1:2))
!Set data points	
	distX = 0.0D0
	distXInterval = 1.0D0/(distributionDataPoints-1)
!Make array
    Do i=1,distributionDataPoints
	  distributionPoints(i,1) = distX
	  distributionPoints(i,2) = exp(-1*((distX-mu)**2)/(2*sigma**2))
	  distX = distX + distXInterval
	End Do	
!Whether to add or subtract variation
    Call RANDOM_NUMBER(randNumber)
	addVariation = .false.
	If((1.0D0*randNumber).ge.0.5D0)Then
	  addVariation = .true.
	End If
!Generate random variation
    Call RANDOM_NUMBER(randNumber)
	yArray = PointInterpolationArr(distributionPoints,randNumber,5)
	variation = yArray(1)
!New point
    If(addVariation.eqv..true.)Then
      xV = x*(1+variation*T)
    Else
	  xV = x*(1-variation*T)
	End If  
  End Function VaryPoint  
  
  
  Function VaryPoints(inputArray,T,pointIn) RESULT (outputArray)  
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i, point
	Real(kind=DoubleReal) :: y,T
	Integer(kind=StandardInteger) :: factor
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: inputArray
	Real(kind=DoubleReal), Dimension(:,:), Allocatable :: outputArray
	Integer(kind=StandardInteger), optional :: pointIn
!Allocate output array	
	Allocate(outputArray(1:size(inputArray,1),1:2))
!Optional in
    point = 0
    If(present(pointIn))Then
	  point = pointIn
	End If
!Loop over points
    Do i=1,size(inputArray,1)	  
	  If(point.eq.0)Then
	    y = inputArray(i,2)
	    outputArray(i,1) = inputArray(i,1)
	    outputArray(i,2) = VaryPoint(y,T)
	  Else
	    If(i.eq.point)Then
		  y = inputArray(i,2)
	      outputArray(i,1) = inputArray(i,1)
	      outputArray(i,2) = VaryPoint(y,T)
		Else
	      outputArray(i,1) = inputArray(i,1)
	      outputArray(i,2) = inputArray(i,2)
		End If
	  End If
	End Do
  End Function VaryPoints 

    
!------------------------------------------------------------------------!
! Laplace Transform Functions
!------------------------------------------------------------------------! 
  
  Function GaverStehfestCoeffs(N) RESULT (coefficients)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: k,m,N,NHalf,kMin,kMax,signVal
	Real(kind=DoubleReal) :: coefficientSum
	Integer(kind=VeryLongInteger) :: factA, factB, factC, factD, factE, factF
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	!Integer(kind=VeryLongInteger)
!set variables
    NHalf = N/2
!Allocate coefficient array
    Allocate(coefficients(1:N))
	
	!Loop through coefficients
	Do m=1,N
	  coefficients(m) = 0.0D0
	  kMin = floor(1.0D0*(m+1)/2)
	  kMax = min(m,NHalf)
	  signVal = (-1.0D0)**(m+NHalf)
	  Do k=kMin,kMax
	    factA = Factorial(2*k)
	    factB = Factorial(NHalf-k)
	    factC = Factorial(k)
	    factD = Factorial(k-1)
	    factE = Factorial(m-k)
	    factF = Factorial(2*k-m)
	    coefficients(m)=coefficients(m)+1.0D0*((k**NHalf)*factA)/&
		  (factB*factC*factD*factE*factF)
	  End Do
	  coefficients(m)=signVal*coefficients(m)
	End Do
  
  End Function GaverStehfestCoeffs  
  
    
  
  
  
  
!------------------------------------------------------------------------!
! Decay Functions
!------------------------------------------------------------------------! 
  
   
  
  Function CalcIsotopeAmount(parentProductionRate,decayDataArray,tStart,&
  tEnd,tBeamEnd) RESULT (isotopeChange)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i,j,k,decaySteps,decayStepCounter
	Real(kind=DoubleReal) :: parentProductionRate,tStart,tEnd,tBeamEnd,dTime
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArrayTemp
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
	Real(kind=DoubleReal) :: w,lP,lA,lB,lC
	Real(kind=DoubleReal) :: lamlp !la-lp
	Real(kind=DoubleReal) :: lpmla !lp-la
	Real(kind=DoubleReal) :: lalamlalp !lala-lalp
	!Real(kind=DoubleReal) :: lalalalp !lala-lalp
	Real(kind=DoubleReal) :: eplP,enlP,eplA,enlA,eplB,enlB,eplC,enlC
	Real(kind=DoubleReal) :: eplAnlP,enlAnlP,eplAplP
	Real(kind=DoubleReal) :: nPstart,nAstart,nBstart,nCstart
	Real(kind=DoubleReal) :: nPend,nAend,nBend,nCend
	Real(kind=DoubleReal) :: ePA, eAB, eBC
	Real(kind=DoubleReal) :: stableLimit
	Real(kind=DoubleReal) :: termA, termB, termC, termD
	
	Integer(kind=StandardInteger) :: m,N,testILP
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: lnstore
	Real(kind=DoubleReal) :: sumN, s, Fs, ft
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: cS
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: nO
	Real(kind=DoubleReal) :: exactResult, iltResult, errorPercentage
	
!set variables
    dTime = tEnd - tStart
	
!Alter decay chain
! - If dTime * decay constant lt 1.0D-11 then assume stable for purposes of simulation
    decayStepCounter = 0
	Do i=1,size(decayDataArray,1)
	  stableLimit = (log(2.0D0)/decayDataArray(i,3))*dTime	  
	  decayStepCounter = decayStepCounter + 1
	  If(stableLimit.lt.1.0D-10)Then
	    decayDataArray(i,3) = -1		!set as stable
	    Exit
	  End If
	End Do
!Resize array
	decayDataArray = ArraySize2DDouble(decayDataArray,decayStepCounter)
 
!set decay steps/isotopes
    decaySteps = size(decayDataArray,1)
!allocate isotopeChange array	
	Allocate(isotopeChange(1:decaySteps,1:10))
	!isotopeChange(i,1)		!Tally key
	!isotopeChange(i,2)		!Change in isotope amount
	!isotopeChange(i,3)		!Start amount
	!isotopeChange(i,4)		!End amount
	!isotopeChange(i,5)		!Isotope Z
	!isotopeChange(i,6)		!Isotope A
	!isotopeChange(i,7)		!T1/2
	!isotopeChange(i,8)		!Decay constant
	!isotopeChange(i,9)		!Branching factor
	!isotopeChange(i,10)	!Parent production rate
	
!Fill with starting data
    Do i=1,decaySteps
	  isotopeChange(i,1) = decayDataArray(i,1)	
	  isotopeChange(i,2) = 0.0D0					!default no change
	  isotopeChange(i,3) = decayDataArray(i,2)	
	  isotopeChange(i,4) = decayDataArray(i,2)		!default no change
	  isotopeChange(i,5) = decayDataArray(i,5)	
	  isotopeChange(i,6) = decayDataArray(i,6)
	  isotopeChange(i,7) = decayDataArray(i,3)	
	  isotopeChange(i,8) = log(2.0D0)/decayDataArray(i,3)	
	  isotopeChange(i,9) = decayDataArray(i,4)	
	  isotopeChange(i,10) = parentProductionRate			
	End Do	
  
	w = parentProductionRate
!------------------
!---1 Decay Step	
!------------------
	If(decaySteps.eq.1)Then
	  nPstart = decayDataArray(1,2)
	  nPend = 1.0D0*(nPstart+1.0D0*dTime*w)	  
!force ge than 0
	  If(nPend.lt.0.0D0)Then
	    nPend = 0.0D0
	  End If
!Store key and tally change data	
	  isotopeChange(1,4) = nPend
	End If  
!------------------
!---2 Decay Steps	
!------------------	  
	If(decaySteps.eq.2)Then
!parent terms
	  lP = log(2.0D0)/decayDataArray(1,3)	 
	  enlP = exp(-1.0D0*dTime*lP)
	  nPstart = decayDataArray(1,2)
!child A terms	  
	  nAstart = decayDataArray(2,2)
	  ePA = decayDataArray(2,4) 
!calc nP
	  nPend = (1.0D0/lP)*(enlP*(nPstart*lP-w)+w)
!calc nA	
	  nAend = nAstart+ePA*((1-enlP)*nPstart+w*(dTime+(enlP-1)/lP))
!force ge than 0
	  If(nPend.lt.0.0D0)Then
	    nPend = 0.0D0
	  End If
	  If(nAend.lt.0.0D0)Then
	    nAend = 0.0D0
	  End If
!Store key and tally change data
	  isotopeChange(1,4) = nPend
	  isotopeChange(2,4) = nAend
	End If
!------------------
!---3 Decay Steps	
!------------------	  
	If(decaySteps.eq.3)Then
	  dTime = tEnd - tStart
	  w = parentProductionRate
!parent terms
	  lP = log(2.0D0)/decayDataArray(1,3)
	  eplP = exp(1.0D0*dTime*lP)
	  enlP = exp(-1.0D0*dTime*lP)
	  nPstart = decayDataArray(1,2)
!child A terms	  
	  nAstart = decayDataArray(2,2)
	  ePA = decayDataArray(2,4) 
	  lA = log(2.0D0)/decayDataArray(2,3)
	  eplA = exp(1.0D0*dTime*lA)
	  enlA = exp(-1.0D0*dTime*lA)
	  eplAnlP = exp(1.0D0*dTime*lA-1.0D0*dTime*lP)
!child B terms	  
	  nBstart = decayDataArray(3,2)
	  eAB = decayDataArray(3,4) 
	  lB = log(2.0D0)/decayDataArray(3,3)
	  eplA = exp(1.0D0*dTime*lA)
	  enlA = exp(-1.0D0*dTime*lA)
	  eplP = exp(1.0D0*dTime*lP)
	  enlP = exp(-1.0D0*dTime*lP)
	  enlAnlP = exp(-1.0D0*dTime*lA-1.0D0*dTime*lP)
	  eplAplP = exp(1.0D0*dTime*lA+1.0D0*dTime*lP)
!mixed terms
      lamlp = lA-lP
      lpmla = lP-lA
      lalamlalp = lA*lA-lA*lP	  
!calc nP
	  nPend = (1.0D0/lP)*(enlP*(nPstart*lP-w)+w)
!calc nA	
	  nAend = (1/(lA*(lA-lP)))*(nAstart*(lA*(lA-lP))+ePA*((enlA-1)*w*lP+&
	  lA*(w-nPstart*lP*enlA+enlP*(nPstart*lP-w))))	
!calc nB
	  nBend = (1/(lA*(lA-lP)*lP))*((lA-lP)*(-1.0D0*w*eAB*ePA*lA+&
	  (-1.0D0*w*eAB*ePA+(nBstart+eAB*(nAstart+(dTime*w+nPstart)*ePA))*&
	  lA)*lP)+eAB*(enlP*ePA*lA*lA*(w-nPstart*lP)+enlA*lP*(ePA*(-1.0D0*&
	  w+nPstart*lA)*lP+nAstart*lA*(lP-lA))))	
!force ge than 0
	  If(nPend.lt.0.0D0)Then
	    nPend = 0.0D0
	  End If
	  If(nAend.lt.0.0D0)Then
	    nAend = 0.0D0
	  End If
	  If(nBend.lt.0.0D0)Then
	    nBend = 0.0D0
	  End If  
!Store key and tally change data		
	  isotopeChange(1,4) = nPend
	  isotopeChange(2,4) = nAend
	  isotopeChange(3,4) = nBend	  
	End If
!------------------
!---4+ Decay Steps	
!------------------	  
	If(decaySteps.ge.4)Then
	  dTime = tEnd - tStart
	  w = parentProductionRate
!parent terms
	  lP = log(2.0D0)/decayDataArray(1,3)
	  eplP = exp(1.0D0*dTime*lP)
	  enlP = exp(-1.0D0*dTime*lP)
	  nPstart = decayDataArray(1,2)
!child A terms	  
	  nAstart = decayDataArray(2,2)
	  ePA = decayDataArray(2,4) 
	  lA = log(2.0D0)/decayDataArray(2,3)
	  eplA = exp(1.0D0*dTime*lA)
	  enlA = exp(-1.0D0*dTime*lA)
	  eplAnlP = exp(1.0D0*dTime*lA-1.0D0*dTime*lP)
!child B terms	  
	  nBstart = decayDataArray(3,2)
	  eAB = decayDataArray(3,4) 
	  lB = log(2.0D0)/decayDataArray(3,3)
	  eplA = exp(1.0D0*dTime*lA)
	  enlA = exp(-1.0D0*dTime*lA)
	  eplP = exp(1.0D0*dTime*lP)
	  enlP = exp(-1.0D0*dTime*lP)
	  enlAnlP = exp(-1.0D0*dTime*lA-1.0D0*dTime*lP)
	  eplAplP = exp(1.0D0*dTime*lA+1.0D0*dTime*lP)
!mixed terms
      lamlp = lA-lP
      lpmla = lP-lA
      lalamlalp = lA*lA-lA*lP	  
!calc nP
	  nPend = (1.0D0/lP)*(enlP*(nPstart*lP-w)+w)
!calc nA	
	  nAend = (1/(lA*(lA-lP)))*(nAstart*(lA*(lA-lP))+ePA*((enlA-1)*w*lP+&
	  lA*(w-nPstart*lP*enlA+enlP*(nPstart*lP-w))))	
!calc nB	  		  
	  termA=(1.0D0*exp(-1.0D0*dTime*(lA+lP)))/(lB*(lA-lB)*(lA-lP)*(lP-lB))
	  termB=exp(dTime*(lA-lB+lP))*nBstart*lB*(lB-lA)*(lA-lP)*(lB-lP)
	  termC=eAB*((exp(dTime*lP)-exp(dTime*(lA-lB+lP)))*nAstart*lA*lB*&
	  (lA-lP)*(lB-lP)+ePA*(exp(dTime*lP)*(exp(dTime*lA)-1.0D0)*w*lB*lP*&
	  (lP-lB)+lA*lA*(exp(dTime*(lA-lB+lP))*((exp(dTime*lB)-1)*w+&
	  nPstart*lB)*lP-exp(dTime*lA)*lB*((-1.0D0+exp(dTime*lP))*w+&
	  nPstart*lP))+lA*(-1.0D0*exp(dTime*(lA-lB+lP))*(-1.0D0+&
	  exp(dTime*lB))*w*lP*lP+(exp(dTime*lP)-exp(dTime*(lA-lB+lP)))*&
	  nPstart*lB*lP*lP+lB*lB*(exp(dTime*lA)*(-1+exp(dTime*lP))*w+&
	  (exp(dTime*lA)-exp(dTime*lP))*nPstart*lP))))
	  nBend = termA*(termB+termC)
!force ge than 0
	  If(nPend.lt.0.0D0)Then
	    nPend = 0.0D0
	  End If
	  If(nAend.lt.0.0D0)Then
	    nAend = 0.0D0
	  End If
	  If(nBend.lt.0.0D0)Then
	    nBend = 0.0D0
	  End If    
!Store key and tally change data		
	  isotopeChange(1,4) = nPend
	  isotopeChange(2,4) = nAend
	  isotopeChange(3,4) = nBend
!use numerical inverse laplace transform for 4 steps and more
!Gaver Stehfest coefficients
      N = 12
      Allocate(coefficients(1:12))
      coefficients(1) = -0.16666666666666666666D-1
      coefficients(2) = 0.160166666666666666666D2
      coefficients(3) = -0.1247D4
      coefficients(4) = 0.27554333333333333333D5
      coefficients(5) = -0.26328083333333333333D6
      coefficients(6) = 0.13241387D7
      coefficients(7) = -0.38917055333333333333D7
      coefficients(8) = 0.7053286333333333333333D7
      coefficients(9) = -0.80053365D7
      coefficients(10) = 0.55528305D7
      coefficients(11) = -0.21555072D7
      coefficients(12) = 0.3592512D6
!Allocate array for laplace transform equation coefficients
      Allocate(cS(0:decaySteps))
!starting nO
      Allocate(nO(1:decaySteps))
!------------------
!---4 Decay Steps	
!------------------	 
!calculate for 4 decay steps (4th isotope is stable)
	  If(decaySteps.eq.4)Then
!set atom start amounts
        Do i=1,decaySteps
          nO(i) = isotopeChange(i,3)
		End Do
!set other variables
	    lnstore = 0.693147180559945/dTime
	    sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  !lambda isotopeChange(i,8), !epsilon isotopeChange(i,9), !omega isotopeChange(i,10)
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = (1.0D0/s)
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(4,4) = lnstore * sumN		!store end atom count		
	  End If !End 4
!------------------
!---5 Decay Steps	
!------------------	 
!calculate for 5 decay steps (4th isotope is unstable, 5th is stable)
	  If(decaySteps.eq.5)Then	
	    Do i=1,decaySteps
          nO(i) = isotopeChange(i,3)
		End Do
!set other variables
	    lnstore = 0.693147180559945/dTime
!4th isotope Unstable
	    sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = (1.0D0/(s+isotopeChange(4,8)))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(4,4) = lnstore * sumN		!store end atom count		
!5th isotope Stable
		sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = 1.0D0*((isotopeChange(4,9)*isotopeChange(4,8))/(s+isotopeChange(4,8)))
          cS(5) = (1.0D0/s)
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(5)*(cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))+nO(5))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(5,4) = lnstore * sumN		!store end atom count	
	  End If !End 5
!------------------
!---6 Decay Steps	
!------------------	 
!calculate for 6 decay steps (4th isotope is unstable, 5th is unstable, 6th is stable)
	  If(decaySteps.eq.6)Then	
	    Do i=1,decaySteps
          nO(i) = isotopeChange(i,3)
		End Do
!set other variables
	    lnstore = 0.693147180559945/dTime
!4th isotope Unstable
	    sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients		  
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = (1.0D0/(s+isotopeChange(4,8)))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)+nO(1))+nO(2))+nO(3))+nO(4))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(4,4) = lnstore * sumN		!store end atom count		
!5th isotope Unstable
		sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = 1.0D0*((isotopeChange(4,9)*isotopeChange(4,8))/(s+isotopeChange(4,8)))
          cS(5) = (1.0D0/(s+isotopeChange(5,8)))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(5)*(cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))+nO(5))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(5,4) = lnstore * sumN		!store end atom count	
!6th isotope Stable
		sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients		  
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = 1.0D0*((isotopeChange(4,9)*isotopeChange(4,8))/(s+isotopeChange(4,8)))
          cS(5) = 1.0D0*((isotopeChange(5,9)*isotopeChange(5,8))/(s+isotopeChange(5,8)))
          cS(6) = (1.0D0/(s))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(6)*(cS(5)*(cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))+nO(5))+nO(6))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(6,4) = lnstore * sumN		!store end atom count	
	  End If !End 6
!------------------
!---7 Decay Steps	
!------------------	 
!calculate for 7 decay steps (4th isotope is unstable, 5th is unstable, 6th is stable)
	  If(decaySteps.eq.7)Then	
	    Do i=1,decaySteps
          nO(i) = isotopeChange(i,3)
		End Do
!set other variables
	    lnstore = 0.693147180559945/dTime
!4th isotope Unstable
	    sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = (1.0D0/(s+isotopeChange(4,8)))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)+nO(1))+nO(2))+nO(3))+nO(4))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(4,4) = lnstore * sumN		!store end atom count		
!5th isotope Unstable
		sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = 1.0D0*((isotopeChange(4,9)*isotopeChange(4,8))/(s+isotopeChange(4,8)))
          cS(5) = (1.0D0/(s+isotopeChange(5,8)))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(5)*(cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))+nO(5))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(5,4) = lnstore * sumN		!store end atom count	
!6th isotope Stable
		sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = 1.0D0*((isotopeChange(4,9)*isotopeChange(4,8))/(s+isotopeChange(4,8)))
          cS(5) = 1.0D0*((isotopeChange(5,9)*isotopeChange(5,8))/(s+isotopeChange(5,8)))
          cS(6) = (1.0D0/(s+isotopeChange(6,8)))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(6)*(cS(5)*(cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))+nO(5))+nO(6))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(6,4) = lnstore * sumN		!store end atom count	
!7th isotope Stable
		sumN = 0.0D0
	    Do m=1,N
	      s = lnstore*m
!Store laplace transform equation coefficients
		  cS(0) = 1.0D0*(isotopeChange(1,10)/s)
          cS(1) = 1.0D0*((isotopeChange(1,9)*isotopeChange(1,8))/(s+isotopeChange(1,8)))
          cS(2) = 1.0D0*((isotopeChange(2,9)*isotopeChange(2,8))/(s+isotopeChange(2,8)))
          cS(3) = 1.0D0*((isotopeChange(3,9)*isotopeChange(3,8))/(s+isotopeChange(3,8)))
          cS(4) = 1.0D0*((isotopeChange(4,9)*isotopeChange(4,8))/(s+isotopeChange(4,8)))
          cS(5) = 1.0D0*((isotopeChange(5,9)*isotopeChange(5,8))/(s+isotopeChange(5,8)))
          cS(6) = 1.0D0*((isotopeChange(6,9)*isotopeChange(6,8))/(s+isotopeChange(6,8)))
          cS(7) = (1.0D0/(s))
!--------------------------------------------------------------------------------------------------- 
	      Fs = 1.0D0*cS(7)*(cS(6)*(cS(5)*(cS(4)*(cS(3)*(cS(2)*(cS(1)*(cS(0)&
		  +nO(1))+nO(2))+nO(3))+nO(4))+nO(5))+nO(6))+nO(7))
!---------------------------------------------------------------------------------------------------
          sumN = sumN + coefficients(m) * Fs 
        End Do	
	    isotopeChange(7,4) = lnstore * sumN		!store end atom count
	  End If !End 7
	End If		
!force ge than 0
    Do i=1,size(isotopeChange,1)
	  If(isotopeChange(i,4).lt.0.0D0)Then
	    isotopeChange(i,4) = 0.0D0
	  End If
	End Do	
!Store changes in isotope amounts
	Do i=1,size(isotopeChange,1)
	  isotopeChange(i,2) = isotopeChange(i,4) - isotopeChange(i,3)
	End Do
  End Function CalcIsotopeAmount  
    
  
  
  
  
!------------------------------------------------------------------------!
! Rounding Functions
!------------------------------------------------------------------------! 
  
  
  Function Ceil (input) RESULT (output)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: tempInt, output
	Real(kind=DoubleReal) :: input
    output = input
	If(input.gt.(1.0D0*output))Then
	  output = output + 1
	End If
  End Function Ceil 
  
  
  
  
!------------------------------------------------------------------------!
! Miscellaneous Functions
!------------------------------------------------------------------------! 
  
  
  Function ArraySize1DDouble (inputArray,arraySize) RESULT (outputArray)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: arraySize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: inputArray
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
!Allocate output array
    Allocate(outputArray(1:arraySize))
!transfer data
    Do i=1,arraySize
	  If(i.le.size(inputArray))Then
	    outputArray(i) = inputArray(i)
	  Else
        outputArray(i) = 0.0D0	  
	  End If 	
	End Do  
  End Function ArraySize1DDouble 
  
  Function ArraySize2DDouble (inputArray,arraySizeHeight,arraySizeWidthIn) &
    RESULT (outputArray)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i, j
	Integer(kind=StandardInteger) :: arraySizeHeight
	Integer(kind=StandardInteger), optional :: arraySizeWidthIn
	Integer(kind=StandardInteger) :: arraySizeWidth
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputArray
!catch optional width
	if(present(arraySizeWidthIn))then
	  arraySizeWidth = arraySizeWidthIn
	else
	  arraySizeWidth = size(inputArray,2)
	endif
!Allocate output array
    Allocate(outputArray(1:arraySizeHeight,1:arraySizeWidth))
!transfer data
    Do i=1,arraySizeHeight
	  Do j=1,arraySizeWidth
	    If(i.le.size(inputArray,1).and.j.le.size(inputArray,2))Then
	      outputArray(i,j) = inputArray(i,j)
	    Else
          outputArray(i,j) = 0.0D0	  
	    End If 
	  End Do	
	End Do  
  End Function ArraySize2DDouble 
  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
! String functions
!------------------------------------------------------------------------! 
  
  Function StrToUpper (input) RESULT (output)
    ! -- Argument and result
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: output
	character( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, n
    ! -- Copy input string
    output = input
    ! -- Loop over string elements
    Do i=1,LEN(output)
      n = INDEX( LOWER_CASE, output( i:i ) )
      If (n/=0) output( i:i ) = UPPER_CASE( n:n )
    End Do
  End Function StrToUpper 
  
  
  

!------------------------------------------------------------------------!
! Miscellaneous Subroutines
!------------------------------------------------------------------------!   
  
  Subroutine swapMatrixRows(matrix,rowA,rowB) 
!Swap rows of square dp matrix
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
	Real(kind=DoubleReal), Dimension( : , :), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
!Set variables
	matH = size(matrix,1)
	matW = size(matrix,2)
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
!Allocate arrays
	  Allocate(rowAArr(1:matW))
	  Allocate(rowBArr(1:matW))
!Swap rows
	  Do i=1,matW
	    rowAArr(i) = matrix(rowA,i)
	    rowBArr(i) = matrix(rowB,i)
	  End Do
	  Do i=1,matW
	    matrix(rowA,i) = rowBArr(i)
	    matrix(rowB,i) = rowAArr(i)
	  End Do
    End If
  End Subroutine swapMatrixRows
  
    
  Subroutine swapMatrixRows1D(matrix,rowA,rowB) 
!Swap rows of square dp matrix
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
!Set variables
	matH = size(matrix,1)
	matW = 1
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
!Allocate arrays
	  Allocate(rowAArr(1:matW))
	  Allocate(rowBArr(1:matW))
!Swap rows
	  Do i=1,matW
	    rowAArr(i) = matrix(rowA)
	    rowBArr(i) = matrix(rowB)
	  End Do
	  Do i=1,matW
	    matrix(rowA) = rowBArr(i)
	    matrix(rowB) = rowAArr(i)
	  End Do
    End If
  End Subroutine swapMatrixRows1D
  
  
  
  Subroutine swapMatrixRows2D(matrix,rowA,rowB) 
!Swap rows of square dp matrix
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
	Real(kind=DoubleReal), Dimension( : , :), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
!Set variables
	matH = size(matrix,1)
	matW = size(matrix,2)
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
!Allocate arrays
	  Allocate(rowAArr(1:matW))
	  Allocate(rowBArr(1:matW))
!Swap rows
	  Do i=1,matW
	    rowAArr(i) = matrix(rowA,i)
	    rowBArr(i) = matrix(rowB,i)
	  End Do
	  Do i=1,matW
	    matrix(rowA,i) = rowBArr(i)
	    matrix(rowB,i) = rowAArr(i)
	  End Do
    End If
  End Subroutine swapMatrixRows2D
  
  
  
End Module maths  