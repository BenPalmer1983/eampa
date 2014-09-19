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
  Public :: ComputePolynomial
  Public :: ComputePolynomialA
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
  Public :: SplineAB
  Public :: MatrixInterpolation
! Matrix Functions  
  Public :: InvertSquareMatrix
  Public :: MatAdd
! Unit Vector Functions  
  Public :: ColToSquare
  Public :: SquareToCol
! Strain Matrices
  Public :: HomogeneousStrain  
  Public :: OrthorhombicStrain  
  Public :: MonoclinicStrain   
  Public :: TetragonalStrain  
! Vector Functions
  Public :: CrossProduct  
  Public :: DotProduct  
! Random Number Related Functions  
  Public :: SetRandomSeedArray !Subroutine
  Public :: VaryPoint
  Public :: VaryPoints
  Public :: VaryPointRand
  Public :: RandomInteger
  Public :: NumberList
  Public :: RandomFloat
! Laplace Transforms  
  Public :: GaverStehfestCoeffs
! Decay Equations    
  Public :: CalcIsotopeAmount  
! Rounding functions
  Public :: Ceil    
  Public :: Rounding
  Public :: ForceZero
! Physics type functions
  Public :: Zbl
  Public :: ZblFull
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
! Rounding functions
! - Ceil 
! - Rounding
!
! Physics type functions
! - Zbl
! 
! Useful Miscellaneous Functions
! - change ArraySize1D
! - change ArraySize2D

  

    
  
!------------------------------------------------------------------------!
! Interpolation and Regression Functions  
!------------------------------------------------------------------------! 
  
  
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
    inputSetLengthIn,interpTypeIn,verboseIn) RESULT (yArray) 
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
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray, yArrayTemp
  Integer(kind=StandardInteger) :: interpType
!Optional vars
  Integer(kind=StandardInteger), optional :: subsetSizeIn
  Integer(kind=StandardInteger), optional :: inputSetStartIn
  Integer(kind=StandardInteger), optional :: inputSetLengthIn
  Integer(kind=StandardInteger), optional :: interpTypeIn
  Integer(kind=StandardInteger), optional :: verboseIn
!Sort vars
  Integer(kind=StandardInteger) :: orderData
  Real(kind=DoubleReal) :: sortXA, sortXB, sortYA, sortYB
  Logical :: sorting
!Subset vars  
  Integer(kind=StandardInteger) :: subsetSize, dataSize, xPos, posOffset
  Real(kind=DoubleReal) :: xLower, xUpper
!Set variables  
    interpType = 1    !1 lagrange, 2 matrix
  If(Present(interpTypeIn))Then
    interpType = interpTypeIn
  End If
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
  If(Present(verboseIn))Then
    verbose = verboseIn
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
  orderData = 0
    If(orderData.eq.1)Then
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
    If(orderData.eq.1)Then
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
    If(interpType.eq.1)Then     !If Lagrange formula
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
  End If   !If Lagrange formula end
    If(interpType.eq.2)Then     !If Lagrange formula
    Allocate(yArray(1:3))
    yArrayTemp = MatrixInterpolation(inputPointsInterp,x)
    yArray(1) = yArrayTemp(2)
    yArray(2) = yArrayTemp(3)
    yArray(3) = yArrayTemp(4)
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
    splinePoints(j,1) = inputPoints(i,1)  !x
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


  
  Function SplineAB(points) RESULT (coefficients) 
! *5th order only*
! Spline between two points with 5th order polynomial
! points 1-8 xa, ya, dya, ddya, xb, yb, dyb, ddyb
!force declaration of all variables
  Implicit None
!declare variables
  Integer(kind=StandardInteger) :: i,j,k
  Integer(kind=StandardInteger) :: expon, row, column
  Integer(kind=StandardInteger) :: point, pointA, pointB
  Real(kind=DoubleReal) :: x
  Real(kind=DoubleReal), Dimension(1:8) :: points
  Real(kind=DoubleReal), Dimension(0:5) :: coefficients
  Real(kind=DoubleReal), Dimension(1:6,1:6) :: xMatrix, xMatrixInverse
  Real(kind=DoubleReal), Dimension(1:6) :: yMatrix, cMatrix
!sorting
    Integer(kind=StandardInteger) :: sortingExpon
    Real(kind=DoubleReal), Dimension(1:6) :: sortingValues  
  Logical :: sorting  
!Make xMatrix
!Row 1
  x = points(1)
  Do i=1,6
    xMatrix(1,i) = 1.0D0*x**(i-1)
  End Do
!Row 2
  x = points(5)
  Do i=1,6
    xMatrix(2,i) = 1.0D0*x**(i-1)
  End Do
!Row 3
  x = points(1)
  xMatrix(3,1) = 0.0D0
  Do i=2,6
    xMatrix(3,i) = 1.0D0*(i-1)*x**(i-2)
  End Do
!Row 4
  x = points(5)
  xMatrix(4,1) = 0.0D0
  Do i=2,6
    xMatrix(4,i) = 1.0D0*(i-1)*x**(i-2)
  End Do
!Row 5
  x = points(1)
  xMatrix(5,1) = 0.0D0
  xMatrix(5,2) = 0.0D0
  Do i=3,6
    xMatrix(5,i) = 1.0D0*(i-2)*(i-1)*x**(i-3)
  End Do
!Row 6
  x = points(5)
  xMatrix(6,1) = 0.0D0
  xMatrix(6,2) = 0.0D0
  Do i=3,6
    xMatrix(6,i) = 1.0D0*(i-2)*(i-1)*x**(i-3)
  End Do
!Make yMatrix
    yMatrix(1) = points(2)
    yMatrix(2) = points(6)
    yMatrix(3) = points(3)
    yMatrix(4) = points(7)
    yMatrix(5) = points(4)
    yMatrix(6) = points(8)
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
      Call swapMatrixRows1DA(sortingValues,i,i+1,size(sortingValues,1))
      Call swapMatrixRows2DA(xMatrix,i,i+1,size(xMatrix,1),size(xMatrix,2))
      Call swapMatrixRows1DA(yMatrix,i,i+1,size(yMatrix,1))
      End If
    End Do  
  End Do
!Solve equation XC=Y
    xMatrixInverse = InvertSquareMatrixA(xMatrix,6)
  cMatrix = matmul(xMatrixInverse,yMatrix)
!Copy to coefficients matrix
    Do i=0,5  
    coefficients(i) = cMatrix(i+1)
  End Do
  End Function SplineAB
  
  
   
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
! Unit Vector Functions  
!------------------------------------------------------------------------!   



!------------------------------------------------------------------------!
! Strain Matrices
!------------------------------------------------------------------------!  
  
  Function HomogeneousStrain(inputMatrix, strainAmount) RESULT (outputMatrix)
!force declaration of all variables
  Implicit None  
!declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal) :: strainAmount
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: inputMatrix, outputMatrix  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainMatrix, identityMatrix
!Strain matrix from Mehl et al 1990 Phys Rev B Vol 41 No. 15 pg  311-323
    strainMatrix = 0.0D0
    strainMatrix(1,1) = strainAmount
    strainMatrix(2,2) = strainAmount  
    strainMatrix(3,3) = strainAmount
!set identity matrix
    Do i=1,3
    Do j=1,3
      If(i.eq.j)Then
          identityMatrix(i,j) = 1.0D0  
        Else
        identityMatrix(i,j) = 0.0D0  
      End If
    End Do
  End Do
  strainMatrix = MatAdd(identityMatrix,strainMatrix,3,3)
  outputMatrix = MatMul(strainMatrix,inputMatrix)
  End Function HomogeneousStrain
    
  Function OrthorhombicStrain(inputMatrix, strainAmount) RESULT (outputMatrix)
!force declaration of all variables
  Implicit None  
!declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal) :: strainAmount
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: inputMatrix, outputMatrix  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainMatrix, identityMatrix
!Strain matrix from Mehl et al 1990 Phys Rev B Vol 41 No. 15 pg  311-323
    strainMatrix = 0.0D0
    strainMatrix(1,1) = 1.0D0*strainAmount
    strainMatrix(2,2) = -1.0D0*strainAmount  
    strainMatrix(3,3) = (strainAmount**2)/(1.0D0-strainAmount**2)
!set identity matrix
    Do i=1,3
    Do j=1,3
      If(i.eq.j)Then
          identityMatrix(i,j) = 1.0D0  
        Else
        identityMatrix(i,j) = 0.0D0  
      End If
    End Do
  End Do
  strainMatrix = MatAdd(identityMatrix,strainMatrix,3,3)
  outputMatrix = MatMul(strainMatrix,inputMatrix)
  End Function OrthorhombicStrain
    
  Function MonoclinicStrain(inputMatrix, strainAmount) RESULT (outputMatrix)
!force declaration of all variables
  Implicit None  
!declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal) :: strainAmount
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: inputMatrix, outputMatrix  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainMatrix, identityMatrix
!Strain matrix from Mehl et al 1990 Phys Rev B Vol 41 No. 15 pg  311-323
    strainMatrix = 0.0D0
    strainMatrix(1,2) = 0.5D0*strainAmount
    strainMatrix(2,1) = 0.5D0*strainAmount  
    strainMatrix(3,3) = (strainAmount**2)/(4-(strainAmount**2))
!set identity matrix
    Do i=1,3
    Do j=1,3
      If(i.eq.j)Then
          identityMatrix(i,j) = 1.0D0  
        Else
        identityMatrix(i,j) = 0.0D0  
      End If
    End Do
  End Do
  strainMatrix = MatAdd(identityMatrix,strainMatrix,3,3)
  outputMatrix = MatMul(strainMatrix,inputMatrix)
  End Function MonoclinicStrain
    
  Function TetragonalStrain(inputMatrix, strainAmount) RESULT (outputMatrix)
!force declaration of all variables
  Implicit None  
!declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal) :: strainAmount
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: inputMatrix, outputMatrix  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainMatrix, identityMatrix
!Strain matrix from Mehl et al 1990 Phys Rev B Vol 41 No. 15 pg  311-323
    strainMatrix = 0.0D0
    strainMatrix(1,1) = strainAmount
    strainMatrix(2,2) = strainAmount  
    strainMatrix(3,3) = 1.0D0/(1.0D0+strainAmount)**2-1.0D0
!set identity matrix
    Do i=1,3
    Do j=1,3
      If(i.eq.j)Then
          identityMatrix(i,j) = 1.0D0  
        Else
        identityMatrix(i,j) = 0.0D0  
      End If
    End Do
  End Do
  strainMatrix = MatAdd(identityMatrix,strainMatrix,3,3)
  outputMatrix = MatMul(strainMatrix,inputMatrix)
  End Function TetragonalStrain
  
!------------------------------------------------------------------------!
! Vector Functions
!------------------------------------------------------------------------!   

 
  
  
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
  
  
  Function VaryPoint(x,T,randLargeVariationIn) RESULT (xV)  
!force declaration of all variables
  Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i, distributionDataPoints
  Integer(kind=StandardInteger) :: factor
  Real(kind=DoubleReal) :: x,T,xV
  Real(kind=DoubleReal) :: distX, distY, distXInterval
  Real(kind=DoubleReal) :: sigma, mu, piNum, normalisationFactor
  Real(kind=DoubleReal) :: randNumber, variation, largeVariation
  Real(kind=DoubleReal), Dimension(:,:), Allocatable :: distributionPoints
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
  Logical :: addVariation
!Optional variables
    Integer(kind=StandardInteger), optional :: randLargeVariationIn
    Integer(kind=StandardInteger) :: randLargeVariation
!Set Variables  
  If(Present(randLargeVariationIn))Then
    randLargeVariation = randLargeVariationIn
  Else
      randLargeVariation = 0  
  End If
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
  largeVariation = 1.0D0
!Occasionally add/subtract a large variation
    If(randLargeVariation.gt.0)Then
    Call RANDOM_NUMBER(randNumber)
    yArray = PointInterpolationArr(distributionPoints,randNumber,5)
    largeVariation = yArray(1)
    Call RANDOM_NUMBER(randNumber)
    If((1.0D0*randNumber).ge.0.95D0)Then
      largeVariation = (2.0D0*largeVariation)/T
    End If
  End If  
!New point
    If(addVariation.eqv..true.)Then
      xV = x*(1+variation*T*largeVariation)
    Else
    xV = x*(1-variation*T*largeVariation)
  End If  
  End Function VaryPoint  
  
  
  Function VaryPoints(inputArray,T,pointIn,gradIn) RESULT (outputArray)  
!force declaration of all variables
  Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i, point
  Real(kind=DoubleReal) :: y,T
  Integer(kind=StandardInteger) :: factor
  Real(kind=DoubleReal), Dimension(:,:), Allocatable :: inputArray
  Real(kind=DoubleReal), Dimension(:,:), Allocatable :: outputArray
  Integer(kind=StandardInteger), optional :: pointIn
  Real(kind=DoubleReal), Dimension(:,:), Allocatable, optional  :: gradIn
  Real(kind=DoubleReal), Dimension(:,:), Allocatable :: grad
!Optional in
    point = 0
    If(present(pointIn))Then
    point = pointIn
  End If
    If(present(gradIn))Then
    grad = gradIn
  End If
!Allocate output array  
  Allocate(outputArray(1:size(inputArray,1),1:2))
!Loop over points
    Do i=1,size(inputArray,1)    
    If(point.eq.0)Then
      y = inputArray(i,2)
      outputArray(i,1) = inputArray(i,1)
      outputArray(i,2) = VaryPoint(y,T,1)
    Else
      If(i.eq.point)Then
      y = inputArray(i,2)
        outputArray(i,1) = inputArray(i,1)
        outputArray(i,2) = VaryPoint(y,T,1)
    Else
        outputArray(i,1) = inputArray(i,1)
        outputArray(i,2) = inputArray(i,2)
    End If
    End If
  End Do
  End Function VaryPoints 
  
  
  Function VaryPointRand(x,T,varMax,randLargeVarMaxIn,&
  randLargeVariationProbIn,varTypeIn) RESULT (xV)  
! x - input value
! T - temperature of gaussian
! varMax - maximum variation 
! randLargeVarMaxIn - random larger variation
! randLargeVariationProbIn - probability of larger random variation 0-1
! varTypeIn - 0 absolute (add on fraction of vMax), 1 relative (frac vmax times x)
!force declaration of all variables
  Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i, distributionDataPoints
  Integer(kind=StandardInteger) :: factor
  Real(kind=DoubleReal) :: x,T,xV,varMax
  Real(kind=DoubleReal) :: distX, distY, distXInterval,yZero, yOne
  Real(kind=DoubleReal) :: sigma, mu, piNum, normalisationFactor
  Real(kind=DoubleReal) :: randNumber, variation, largeVariation
  Real(kind=DoubleReal), Dimension(:,:), Allocatable :: distributionPoints
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
  Logical :: addVariation
!Optional variables
    Real(kind=DoubleReal), optional :: randLargeVarMaxIn
    Real(kind=DoubleReal) :: randLargeVarMax
    Real(kind=DoubleReal), optional :: randLargeVariationProbIn
    Real(kind=DoubleReal) :: randLargeVariationProb
    Integer(kind=StandardInteger), optional :: varTypeIn
    Integer(kind=StandardInteger) :: varType
!Set Variables  
  If(Present(randLargeVarMaxIn))Then
    randLargeVarMax = randLargeVarMaxIn
  Else
      randLargeVarMax = 0.0D0  
  End If
  If(Present(randLargeVariationProbIn))Then
    randLargeVariationProb = randLargeVariationProbIn
  Else
      randLargeVariationProb = 0.01D0  
  End If
  If(Present(varTypeIn))Then
    varType = varTypeIn
  Else
      varType = 1  
  End If
!Temperature start T=0 sigma=0.05, T=1000 sigma=5    20-30 good range
    sigma = 0.05+4.95D-3*T
!Make distribution set
    distributionDataPoints = 21
  mu = 0.0D0  
!Allocate distribution data array
    Allocate(distributionPoints(1:distributionDataPoints,1:2))
!Set data points  
  distXInterval = 1.0D0/(distributionDataPoints-1)
  yZero = 1.0D0
  yOne = exp(-1*((1.0D0-mu)**2)/(2*sigma**2))
!Make array
    distributionPoints(1,1) = 0.0D0
    distributionPoints(1,2) = 1.0D0
    distributionPoints(distributionDataPoints,1) = 1.0D0
    distributionPoints(distributionDataPoints,2) = 0.0D0
    Do i=2,distributionDataPoints-1
    distX = (i-1)*1.0D0*distXInterval
    distributionPoints(i,1) = distX
    distributionPoints(i,2) = (exp(-1*((distX-mu)**2)/(2*sigma**2))-yOne)/(1-yOne)
  End Do
!Generate random variation
    Call RANDOM_NUMBER(randNumber)
  yArray = PointInterpolationArr(distributionPoints,randNumber,3)
  variation = yArray(1)
!Whether to add or subtract variation
    Call RANDOM_NUMBER(randNumber)
  If((1.0D0*randNumber).ge.0.5D0)Then
    variation = -1.0D0*variation
  End If
!set xV with standard variation
  If(varType.eq.0)Then
    xV= 1.0D0 * (x + variation*varMax)
  Else If(varType.eq.1)Then
    xV= x * (1.0D0 + variation*varMax)
  End If
!set xV with random large variation
  If(randLargeVarMax.gt.0.0D0)Then
    Call RANDOM_NUMBER(randNumber)
    If(randNumber.le.randLargeVariationProb)Then
      If(varType.eq.0)Then
        xV= 1.0D0 * (x + variation*randLargeVarMax)
      Else If(varType.eq.1)Then
        xV= x * (1.0D0 + variation*randLargeVarMax)
      End If
    End If
  End If   
  End Function VaryPointRand  
  
  Function RandomInteger(lower,upper) RESULT (randInt)
!force declaration of all variables
  Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i  
  Integer(kind=StandardInteger) :: lower, upper, randInt, diff, tempInt
  Real(kind=DoubleReal) :: randDouble, fractionInterval
!Make random integer
    If(lower.gt.upper)Then
    tempInt = lower
    lower = upper
    upper = tempInt
  End If
    diff = (upper - lower) + 1
  Call RANDOM_NUMBER(randDouble)
  randInt = lower+floor(1.0D0*diff*randDouble)  
  End Function RandomInteger  
  

  
  Function RandomFloat(lower,upper,distType,sigma) RESULT (output)
!force declaration of all variables
  Implicit None
!declare variables   
    Real(kind=DoubleReal) :: lower, upper
  Character(len=1) :: distType   !F flat, G Gaussian type, M Maxwell type
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:), Allocatable :: distributionPoints
    Real(kind=DoubleReal) :: yOne
    Real(kind=DoubleReal) :: sigma, mu
    Real(kind=DoubleReal) :: distX, distY, distXInterval
    Real(kind=DoubleReal) :: randNumber, output
  Real(kind=DoubleReal), Dimension(:), Allocatable :: yArray 
!Flat distribution  
    If(distType.eq."F")Then
!Generate random number
      Call RANDOM_NUMBER(randNumber)
      output = lower + (upper-lower) * randNumber
    End If
!Gaussian type distribution
    If(distType.eq."G")Then
    Allocate(distributionPoints(1:21,1:2))
    mu = 0.0D0
!Set data points  
    distXInterval = 1.0D0/(21-1)
    yOne = exp(-1*((1.0D0-mu)**2)/(2*sigma**2))
!Make array
      distributionPoints(1,1) = 0.0D0
      distributionPoints(1,2) = 1.0D0
      distributionPoints(21,1) = 1.0D0
      distributionPoints(21,2) = 0.0D0
      Do i=2,20
      distX = (i-1)*1.0D0*distXInterval
      distributionPoints(i,1) = distX
      distributionPoints(i,2) = (exp(-1*((distX-mu)**2)/(2*sigma**2))-yOne)/(1-yOne)
    End Do
!Generate random number
      Call RANDOM_NUMBER(randNumber)
    yArray = PointInterpolationArr(distributionPoints,randNumber,3)
    output = yArray(1)
      output = lower + (upper-lower) * output
    End If
  
  End Function RandomFloat  
  
    
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
      decayDataArray(i,3) = -1    !set as stable
      Exit
    End If
  End Do
!Resize array
  decayDataArray = ArraySize2DDouble(decayDataArray,decayStepCounter)
 
!set decay steps/isotopes
    decaySteps = size(decayDataArray,1)
!allocate isotopeChange array  
  Allocate(isotopeChange(1:decaySteps,1:10))
  !isotopeChange(i,1)    !Tally key
  !isotopeChange(i,2)    !Change in isotope amount
  !isotopeChange(i,3)    !Start amount
  !isotopeChange(i,4)    !End amount
  !isotopeChange(i,5)    !Isotope Z
  !isotopeChange(i,6)    !Isotope A
  !isotopeChange(i,7)    !T1/2
  !isotopeChange(i,8)    !Decay constant
  !isotopeChange(i,9)    !Branching factor
  !isotopeChange(i,10)  !Parent production rate
  
!Fill with starting data
    Do i=1,decaySteps
    isotopeChange(i,1) = decayDataArray(i,1)  
    isotopeChange(i,2) = 0.0D0          !default no change
    isotopeChange(i,3) = decayDataArray(i,2)  
    isotopeChange(i,4) = decayDataArray(i,2)    !default no change
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
      isotopeChange(4,4) = lnstore * sumN    !store end atom count    
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
      isotopeChange(4,4) = lnstore * sumN    !store end atom count    
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
      isotopeChange(5,4) = lnstore * sumN    !store end atom count  
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
      isotopeChange(4,4) = lnstore * sumN    !store end atom count    
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
      isotopeChange(5,4) = lnstore * sumN    !store end atom count  
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
      isotopeChange(6,4) = lnstore * sumN    !store end atom count  
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
      isotopeChange(4,4) = lnstore * sumN    !store end atom count    
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
      isotopeChange(5,4) = lnstore * sumN    !store end atom count  
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
      isotopeChange(6,4) = lnstore * sumN    !store end atom count  
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
      isotopeChange(7,4) = lnstore * sumN    !store end atom count
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
  
  
  Function Rounding (input, precisionVal) RESULT (output)
!force declaration of all variables
  Implicit None
!declare variables
  Integer(kind=StandardInteger) :: precisionVal
  Real(kind=DoubleReal) :: input, output
    output = ANINT(input*10**precisionVal)/(10**precisionVal)  
  End Function Rounding 
  
  
  Function ForceZero (input, threshold) RESULT (output)
!force declaration of all variables
  Implicit None
!declare variables
  Real(kind=DoubleReal) :: input, output, threshold
  If(abs(input).le.abs(threshold))Then
      output = 0.0D0
    Else
      output = input
    End If  
  End Function ForceZero 
  
  
    
  
!------------------------------------------------------------------------!
! Physics type functions
!------------------------------------------------------------------------! 
  
  
  Function Zbl (x, qA, qB) RESULT (y)
!force declaration of all variables
  Implicit None
!declare variables
    Integer(kind=StandardInteger) :: qA, qB
  Real(kind=DoubleReal) :: xVal, x, y, xa, xs, exa
!Force none infinite result for 0
    If(x.eq.0.0D0)Then
    xVal = 0.001D0
  Else 
    xVal = x
  End If
!Calculate y
  xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
  xa = 1.0D0*xVal/xs
  exa = 0.1818D0*exp(-3.2D0*xa)+0.5099D0*exp(-0.9423D0*xa)+&
        0.2802D0*exp(-0.4029*xa)+0.02817*exp(-0.2016D0*xa)
  y = ((1.0D0*qA*qB)/xVal)*exa  
  End Function Zbl 
  

  
  
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
  
  
  
  
  
  Function IntegerArray (totalIntegers) &
    RESULT (outputArray)
!force declaration of all variables
  Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i, totalIntegers
  Integer(kind=StandardInteger), Dimension(1:totalIntegers) :: outputArray
!Store numbers
    Do i=1,totalIntegers
      outputArray(i) = i
    End Do
  End Function IntegerArray 
  
  
  
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
  
  
  
!Non allocatable arrays
  Subroutine swapMatrixRows1DA(matrix,rowA,rowB,arraySizeH) 
!Swap rows of square dp matrix
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW, arraySizeH
  Real(kind=DoubleReal), Dimension(1:arraySizeH) :: matrix
    Real(kind=DoubleReal), Dimension(1:1) :: rowAArr
    Real(kind=DoubleReal), Dimension(1:1) :: rowBArr
!Set variables
  matH = arraySizeH
  matW = 1
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
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
  End Subroutine swapMatrixRows1DA 
  
  Subroutine swapMatrixRows2DA(matrix,rowA,rowB,matH,matW) 
!Swap rows of square dp matrix
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
  Real(kind=DoubleReal), Dimension(1:matH,1:matW) :: matrix
    Real(kind=DoubleReal), Dimension(1:matW) :: rowAArr
    Real(kind=DoubleReal), Dimension(1:matW) :: rowBArr
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
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
  End Subroutine swapMatrixRows2DA
  
  Subroutine swapMatrixRows1DInt(matrix,rowA,rowB,matH) 
!Swap rows of square dp matrix
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i, rowA, rowB, matH, tempA, tempB
  Integer(kind=StandardInteger), Dimension(1:matH) :: matrix
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
!Swap rows
      tempA = matrix(rowA)
      tempB = matrix(rowB)
    matrix(rowB) = tempA
    matrix(rowA) = tempB
    End If
  End Subroutine swapMatrixRows1DInt
  
  
  
End Module maths  