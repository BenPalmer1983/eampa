Module maths

! --------------------------------------------------------------!
! Maths functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! ----------------------------------------
! Updated: 8th July 2015
! ----------------------------------------

! Force declaration of all variables
  Implicit None
! Data types
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: QuadrupleReal = Selected_Real_Kind(33,4931)   ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer      
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer      0 to 32767
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer   0 to 2147483647
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer       0 to 9223372036854775807
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer  0 to 170141183460469231731687303715884105727
! Constants
  Real(kind=DoubleReal), Parameter :: lnTwo = 0.6931471805599453094172321214581765680755D0
  Real(kind=QuadrupoleReal), Parameter :: lnTwoQ = 0.6931471805599453094172321214581765680755D0
  Real(kind=DoubleReal), Parameter :: pi = 3.1415926535898D0
  Real(kind=DoubleReal), Parameter :: sqrtTwo = 1.4142135623731D0
  Real(kind=DoubleReal), Parameter :: sqrtTwoPi = 2.506628274631D0
! Variables
  Integer(kind=LongInteger) :: randomLCG_n=0  
  Integer(kind=LongInteger) :: randomLCG_xn  
  Real(kind=DoubleReal) :: randomDist_randNumberG = -1.0D0
  Real(kind=DoubleReal) :: randomDist_randNumberH = -1.0D0
  Real(kind=DoubleReal), Dimension(0:100,1:2) :: randomDist_inverseInt = 0.0D0
  Integer(kind=LongInteger) :: ddfRssCount

! Make private
  Private
! --Functions--!
! General Maths Functions
  Public :: Factorial
  Public :: BinomialCoefficient
  Public :: RoundDP
  Public :: logN2
  Public :: Gaussian
  Public :: MaxwellBoltzman
! Polynomial Related Functions
  Public :: SolvePolynomial
  Public :: CalcPolynomial
  Public :: DerivativePolynomial
! Fitting, Regression, Interpolation
  Public :: PolyFit             ! Fit polynomial to large set of data - Vandermonde
  Public :: SP_PolyFit          ! Fit polynomial to large set of data - Superposition
  Public :: PolyFitQ            ! Fit polynomial to large set of data - Vandermonde (use quadruple precision)
  Public :: MinPolyFit          ! Fit polynomial, find minimum (if in region of data points)
  Public :: InterpMatrix        ! Returns coefficients, interpolating set of points
  Public :: InterpMatrixPoint   ! Returns value of a point
  Public :: InterpMatrixPointInaccurate
  Public :: InterpLagrange      !find y(x) or y'(x) using lagrange interp from set of x-y
  Public :: PointInterp         !y(x) from large set of data x-y, finds region and uses lagrange interp
  Public :: CalcResidualSquareSum
  Public :: MurnFit
  Public :: MurnCalc
  Public :: MurnRSS
  Public :: MurnFitBP
  Public :: BirchMurnFit
  Public :: BirchMurnFit_R
  Public :: BirchMurnCalc
  Public :: BirchMurnRSS
  Public :: BirchMurnFitBP
  Public :: PolyPoints
  Public :: BirchMurnPoints
  Public :: SingleDecayFit
  Public :: DoubleDecayFit
! Optimization functions  
  Public :: NewtonGaussOpt
! Vector Functions
  Public :: CrossProduct
  Public :: DotProduct
  Public :: TripleProduct
  Public :: TripleProductSq
! Numbers  
  Public :: Odd
  Public :: Even
! Random Number Related Functions
  Public :: RandomInteger
  Public :: RandomFloat
  Public :: RandomLCG
  Public :: RandomDist
  Public :: RandomDist_GP
  Public :: RandomVaryPoint
  Public :: IntegerList
! Matrix Functions
  Public :: InvertMatrix
  Public :: TransposeMatrix
  Public :: DiagMatrix
  Public :: MatAdd
  Public :: MatMult
  Public :: ScalarMult
  Public :: InvertMatrixQ
! Coordinates
  Public :: TransformCoords
  Public :: RdCoords
  Public :: HeatCoords
! Spline Functions
  Public :: SplineAB
  Public :: SplineNodes
  Public :: SplineNodesV
  Public :: Spline
  Public :: VaryNode
  Public :: FillSplineResponse
! Other Physical Functions
  Public :: Zbl
  Public :: ZblFull
  Public :: CalcIsotopeAmount
  Public :: GaverStehfestWeighting
  Public :: GaverStehfestWeightingQ
  Public :: MaxTrajDepth

! --Subroutines--!
  Public :: setRandomSeedArray
  Public :: CompleteNodeData
  Public :: swapRows_Double_1D
  Public :: swapRows_Double_2D
  
! --Test Functions--!  
  Public :: expFit_NG  ! exp fit Newton Gauss
  Public :: expFit_LMA  ! exp fit Newton Gauss
  Public :: test_NG
  
! --variables--!  
  Public :: randomLCG_n
  Public :: randomLCG_xn 
  Public :: ddfRssCount

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE INTERFACE                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! Interface swapRows
!  Module Procedure swapRows_Integer_1D
! swapRows_Real, swapRows_Double
! End Interface swapRows

  Contains

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

! ------------------------------------------------------------------------!
! General Maths Functions
! ------------------------------------------------------------------------!

  Function Factorial(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=StandardInteger) :: output
! calculate factorial
    output = 1
    Do i=1,input
      output = i * output
    End Do
  End Function Factorial

  Function BinomialCoefficient(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: c,n,k
! calculate factorial
    c = Factorial(n)/(Factorial(n-k)*Factorial(k))
  End Function BinomialCoefficient

  Function FactorialDP(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: tempInt
    Real(kind=DoubleReal) :: output
! calculate factorial
    tempInt = 1
    Do i=1,input
      tempInt = i * tempInt
    End Do
    output = 1.0D0*tempInt
  End Function FactorialDP

  Function BinomialCoefficientDP(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=DoubleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientDP

  Function FactorialQ(input) RESULT (output)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=VeryLongInteger) :: tempInt
    Real(kind=QuadrupoleReal) :: tempQ
    Real(kind=QuadrupoleReal) :: output
! calculate factorial
    tempInt = 1
    tempQ = 1
    Do i=1,input
      If(i.le.33)Then
        tempInt = i * tempInt
        tempQ = 1.0D0*tempInt
      End If
      If(i.eq.34)Then
        tempQ = 1.0D0*i*tempInt
      End If
      If(i.ge.35)Then
        tempQ = 1.0D0*i*tempQ
      End If
    End Do
    output = tempQ
  End Function FactorialQ

  Function BinomialCoefficientQ(n,k) RESULT (c)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: n,k
    Real(kind=QuadrupoleReal) :: c, nDP, nkDP, kDP
! calculate factorial
    nDP = FactorialDP(n)
    nkDP = Factorial(n-k)
    kDP = FactorialDP(k)
    c = 1.0D0*nDP/(nkDP*kDP)
  End Function BinomialCoefficientQ

  Function RoundDP(dpIn) RESULT (intOut)
! Round DP to nearest int
    Implicit None ! Force declaration of all variables
    Real(kind=DoubleReal) :: dpIn
    Integer(kind=StandardInteger) :: intOut
    intOut = Floor(dpIn+0.5D0)
  End Function RoundDP

  Function logN2(loops) RESULT (output)
! Return quadruple value of ln(2)
    Implicit None ! Force declaration of all variables
    Integer(kind=StandardInteger) :: loops, k
    Real(kind=QuadrupoleReal) :: output
    output = 0.0D0
    Do k=0,loops
      output = output + (1.0D0/((2.0D0*k+1.0D0)*(2.0D0*k+2.0D0)))
    End Do
! Calculate Bailey et al 2007 method
! output = 0.0D0
! Do k=0,loops
!   output = output + (1.0D0/(9**k*(2.0D0*k+1.0D0)))
! End Do
! output = (2.0D0/3.0D0)*output
  End Function logN2
  
  
  Function Gaussian(x, sigma, mu) RESULT (y)
    Implicit None  !Force declaration of all variables
! Declare variables  
    Real(kind=DoubleReal) :: x,y,mu,sigma
! Calculation
    y = (1.0D0/(sigma*sqrtTwoPi))*exp(-1*((x-mu)**2/(2*sigma**2)))
  End Function Gaussian
  
  Function MaxwellBoltzman(x, a) RESULT (y)
    Implicit None  !Force declaration of all variables
! Declare variables  
    Real(kind=DoubleReal) :: x,y,a
! Calculation
    y = sqrtTwoPi*((x**2*exp(-1.0D0*((x*x)/(2.0D0*a*a))))/(a**3))
  End Function MaxwellBoltzman
  

! ------------------------------------------------------------------------!
! Polynomial Related Functions
! ------------------------------------------------------------------------!

  Function SolvePolynomial(coefficients, lower, upper, convergenceThresholdIn) RESULT (output)
! Solves the polynomial p(x) = 0, in region close to p(x) = 0
! For 3rd order polynomials or above (quadratic equation for 2nd order)
! Newton's method
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: upper, lower, output
    Real(kind=DoubleReal) :: xN,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Integer(kind=StandardInteger) :: maxLoops
! Optional argument
    convergenceThreshold = 1.0D-7
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
! Set values
    convergenceTarget = 0.0D0
    convergence = 1000
    maxLoops = 0
! set start value for x
    xN = lower + RandomLCG()*(upper-lower)    
    Do while(convergence.gt.convergenceThreshold.and.maxLoops.le.10000)
      maxLoops = maxLoops + 1
      y = CalcPolynomial(coefficients, xN, 0)
      dydx = CalcPolynomial(coefficients, xN, 1)
      xN = xN - (y/dydx)
      convergence = abs(convergenceTarget - y)
    End Do
    print *,maxLoops
    output = xN
  End Function solvePolynomial

  Function CalcPolynomial(polyCoefficientsIn, x, derivativeIn) RESULT (y)
! Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal), Dimension( : ) :: polyCoefficientsIn
    Real(kind=DoubleReal), Dimension(1:Size(polyCoefficientsIn,1)) :: polyCoefficients
    Real(kind=DoubleReal) :: x, y
! Handle Optional Arguments
    If(Present(derivativeIn))Then
      derivative = derivativeIn
      If(derivative.lt.0)Then
        derivative = 0
      End If
    Else
      derivative = 0
    End If
! Init variables
    polyCoefficients = polyCoefficientsIn
! Change coefficients (if required)
    If(derivative.gt.0)Then
      Do i=1,derivative
        Do j=1,size(polyCoefficients,1)-1
          polyCoefficients(j) = j * polyCoefficients(j+1)
        End Do
        polyCoefficients(size(polyCoefficients,1)+1-derivative) = 0.0D0
      End Do
    End If
! Calculate
    y = 0.0D0
    Do i=1,size(polyCoefficients,1)
      If(polyCoefficients(i).ne.0.0D0)Then
        If(i.eq.1)Then
          y = y + polyCoefficients(i)
        Else
          y = y + polyCoefficients(i)*x**(i-1)
        End If
      End If
    End Do
  End Function CalcPolynomial
  
  Function DerivativePolynomial (coefficientsIn) RESULT (coefficientsOut)
! Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : ) :: coefficientsIn
    Real(kind=DoubleReal), Dimension(1:(size(coefficientsIn)-1)) :: coefficientsOut
    Integer(kind=StandardInteger) :: j
    Do j=1,size(coefficientsIn,1)-1
      coefficientsOut(j) = j * coefficientsIn(j+1)
    End Do
  End Function DerivativePolynomial
  
  Function SolvePolynomialOld (coefficients, lower, upper, convergenceThresholdIn) RESULT (output)
! Solves the polynomial p(x) = 0, in region close to p(x) = 0
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: upper, lower, output
    Real(kind=DoubleReal) :: x,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget, factor, difference
    Real(kind=DoubleReal), Optional :: convergenceThresholdIn
    Integer(kind=StandardInteger) :: i,maxLoops
! Optional argument
    convergenceThreshold = 1.0D-5
    If(Present(convergenceThresholdIn))Then
      convergenceThreshold = convergenceThresholdIn
    End If
! Set values
    convergenceTarget = 0
    convergence = 1000
    maxLoops = 0
! set start value for x
    factor = 0.5
    difference = upper - lower
    x = lower + factor * difference
    Do while(convergence.gt.convergenceThreshold.and.maxLoops.le.10000)
      maxLoops = maxLoops + 1
      difference = factor * difference
      y = 0
      Do i=1,size(coefficients)
        y = y + x**(i-1) * coefficients(i)
      End Do
      dydx = 0
      Do i=2,size(coefficients)
        dydx = dydx + i * x**(i-2) * coefficients(i)
      End Do
      convergence = abs(convergenceTarget - y)
      If(convergence.gt.convergenceThreshold)Then
        If((dydx.lt.0.and.y.ge.0).or.(dydx.ge.0.and.y.lt.0))Then
          x = x + difference
        Else
          x = x - difference
        End If
      End If
    End Do
    output = x
  End Function solvePolynomialOld

  Function MinimaPolynomial (coefficients, lower, upper) RESULT (x)
! Finds minima in section of polynomial
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal), Dimension(1:size(coefficients,1)-1) :: coefficientsD
    Real(kind=DoubleReal) :: upper, lower
    Real(kind=DoubleReal) :: x, y, xInc, minX, minY
! Find section with minima
    xInc = (upper-lower)/100.0D0
    Do i=0,100
      x = lower + i * xInc
      y = CalcPolynomial(coefficients,x)
      If(i.eq.0)Then
        minX = x
        minY = y
      Else
        If(y.lt.minY)Then
          minY = y
          minX = x
        End If
      End If
    End Do
    coefficientsD = DerivativePolynomial(coefficients)
! Find minimum
    x = SolvePolynomial(coefficientsD,minX-xInc,minX+xInc)
  End Function MinimaPolynomial

! ------------------------------------------------------------------------!
! Fitting, Regression, Interpolation
! ------------------------------------------------------------------------!
  Function PolyFit(points,order,extendedFitIn) RESULT (coefficients)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Integer(kind=StandardInteger), Optional :: extendedFitIn
    Integer(kind=StandardInteger) :: extendedFit
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal) :: vRSS, spRSS
    Integer(kind=StandardInteger) :: order
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients, coefficientsSP
    Real(kind=DoubleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: yMatrix
! Optional argument
    extendedFit = 0
    If(Present(extendedFitIn))Then
      extendedFit = extendedFitIn
    End If
! Step 1 - Standard fit with Vandermonde matrix
! Build Least Squares Fitting Vandermonde matrix
    Do row=1,(order+1)
      Do col=1,(order+1)
        exponentValue = row+col-2
        xMatrix(row,col) = 0.0D0
        Do k=1,size(points,1)
          xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*points(k,1)&
          **exponentValue
        End Do
      End Do
    End Do
    Do row=1,(order+1)
      exponentValue = row-1
      yMatrix(row) = 0.0D0
      Do k=1,size(points,1)
        yMatrix(row) = 1.0D0*yMatrix(row)+1.0D0*points(k,2)*&
        points(k,1)**exponentValue
      End Do
    End Do
! invert xMatrix
    xMatrix = InvertMatrix(xMatrix)
! multiply inverse by y to get coefficients
    coefficients = matMul(xMatrix,yMatrix)
! Step 2 - extended fit
    If(extendedFit.ne.0)Then
! If required, use superposition fit and check if better than vandermonde
      vRSS = CalcResidualSquareSum(points,coefficients)
      coefficientsSP = SP_PolyFit(points,order,50)
      spRSS = CalcResidualSquareSum(points,coefficientsSP)
      If(spRSS.lt.vRSS)Then
        coefficients = coefficientsSP
      End If
    End If
  End Function PolyFit
  Function SP_PolyFit(points,order,iterationsIn) RESULT (coefficients)
! Fits multiple polynomials exactly to random data points and takes superposition of these as the result
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, order, iterations, spCount
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:(order+1),1:2) :: interpPoints
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients, coefficientsTemp, &
    coefficientsTrial
    Real(kind=DoubleReal) :: rss, lastRSS
    Integer(kind=StandardInteger), Dimension(1:(size(points,1)-2)) :: pointList
    Integer(kind=StandardInteger), Optional :: iterationsIn
! Initialise variables
    coefficients = 0.0D0
    coefficientsTemp = 0.0D0
    coefficientsTrial = 0.0D0
    interpPoints = 0.0D0
    spCount = 1
! Optional arguments
    iterations = 50
    If(Present(iterationsIn))Then
      iterations = iterationsIn
    End If
! Choose starting coefficients, best out of 3
    Do i=1,iterations
! Make set of interp points
      interpPoints(1,1) = points(1,1)   ! Tether first point
      interpPoints(1,2) = points(1,2)
      interpPoints((order+1),1) = points(size(points,1),1)  ! Tether last point
      interpPoints((order+1),2) = points(size(points,1),2)
! Make shuffled array of other points
      pointList = IntegerList(2,size(points,1)-1,20)
! Set remaining points
      Do j=2,order
        interpPoints(j,1) = points(pointList(j-1),1)
        interpPoints(j,2) = points(pointList(j-1),2)
      End Do
! Exactly fit polynomial
      coefficientsTemp = PolyFitExact(interpPoints)
      rss = CalcResidualSquareSum(points,coefficientsTemp)
! Store if an improvement
      If(i.eq.1)Then
        coefficients = coefficientsTemp
      Else
        If(rss.lt.lastRSS)Then
          coefficients = coefficientsTemp
        End If
      End If
! Store last RSS
      lastRSS = rss
    End Do
! Loop and make coefficients - with tethers
    Do i=1,iterations
! Make set of interp points
      interpPoints(1,1) = points(1,1)   ! Tether first point
      interpPoints(1,2) = points(1,2)
      interpPoints((order+1),1) = points(size(points,1),1)  ! Tether last point
      interpPoints((order+1),2) = points(size(points,1),2)
! Make shuffled array of other points
      pointList = IntegerList(2,size(points,1)-1,20)
! Set remaining points
      Do j=2,order
        interpPoints(j,1) = points(pointList(j-1),1)
        interpPoints(j,2) = points(pointList(j-1),2)
      End Do
! Exactly fit polynomial
      coefficientsTemp = PolyFitExact(interpPoints)
! Trial merge with existing coefficients
      Do j=1,size(coefficients,1)
        coefficientsTrial(j) = (spCount*coefficients(j)+coefficientsTemp(j))/(1.0D0*(spCount+1))
      End Do
      rss = CalcResidualSquareSum(points,coefficientsTrial)
! Store if an improvement
      If(rss.lt.lastRSS)Then
        coefficients = coefficientsTrial
        spCount = spCount + 1
      End If
      lastRSS = rss
    End Do
! Loop and make coefficients - without tethers
    Do i=1,iterations
! Make set of interp points
! Make shuffled array of points
      pointList = IntegerList(2,size(points,1)-1,20)
! Set remaining points
      Do j=1,order+1
        interpPoints(j,1) = points(pointList(j),1)
        interpPoints(j,2) = points(pointList(j),2)
      End Do
! Exactly fit polynomial
      coefficientsTemp = PolyFitExact(interpPoints)
! Trial merge with existing coefficients
      Do j=1,size(coefficients,1)
        coefficientsTrial(j) = (spCount*coefficients(j)+coefficientsTemp(j))/(1.0D0*(spCount+1))
      End Do
      rss = CalcResidualSquareSum(points,coefficientsTrial)
! Store if an improvement
      If(rss.lt.lastRSS)Then
        coefficients = coefficientsTrial
        spCount = spCount + 1
      End If
      lastRSS = rss
    End Do
  End Function SP_PolyFit
! -------------------------------------------------------------------------------------------------!
  Function PolyFitQ(points,order,extendedFitIn) RESULT (coefficientsOut)
! Fits a polynomial of order to the points input
! Input points and output coeffs DP, internal matrix QP
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Integer(kind=StandardInteger), Optional :: extendedFitIn
    Integer(kind=StandardInteger) :: extendedFit
    Integer(kind=StandardInteger) :: order
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficientsOut
! Real(kind=QuadrupleReal) :: vRSS, spRSS
    Real(kind=QuadrupleReal), Dimension(1:(order+1),1:1) :: coefficients
    Real(kind=QuadrupleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=QuadrupleReal), Dimension(1:(order+1),1:1) :: yMatrix
! Optional argument
    extendedFit = 0
    If(Present(extendedFitIn))Then
      extendedFit = extendedFitIn
    End If
! Step 1 - Standard fit with Vandermonde matrix
! Init matrices
    xMatrix = 0.0D0
    yMatrix = 0.0D0
    coefficientsOut = 0.0D0
! Build Least Squares Fitting Vandermonde matrix
    Do row=1,(order+1)
      Do col=1,(order+1)
        exponentValue = row+col-2
        Do k=1,size(points,1)
          xMatrix(row,col) = 1.0D0*xMatrix(row,col)+1.0D0*points(k,1)&
          **exponentValue
        End Do
      End Do
    End Do
    Do row=1,(order+1)
      exponentValue = row-1
      Do k=1,size(points,1)
        yMatrix(row,1) = 1.0D0*yMatrix(row,1)+1.0D0*points(k,2)*&
        points(k,1)**exponentValue
      End Do
    End Do
! invert xMatrix
    xMatrix = InvertMatrixQ(xMatrix)
! multiply inverse by y to get coefficients
    coefficients = MatMultQ(xMatrix,yMatrix)
! move values to output matrix
    Do row=1,(order+1)
      coefficientsOut(row) = Dble(coefficients(row,1))
    End Do
! Step 2 - extended fit
! If(extendedFit.ne.0)Then
! If required, use superposition fit and check if better than vandermonde
! vRSS = CalcResidualSquareSum(points,coefficients)
! coefficientsSP = SP_PolyFit(points,order,50)
! spRSS = CalcResidualSquareSum(points,coefficientsSP)
! If(spRSS.lt.vRSS)Then
!  coefficients = coefficientsSP
! End If
! End If
  End Function PolyFitQ

  Function PolyFitExact(points) RESULT (coefficients)
! Exactly fit polynomial to points, matrix method
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: yMatrix
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(points,1)) :: xMatrix
    Integer(kind=StandardInteger) :: i, j
! Initialise variables
    coefficients = 0.0D0
    xMatrix = 0.0D0
    yMatrix = 0.0D0
! Make xMatrix and yMatrix
    Do i=1,size(points,1)
      Do j=1,size(points,1)
        xMatrix(i,j) = points(i,1)**(j-1)
      End Do
      yMatrix(i) = points(i,2)
    End Do
    xMatrix = InvertMatrix(xMatrix)
    coefficients = MatMul(xMatrix,yMatrix)
  End Function PolyFitExact

  Function CalcResidualSquareSum(points,coefficients) RESULT (rss)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension( : ) :: coefficients
    Real(kind=DoubleReal) :: rss, x, y
    rss = 0.0D0
    Do i=1,size(points,1)
      x = 1.0D0*points(i,1)
      y = CalcPolynomial(coefficients,x)
      rss = rss + (y-points(i,2))**2
    End Do
  End Function CalcResidualSquareSum

  Function MinPolyFit(points,order) RESULT (x)
! Fit poly to points, then calculate minimum of the curve, assuming it's in the region of the points (+- 25%)
! and there isn't too much wobble
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Integer(kind=StandardInteger) :: order           ! Largest poly term e.g. x3+x2+x+1 3rd order
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients
    Real(kind=DoubleReal) :: xMin, xMax, x
    Integer(kind=StandardInteger) :: i
! Init values
    x = 0.0D0
    Do i=1,size(points,1)
      If(i.eq.1)Then
        xMin = points(i,1)
        xMax = points(i,1)
      Else
        If(points(i,1).lt.xMin)Then
          xMin = points(i,1)
        End If
        If(points(i,1).gt.xMax)Then
          xMax = points(i,1)
        End If
      End If
    End Do
    coefficients = PolyFit(points,order)
    x = MinimaPolynomial (coefficients, xMin, xMax)
  End Function MinPolyFit

  Function MurnFit(points, varianceIn, loopsIn, refinementsIn) RESULT (coefficients)
! Fit Murnaghan EoS to data
! Fitting method adapted from http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
! Murnaghan equation from Murnaghan 1944 described by Fu 1983
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, pointToVary
! Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS, loopFactor
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsQ
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), optional :: varianceIn
    Integer(kind=StandardInteger), optional :: loopsIn, refinementsIn
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops, refinements
! Optional argument variables
    variance = 0.01D0
    loops = 1000
    refinements = 5
    If(present(varianceIn))Then
      variance = varianceIn
    End If
    If(present(loopsIn))Then
      loops = loopsIn
    End If
    If(present(refinementsIn))Then
      refinements = refinementsIn
    End If
! Quadratic fit  (as a starting point)
    coefficientsQ = PolyFit(points,2)
! Random number
    Call RANDOM_NUMBER(randDouble)
! Starting values for fit
    coefficients(2) = (-1.0D0*coefficientsQ(2))/(2.0D0*coefficientsQ(3))    !V0
    coefficients(1) = coefficientsQ(3)*coefficients(2)**2+&                 !E0
    coefficientsQ(2)*coefficients(2)+&
    coefficientsQ(1)
    coefficients(3) = 2.0D0 * coefficientsQ(3) * coefficients(2)            !B0
    coefficients(4) = 2.0D0 + 2.0D0 * randDouble                            !B'0
! Starting RSS
    optRSS = MurnRSS(points,coefficients)
! Adjust points
    Do n=0,refinements
      loopFactor = exp(-0.5D0*n)
      Do i=1,loops
        pointToVary = mod(i-1,4)+1
        coefficientsTemp = coefficients
        Call RANDOM_NUMBER(randDouble)
        varyAmount = variance*2.0D0*(-0.5D0+randDouble)*loopFactor
        coefficientsTemp(pointToVary) = &
        (1.0D0 + varyAmount)*coefficientsTemp(pointToVary)
        testRSS = MurnRSS(points,coefficientsTemp)
        If(testRSS.lt.optRSS)Then
          optRSS = testRSS
          coefficients = coefficientsTemp
          If(optRSS.lt.1.0D-5)Then
            Exit
          End If
        End If
      End Do
      If(optRSS.lt.1.0D-5)Then
        Exit
      End If
    End Do
  End Function MurnFit

  Function MurnCalc(volume,coefficients) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: volume, energy
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Calculate energy
    energy = coefficients(1) + &
    ((coefficients(3)*volume)/coefficients(4))*&
    (((coefficients(2)/volume)**coefficients(4))/&
    (coefficients(4)-1)+1.0D0)-&
    ((coefficients(2)*coefficients(3))/(coefficients(4)-1.0D0))
  End Function MurnCalc

  Function MurnRSS(points,coefficients) RESULT (rss)
! Fit Murnaghan EoS to data
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
    Real(kind=DoubleReal) :: volume, energy, energyC, rss
! calculate RSS
    rss = 0.0D0
    Do i=1,size(points,1)
      volume = points(i,1)
      energy = points(i,2)
      energyC = MurnCalc(volume,coefficients)
      rss = rss + (energyC-energy)**2
    End Do
  End Function MurnRSS

  Function MurnFitBP(points, coefficientsIn) RESULT (coefficients)
! Fits B'0 and holds other coefficients
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
! Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsIn  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops
! Optional argument variables
    variance = 0.001D0
    loops = 100
    coefficients = coefficientsIn
! Starting RSS
    optRSS = MurnRSS(points,coefficients)
! Adjust points
    Do i=1,loops
      coefficientsTemp = coefficients
      Call RANDOM_NUMBER(randDouble)
      varyAmount = variance*2.0D0*(-0.5D0+randDouble)
      coefficientsTemp(4) = &
      (1.0D0 + varyAmount)*coefficientsTemp(4)
      testRSS = MurnRSS(points,coefficientsTemp)
      If(testRSS.lt.optRSS)Then
        optRSS = testRSS
        coefficients = coefficientsTemp
        If(optRSS.lt.1.0D-5)Then
          Exit
        End If
      End If
    End Do
  End Function MurnFitBP

  Function BirchMurnFit(points, varianceIn, loopsIn, refinementsIn, coeffsIn) RESULT (coefficients)
! Fit Murnaghan EoS to data
! Fitting method adapted from http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
! Birch-Murnaghan equation described by Hebbache 2004
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, pointToVary
! Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS, loopFactor
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsQ
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), optional :: varianceIn
    Integer(kind=StandardInteger), optional :: loopsIn, refinementsIn
    Real(kind=DoubleReal), Dimension(1:4), optional :: coeffsIn
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops, refinements
! Optional argument variables
    variance = 0.01D0
    loops = 1000
    refinements = 5
    If(present(varianceIn))Then
      variance = varianceIn
    End If
    If(present(loopsIn))Then
      loops = loopsIn
    End If
    If(present(refinementsIn))Then
      refinements = refinementsIn
    End If
! Init variables
    loopFactor = 1.0D0
! Starting point
    If(Present(coeffsIn))Then
      coefficients = coeffsIn
    Else
! Quadratic fit  (as a starting point)
      coefficientsQ = PolyFit(points,2)
! Random number
      Call RANDOM_NUMBER(randDouble)
! Starting values for fit
      coefficients(2) = (-1.0D0*coefficientsQ(2))/(2.0D0*coefficientsQ(3))    !V0
      coefficients(1) = coefficientsQ(3)*coefficients(2)**2+&                 !E0
      coefficientsQ(2)*coefficients(2)+&
      coefficientsQ(1)
      coefficients(3) = 2.0D0 * coefficientsQ(3) * coefficients(2)            !B0
      coefficients(4) = 2.0D0 + 2.0D0 * randDouble                            !B'0
    End If
! Starting RSS
    optRSS = BirchMurnRSS(points,coefficients)
! Adjust points - second pass
    Do n=0,refinements
      loopFactor = exp(-0.5D0*n)
      Do i=1,loops
        If(n.lt.3)Then
          pointToVary = mod(i-1,3)+1
        Else
          pointToVary = mod(i-1,4)+1
        End If
        coefficientsTemp = coefficients
        Call RANDOM_NUMBER(randDouble)
        varyAmount = variance*2.0D0*(-0.5D0+randDouble)*loopFactor
        coefficientsTemp(pointToVary) = &
        (1.0D0 + varyAmount)*coefficientsTemp(pointToVary)
        testRSS = BirchMurnRSS(points,coefficientsTemp)
        If(testRSS.lt.optRSS)Then
          optRSS = testRSS
          coefficients = coefficientsTemp
          If(optRSS.lt.1.0D-7)Then
            Exit
          End If
        End If
      End Do
      If(optRSS.lt.1.0D-7)Then
        Exit
      End If
    End Do
  End Function BirchMurnFit
  
  Function BirchMurnFit_R(points, varianceIn, loopsIn, refinementsIn, coeffsIn) RESULT (coefficients)
! Fit Murnaghan EoS to data
! Fitting method adapted from http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
! Birch-Murnaghan equation described by Hebbache 2004
! Uses linear congruential generator with same start seed for repeatable fitting
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, pointToVary
! Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS, loopFactor
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsQ
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), optional :: varianceIn
    Integer(kind=StandardInteger), optional :: loopsIn, refinementsIn
    Real(kind=DoubleReal), Dimension(1:4), optional :: coeffsIn
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops, refinements
! Set seed
    randDouble = RandomLCG(12791244)
! Optional argument variables
    variance = 0.01D0
    loops = 1000
    refinements = 5
    If(present(varianceIn))Then
      variance = varianceIn
    End If
    If(present(loopsIn))Then
      loops = loopsIn
    End If
    If(present(refinementsIn))Then
      refinements = refinementsIn
    End If
! Init variables
    loopFactor = 1.0D0
! Starting point
    If(Present(coeffsIn))Then
      coefficients = coeffsIn
    Else
! Quadratic fit  (as a starting point)
      coefficientsQ = PolyFit(points,2)
! Random number
      randDouble = RandomLCG()
! Starting values for fit
      coefficients(2) = (-1.0D0*coefficientsQ(2))/(2.0D0*coefficientsQ(3))    !V0
      coefficients(1) = coefficientsQ(3)*coefficients(2)**2+&                 !E0
      coefficientsQ(2)*coefficients(2)+&
      coefficientsQ(1)
      coefficients(3) = 2.0D0 * coefficientsQ(3) * coefficients(2)            !B0
      coefficients(4) = 2.0D0 + 2.0D0 * randDouble                            !B'0
    End If
! Starting RSS
    optRSS = BirchMurnRSS(points,coefficients)
! Adjust points - second pass
    Do n=0,refinements
      loopFactor = exp(-0.5D0*n)
      Do i=1,loops
        If(n.lt.3)Then
          pointToVary = mod(i-1,3)+1
        Else
          pointToVary = mod(i-1,4)+1
        End If
        coefficientsTemp = coefficients
        randDouble = RandomLCG()
        varyAmount = variance*2.0D0*(-0.5D0+randDouble)*loopFactor
        coefficientsTemp(pointToVary) = &
        (1.0D0 + varyAmount)*coefficientsTemp(pointToVary)
        testRSS = BirchMurnRSS(points,coefficientsTemp)
        If(testRSS.lt.optRSS)Then
          optRSS = testRSS
          coefficients = coefficientsTemp
          If(optRSS.lt.1.0D-7)Then
            Exit
          End If
        End If
      End Do
      If(optRSS.lt.1.0D-7)Then
        Exit
      End If
    End Do
  End Function BirchMurnFit_R

  Function BirchMurnCalc(volume,coefficients) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: volume, energy, eta
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients    ! E0, V0, B0, B'0
! Calculate energy
    eta = (volume/coefficients(2))**(1.0D0/3.0D0)
    energy = coefficients(1) + &
    ((9.0D0*coefficients(3)*coefficients(2))/(16.0D0))*&
    ((eta**2-1.0D0)**2)*&
    (6.0D0+coefficients(4)*(eta**2-1.0D0)-4.0D0*eta**2)
  End Function BirchMurnCalc

  Function BirchMurnRSS(points,coefficients) RESULT (rss)
! Fit Murnaghan EoS to data
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients    ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: volume, energy, energyC, rss
! calculate RSS
    rss = 0.0D0
    Do i=1,size(points,1)
      volume = points(i,1)
      energy = points(i,2)
      energyC = BirchMurnCalc(volume,coefficients)
      rss = rss + (energyC-energy)**2
    End Do
  End Function BirchMurnRSS

  Function BirchMurnFitBP(points, coefficientsIn) RESULT (coefficients)
! Fits B'0 and holds other coefficients
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
! Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsIn  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops
! Optional argument variables
    variance = 0.02D0
    loops = 100
    coefficients = coefficientsIn
! Starting RSS
    optRSS = BirchMurnRSS(points,coefficients)
! Adjust points
    Do i=1,loops
      coefficientsTemp = coefficients
      Call RANDOM_NUMBER(randDouble)
      varyAmount = variance*2.0D0*(-0.5D0+randDouble)
      coefficientsTemp(4) = &
      (1.0D0 + varyAmount)*coefficientsTemp(4)
      testRSS = BirchMurnRSS(points,coefficientsTemp)
      If(testRSS.lt.optRSS)Then
        optRSS = testRSS
        coefficients = coefficientsTemp
        If(optRSS.lt.1.0D-5)Then
          Exit
        End If
      End If
    End Do
  End Function BirchMurnFitBP

! Interpolation


  Function InterpMatrix(points) RESULT (coefficients)
! Calculates coefficients of interpolation polynomial
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, matrixSize
    Real(kind=DoubleReal), Dimension( : , : ) :: points  
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients 
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(points,1)) :: xMatrix 
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: yMatrix 
! Init
    matrixSize = size(points,1)  
    coefficients = 0.0D0    
    xMatrix = 0.0D0    
    yMatrix = 0.0D0    
! Build Y and X matrix  
    Do i=1,matrixSize
      Do j=1,matrixSize
        xMatrix(i,j) = 1.0D0*points(i,1)**(j-1)
      End Do
      yMatrix(i) = 1.0D0*points(i,2)
    End Do
! Invert xMatrix (reuse xMatrix)
    xMatrix = InvertMatrix(xMatrix)
! calculate coefficients c = x^(-1)y
    coefficients = matmul(xMatrix,yMatrix)
  End Function InterpMatrix
  
  Function InterpMatrixPoint(x, points, derivativeIn) RESULT (y)
! Calculates values f(x), f'(x), f''(x) of a point using polynomial interpolation
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:size(points,2)) :: pointsTemp
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients 
    Integer(kind=StandardInteger) :: i, derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal) :: x, y
! Initialise variables
    y = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! Get coefficients
    coefficients = InterpMatrix(points)
! x,f(x)
    If(derivative.lt.2)Then
! Calculate value of point    
      y = CalcPolynomial(coefficients, x, derivative)
    Else
! Original polynomial not accurate enough - update points from x,f(x) to x,f'(x)    
      Do i=1,size(points,1)
        pointsTemp(i,1) = points(i,1)
        pointsTemp(i,2) = CalcPolynomial(coefficients, points(i,1), 1)
      End Do
! Get coefficients
      coefficients = InterpMatrix(pointsTemp)   
! Calculate value of point    
      y = CalcPolynomial(coefficients, x, 1)      
    End If
  End Function InterpMatrixPoint
    
  Function InterpMatrixPointInaccurate(x, points, derivativeIn) RESULT (y)
! Calculates values f(x), f'(x), f''(x) of a point using polynomial interpolation
! Purposeley inaccurate function
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients 
    Integer(kind=StandardInteger) :: derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal) :: x, y
! Initialise variables
    y = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! Get coefficients
    coefficients = InterpMatrix(points)   
! Calculate value of point    
    y = CalcPolynomial(coefficients, x, derivative)      
  End Function InterpMatrixPointInaccurate
  

  Function InterpLagrange(x, points, derivativeIn) RESULT (output)
! Calculates y(x), y'(x) or y''(x) using Lagrange interpolation
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension( : , : ) :: points
    Real(kind=DoubleReal), Dimension(1:size(points,1),1:2) :: pointsTemp
    Real(kind=DoubleReal), Dimension(1:size(points,1)) :: coefficients
    Integer(kind=StandardInteger) :: j, i, n, k, derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal) :: numerator, denominator, numeratorPart
    Real(kind=DoubleReal) :: x, y, dy, output, xTemp, dyy
! Initialise variables
    output = 0.0D0
! Handle optional argument
    derivative = 0
    If(Present(derivativeIn))Then
      derivative = derivativeIn
    End If
! y(x)
    If(derivative.eq.0)Then
! Make coefficients
      Do n=1,size(points,1)
        numerator = 1.0D0
        denominator = 1.0D0
        Do k=1,size(points,1)
          If(k.ne.n)Then
            numerator=numerator*(x-points(k,1))
            denominator=denominator*(points(n,1)-points(k,1))
          End If
        End Do
        coefficients(n)=1.0D0*(numerator/denominator)
      End Do
! Calculate y
      y = 0.0D0
      Do n=1,size(points,1)
        y=y+points(n,2)*coefficients(n)
      End Do
      output = y
    End If
! y'(x)
    If(derivative.eq.1)Then
      Do n=1,size(points,1)
        numerator = 0.0D0
        denominator = 1.0D0
        Do k=1,size(points,1)
          If(k.ne.n)Then
            denominator=denominator*(points(n,1)-points(k,1))
            numeratorPart = 1.0D0
            Do i=1,size(points,1)
              If(i.ne.n.and.i.ne.k)Then
                numeratorPart=numeratorPart*(x-points(i,1))
              End If
            End Do
            numerator=numerator+numeratorPart
          End If
        End Do
        coefficients(n)=1.0D0*(numerator/denominator)
      End Do
! Calculate dy
      dy = 0.0D0
      Do n=1,size(points,1)
        dy = dy + points(n,2)*coefficients(n)
      End Do
      output = dy
    End If
! y''(x)
    If(derivative.eq.2)Then
! Could use recursive functions for higher orders, but y''(x) high enough for now
! Calculate y'(x) for each x input
      Do j = 1,size(points,1)
        xTemp = points(j,1)
        Do n=1,size(points,1)
          numerator = 0.0D0
          denominator = 1.0D0
          Do k=1,size(points,1)
            If(k.ne.n)Then
              denominator=denominator*(points(n,1)-points(k,1))
              numeratorPart = 1.0D0
              Do i=1,size(points,1)
                If(i.ne.n.and.i.ne.k)Then
                  numeratorPart=numeratorPart*(xTemp-points(i,1))
                End If
              End Do
              numerator=numerator+numeratorPart
            End If
          End Do
          coefficients(n)=1.0D0*(numerator/denominator)
        End Do
! Calculate dy
        dy = 0.0D0
        Do n=1,size(points,1)
          dy = dy + points(n,2)*coefficients(n)
        End Do
        pointsTemp(j,1) = xTemp
        pointsTemp(j,2) = dy
      End Do
! Use the x, y'(x) points to interpolate y''(x)
      Do n=1,size(points,1)
        numerator = 0.0D0
        denominator = 1.0D0
        Do k=1,size(points,1)
          If(k.ne.n)Then
            denominator=denominator*(pointsTemp(n,1)-pointsTemp(k,1))
            numeratorPart = 1.0D0
            Do i=1,size(pointsTemp,1)
              If(i.ne.n.and.i.ne.k)Then
                numeratorPart=numeratorPart*(x-pointsTemp(i,1))
              End If
            End Do
            numerator=numerator+numeratorPart
          End If
        End Do
        coefficients(n)=1.0D0*(numerator/denominator)
      End Do
! Calculate dy
      dyy = 0.0D0
      Do n=1,size(points,1)
        dyy = dyy + pointsTemp(n,2)*coefficients(n)
      End Do
      output = dyy
    End If
  End Function InterpLagrange
    
  Function PointInterp(points,x,subsetSize,derivativeIn,inputSetStartIn,inputSetLengthIn) RESULT (yArray)
! Takes large set of data points, finds region of points around the input "x", and interps with lagrange
    Implicit None  !Force declaration of all variables
! Declare variables - In/Out
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger), optional :: derivativeIn, inputSetStartIn, inputSetLengthIn
    Integer(kind=StandardInteger) :: subsetSize, inputStart, inputLength, derivative
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Declare variables
    Real(kind=DoubleReal), Dimension(1:subsetSize,1:2) :: pointsInterp
    Real(kind=DoubleReal) :: xLower, xUpper
    Integer(kind=StandardInteger) :: i, j, dataSetSize, xPos
    Integer(kind=StandardInteger) :: xPosOffset, xPosOffsetR, xPosStart
    Integer(kind=StandardInteger) :: inputEnd, xPosUpper, xPosLower
! Initialise Variables
    xPos = 1
    dataSetSize = size(points,1)
    yArray = 0.0D0
! Handle optional arguments
    inputStart = 1
    inputLength = dataSetSize
    derivative = 0
    If(Present(inputSetStartIn))Then
      inputStart = inputSetStartIn
    End If
    If(Present(inputSetLengthIn))Then
      inputLength = inputSetLengthIn
    End If
    inputEnd = inputStart+inputLength-1
    If(Present(derivativeIn))Then
      derivative = derivativeIn
      If(derivative.lt.0)Then
        derivative = 0
      End If
      If(derivative.gt.2)Then
        derivative = 2
      End If
    End If
! Check data set size
    If(subsetSize.lt.2)Then
      subsetSize = 2
    End If
    If(subsetSize.gt.inputLength)Then
      subsetSize = inputLength
    End If
! Reduce set of data points
    If(subsetSize.eq.inputLength)Then
      j = 0
      Do i=inputStart,inputEnd
        j = j + 1
        pointsInterp(j,1) = points(i,1)
        pointsInterp(j,2) = points(i,2)
      End Do
    ElseIf(subsetSize.lt.inputLength)Then
! Reduce set of data points
      xLower = points(inputStart,1)
      xUpper = points(inputEnd,1)
! Find xPos
      If(x.lt.xLower)Then  !If x lower than data set, use lowest possible points
        xPos = inputStart
      ElseIf(x.gt.xUpper)Then  !If x higher than data set, use highest possible points
        xPos = inputEnd
      Else
! Estimate position
        xPos = INT(Floor(((x - xLower) / (xUpper - xLower)) * 1.0D0 * inputLength) + inputStart)
        If(xPos.lt.inputStart)Then
          xPos = inputStart
        End If
        If((xPos+1).gt.inputEnd)Then
          xPos = inputEnd-1
        End If
        xLower = points(xPos,1)
        xUpper = points(xPos+1,1)
! If estimate is incorrect, search for better value
        If(x.lt.xLower)Then
          xPosStart = xPos
          Do xPos=xPosStart,inputStart,-1    !Search down
            xLower = points(xPos,1)
            xUpper = points(xPos+1,1)
            If(x.le.xUpper.and.x.ge.xLower)Then
              Exit  !xPos found
            End If
          End Do
        End If
        If(x.gt.xUpper)Then
          xPosStart = xPos
          Do xPos=xPosStart,inputEnd,+1    !Search down
            xLower = points(xPos,1)
            xUpper = points(xPos+1,1)
            If(x.le.xUpper.and.x.ge.xLower)Then
              Exit  !xPos found
            End If
          End Do
        End If
      End If
! Adjust xPos to center of subset
      xPosOffset = INT(Floor(1.0D0*subsetSize/2))
      xPosOffsetR = subsetSize - xPosOffset
      xPosLower = xPos - xPosOffset
      xPosUpper = xPos + xPosOffsetR - 1
! Adjust xPos start, so it fits in the range of subset of selected data points
      If(xPosLower.lt.inputStart)Then
        xPosLower = inputStart
        xPosUpper = inputStart + subsetSize - 1
      End If
      If(xPosUpper.gt.inputEnd)Then
        xPosLower = inputEnd - subsetSize + 1
        xPosUpper = inputEnd
      End If
! Transfer data points to pointsInterp
      j = 0
      Do i=xPosLower,xPosUpper
        j = j + 1
        pointsInterp(j,1) = points(i,1)
        pointsInterp(j,2) = points(i,2)
      End Do
    End If
! Store interpolation results
    If(derivative.ge.0)Then
      yArray(1) = InterpLagrange(x, pointsInterp)
    End If
    If(derivative.ge.1)Then
      yArray(2) = InterpLagrange(x, pointsInterp, 1)
    End If
    If(derivative.ge.2)Then
      yArray(3) = InterpLagrange(x, pointsInterp, 2)
    End If
  End Function PointInterp
  
  Function PolyPoints(coefficients,xStart,xEnd,points) RESULT (polyPointsArr)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, points
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal), Dimension(1:points,1:2) :: polyPointsArr
    Real(kind=DoubleReal) :: xStart, xEnd, xInc, xVal, yVal
! Set values
    xInc = (xEnd-xStart)/(1.0D0*(points-1.0D0))
! Loop through points
    Do i=1,points
      xVal = xStart+1.0D0*(i-1)*xInc
      yVal = CalcPolynomial(coefficients, xVal)
      polyPointsArr(i,1) = xVal
      polyPointsArr(i,2) = yVal
    End Do    
  End Function PolyPoints  

  Function BirchMurnPoints(coefficients,xStart,xEnd,points) RESULT (bmPointsArr)
! Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
    Real(kind=DoubleReal), Dimension(1:points,1:2) :: bmPointsArr
    Real(kind=DoubleReal) :: xStart, xEnd, xInc, xVal, yVal
! Set values
    xInc = (xEnd-xStart)/(1.0D0*(points-1.0D0))
! Loop through points
    Do i=1,points
      xVal = xStart+(i-1)*xInc
      yVal = BirchMurnCalc(xVal, coefficients)
      bmPointsArr(i,1) = xVal
      bmPointsArr(i,2) = yVal
    End Do    
  End Function BirchMurnPoints    
  
  Function SingleDecayFit(dataPoints) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, iA, m, n, maxLoops
    Real(kind=DoubleReal) :: rss, lastRSS, optRSS, convergence, maxRSSVal
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:2) :: J
    Real(kind=DoubleReal), Dimension(1:2) :: parameters, parametersOpt, change
    Real(kind=DoubleReal), Dimension(1:3) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: gridA
    Real(kind=DoubleReal) :: a_T, lA_T
    Real(kind=DoubleReal), Dimension(1:100,1:3) :: topBoxes 
    Integer(kind=StandardInteger) :: topBoxCount, maxRSS 
    Logical :: storeFlag
!--------------------------------------------------
! Find Starting Parameters  
!--------------------------------------------------
! Init
    topBoxCount = 4
    topBoxes = 2.0D20
! Set a ranges
    gridA = 12
    Do i=1,gridA
      aRange(i,1) = -1.0D0*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do n=1,1000
        Do m=1,3
          a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
          lA_T = 40.0D0*(RandomLCG()-0.5D0)
! Calc rss
          rss = SingleDecayFitRSS(dataPoints, a_T, lA_T)  
! Check if better than boxed - update if better than already stored for this box         
          storeFlag = .true.
          Do i=1,topBoxCount
            If(topBoxes(i,2).eq.a_T.and.rss.lt.topBoxes(i,1))Then
              topBoxes(i,1) = rss
              topBoxes(i,2) = a_T
              topBoxes(i,3) = lA_T
              storeFlag = .false.
              Exit
            End If  
            If(i.eq.1)Then
              maxRSS = 1
              maxRSSVal = topBoxes(1,1)
            Else
              If(topBoxes(i,1).gt.maxRSSVal)Then
                maxRSS = i
                maxRSSVal = topBoxes(i,1)
              End If            
            End If
          End Do
! If better than any in box
          If(storeFlag)Then
            If(rss.lt.maxRSSVal)Then
              topBoxes(maxRSS,1) = rss
              topBoxes(maxRSS,2) = a_T
              topBoxes(maxRSS,3) = lA_T
            End If
          End If          
        End Do
      End Do
    End Do
!--------------------------------------------------
! Newton Gauss Elimination
!--------------------------------------------------
    Do n=1,topBoxCount
      parameters(1) = topBoxes(n,2)
      parameters(2) = topBoxes(n,3) 
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      maxLoops = 0
      Do While(maxLoops.le.100.and.convergence.gt.1.0D-7)
        maxLoops = maxLoops + 1
        rss = 0.0D0
        Do i=1,size(dataPoints,1)
          R(i) = parameters(1)*exp(parameters(2)*dataPoints(i,1))-dataPoints(i,2)   ! f(x)-y
          J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          rss = rss + R(i)**2
        End Do
        change = NewtonGaussOpt(J,R)
        Do i=1,size(change)
          parameters(i) = parameters(i) + change(i)
        End Do
        convergence = abs(lastRSS-rss)
        lastRSS = rss
      End Do
      If(n.eq.1)Then
        optRSS = rss
        Do i=1,size(change)
          parametersOpt(i) = parameters(i)
        End Do
      Else
        If(rss.lt.optRSS)Then
          optRSS = rss
          Do i=1,size(change)
            parametersOpt(i) = parameters(i)
          End Do
        End If
      End If  
    End Do  
    output(1) = parametersOpt(1)
    output(2) = parametersOpt(2)
    output(3) = optRSS
  End Function SingleDecayFit    
    
  Function DoubleDecayFit(dataPoints,gridAIn,gridAFactorIn,&
  searchLoopsIn) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA) + b exp(lB)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, m, k
    Integer(kind=StandardInteger) :: iA,iB
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, lA, lB
    Real(kind=DoubleReal) :: a_T, b_T, lA_T, lB_T
    Real(kind=DoubleReal) :: rss, startRSS, optRSS, lastRSS, convergence, lambda  
    Real(kind=DoubleReal), Dimension(1:6) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: lRange
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: searchLoops, gridA
    Integer(kind=StandardInteger), Optional :: searchLoopsIn, gridAIn
    Real(kind=DoubleReal) :: gridAFactor
    Real(kind=DoubleReal), Optional :: gridAFactorIn
    Integer(kind=StandardInteger) :: maxRSS
    Real(kind=DoubleReal) :: maxRSSVal
    Integer(kind=StandardInteger) :: topBoxCount
    Real(kind=DoubleReal), Dimension(1:100,1:5) :: topBoxes 
    Logical :: storeFlag
! Newton Gauss Opt Vars
    Real(kind=DoubleReal), Dimension(1:4) :: parameters, parameters_Last, parametersOpt  
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R, R_Last
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:4) :: J, J_Last
    Real(kind=DoubleReal), Dimension(1:4,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:4) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:4) :: P      ! Change
! User options:
! gridA - size of grid -10^(grid-6) to 10^(grid-6) [grid size = gridA^2]
! gridAFactor - span of grid increased/decreased by factor, bt keeps grid size same
! searchLoops - no. random lA and lB points for each (a,b) box
! vary loops - number of refinement loops for final result
! filter loops - number of loops to filter out top results   
!     
! Init
    ddfRssCount = 0
! Optional   
    gridA = 10
    If(Present(gridAIn))Then
      gridA = gridAIn
    End If  
    gridAFactor = 1.0D0
    If(Present(gridAFactorIn))Then
      gridAFactor = gridAFactorIn
    End If  
    searchLoops = 50
    If(Present(searchLoopsIn))Then
      searchLoops = searchLoopsIn
    End If      
! Init
    topBoxCount = 5
    topBoxes = 2.0D20
! Set first test values
    a = 1.0D0
    b = 1.0D0
    lA = 1.0D0
    lB = 1.0D0
    startRSS = DoubleDecayFitRSS(dataPoints, a, b, lA, lB)
    optRSS = startRSS
! Assume:
! -10^5 < a < 10^5 
! -10^5 < b < 10^5 
! -20 < lA < 20 
! -20 < lB < 20 
    Do i=1,5
      lRange(i,1) = -2.0D0*10D0**(2-i)
      lRange(i,2) = -2.0D0*10D0**(1-i)
    End Do    
    Do i=1,5
      lRange(i+5,1) = 2.0D0*10D0**(i-5)
      lRange(i+5,2) = 2.0D0*10D0**(i-4)
    End Do
    lRange(5,2) = 0.0D0
    lRange(6,1) = 0.0D0    
! Set a and b ranges
    Do i=1,gridA
      aRange(i,1) = -1.0D0*gridAFactor*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*gridAFactor*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*gridAFactor*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*gridAFactor*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do iB=iA,(gridA+gridA)  
! Set values
        Do m=1,3
          Do n=1,searchLoops
            a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
            lA_T = 20.0D0*(RandomLCG()-0.5D0)
            b_T = (0.25D0*m)*(aRange(iB,1)+aRange(iB,2))
            lB_T = 20.0D0*(RandomLCG()-0.5D0)
! Calculate RSS
            rss = DoubleDecayFitRSS(dataPoints, a_T, b_T, lA_T, lB_T)
            storeFlag = .true.
            Do i=1,topBoxCount
              If(topBoxes(i,2).eq.a_T.and.topBoxes(i,4).eq.b_T)Then
                storeFlag = .false.
                If(rss.lt.topBoxes(i,1))Then
                  topBoxes(i,1) = rss
                  topBoxes(i,2) = a_T
                  topBoxes(i,3) = lA_T
                  topBoxes(i,4) = b_T
                  topBoxes(i,5) = lB_T
                End If
                Exit
              End If
            End Do
! Store if in best 10
            If(storeFlag)Then
              maxRSS = 1
              maxRSSVal = topBoxes(1,1)
              Do i=2,topBoxCount
                If(topBoxes(i,1).gt.maxRSSVal)Then
                  maxRSS = i
                  maxRSSVal = topBoxes(i,1)
                End If
              End Do
              If(rss.lt.maxRSSVal)Then
                topBoxes(maxRSS,1) = rss
                topBoxes(maxRSS,2) = a_T
                topBoxes(maxRSS,3) = lA_T
                topBoxes(maxRSS,4) = b_T
                topBoxes(maxRSS,5) = lB_T
              End If
            End If
          End Do
        End Do
      End Do
    End Do  
!--------------------------------------------------
! LMA
!--------------------------------------------------    
    Do k=1,topBoxCount
      parameters(1) = topBoxes(k,2)
      parameters(2) = topBoxes(k,3)
      parameters(3) = topBoxes(k,4)
      parameters(4) = topBoxes(k,5)
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      lambda = 1.0D0
      Do n=1,100
        rss = 0.0D0
  ! Make Jacobian and Residuals matrix
        Do i=1,size(dataPoints,1)
          R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
          parameters(3)*exp(parameters(4)*dataPoints(i,1)))-dataPoints(i,2)   ! f(x)-y
          J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          J(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
          J(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4
          rss = rss + R(i)**2
        End Do
  ! Choose whether to accept update or increase/decrease lambda      
        If(n.gt.1)Then
  ! Delayed gratification scheme - 1.5*lambda or 0.2*lambda     
          If(rss.gt.lastRSS)Then  ! If worse...reject, and increase lambda
  ! Discard changes, increase lambda
            J = J_Last
            R = R_Last
            parameters = parameters_Last
            rss = lastRSS
            lastRSS = -1.0D0
            lambda = lambda * 1.5D0
          End If
          If(rss.lt.lastRSS)Then  ! If better...accept, and decrease lambda
            lambda = lambda * 0.2D0
          End If
        End If
! calculate change matrix
        !***********     
        ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
        !***********      
  ! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
        JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)  
  ! convergence of RSS
        convergence = abs(lastRSS-rss)  
  ! Store last loop values
        parameters_Last = parameters
        lastRSS = rss 
        J_Last = J
        R_Last = R
  ! Update parameters      
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do   
  ! Breakout if convergence threshold met      
        If(convergence.lt.1.0D-7)Then
          Exit
        End If
      End Do
      If(rss.lt.optRSS)Then
        optRSS = rss
        Do i=1,size(P)
          parametersOpt(i) = parameters(i)
        End Do  
      End If
  ! Breakout if convergence threshold met      
      If(convergence.lt.1.0D-7)Then
        Exit
      End If
    End Do    
! Output
    output(1) = parametersOpt(1)
    output(2) = parametersOpt(2)
    output(3) = parametersOpt(3)
    output(4) = parametersOpt(4)
    output(5) = optRSS
    output(6) = ddfRssCount
  End Function DoubleDecayFit   
  
  Function DoubleDecayFitRSS(dataPoints, a, b, lA, lB) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, b, lA, lB, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)+b*exp(lB*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
    ddfRssCount = ddfRssCount + 1
  End Function DoubleDecayFitRSS 
  
  Function SingleDecayFitRSS(dataPoints, a, lA) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, lA, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function SingleDecayFitRSS 
  
  
! ------------------------------------------------------------------------!
! Optimization Functions
! ------------------------------------------------------------------------!  
  
  Function NewtonGaussOpt(J,R) RESULT (P)!   
    Implicit None  !Force declaration of all variables
! Declare variables 
    Real(kind=DoubleReal), Dimension(:,:) :: J   ! Jacobian
    Real(kind=DoubleReal), Dimension(:) :: R      ! Residuals
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:size(J,2),1:size(J,2)) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:size(J,2)) :: P      ! Change
!***********     
! P = (JTJ)^(-1)(-1*JTR)   
!***********      
! Transpose Jacobian
    JT = TransposeMatrix(J)
    JTJ = matmul(JT,J)
    JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)
    JTR = matmul(JT,R)
    JTR = -1.0D0*JTR ! Recycle JTR var
    P = matmul(JTJ,JTR)  
  End Function NewtonGaussOpt
  

  
! ------------------------------------------------------------------------!
! Spline Functions
! ------------------------------------------------------------------------!

  Function SplineAB(pointA, pointB) RESULT (coefficients)
! Polynomial to spline between points A and B - x, f(x), f'(x) and f''(x) supplied
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:) :: pointA  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
    Real(kind=DoubleReal), Dimension(:) :: pointB  !1 x, 2 f(x), 3 f'(x), 4 f''(x)....
    Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1))) :: coefficients, yMatrix
    Real(kind=DoubleReal), Dimension(1:(2*(size(pointA,1)-1)),1:(2*(size(pointA,1)-1))) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(size(pointA,1)-1),1:(2*(size(pointA,1)-1))) :: &
    xMatrixCoeffs, xMatrixExp
    Integer(kind=StandardInteger) :: row, col, matrixSize, matrixHalfSize
! Init variables
    coefficients = 0.0D0
    yMatrix = 0.0D0
    xMatrix = 0.0D0
    matrixSize = 2*(size(pointA,1)-1)
    matrixHalfSize = (size(pointA,1)-1)
! Make xMatrixExp
    Do row=1,matrixHalfSize
      Do col=1,matrixSize
        xMatrixExp(row,col) = col-row
        If(xMatrixExp(row,col).lt.0)Then
          xMatrixExp(row,col) = 0
        End If
      End Do
    End Do
! Make xMatrixCoeffs
    Do col=1,matrixSize
      xMatrixCoeffs(1,col) = 1.0D0
    End Do
    Do row=2,matrixHalfSize
      Do col=1,matrixSize
        xMatrixCoeffs(row,col) = 1.0D0*xMatrixExp(row-1,col)*xMatrixCoeffs(row-1,col)
      End Do
    End Do
! Make xMatrix
    Do row=1,matrixHalfSize
      Do col=1,matrixSize
        xMatrix(row,col) = &
        1.0D0*xMatrixCoeffs(row,col)*pointA(1)**xMatrixExp(row,col)
      End Do
    End Do
    Do row=1,matrixHalfSize
      Do col=1,matrixSize
        xMatrix(row+matrixHalfSize,col) = &
        1.0D0*xMatrixCoeffs(row,col)*pointB(1)**xMatrixExp(row,col)
      End Do
    End Do
! Make yMatrix
    Do row=1,matrixHalfSize
      yMatrix(row) = pointA(row+1)
    End Do
    Do row=1,matrixHalfSize
      yMatrix(row+matrixHalfSize) = pointB(row+1)
    End Do
! Invert xMatrix
    xMatrix = InvertMatrix(xMatrix)
! Solve equation
    coefficients = matmul(xMatrix,yMatrix)
  End Function SplineAB

  Function SplineNodes(inputNodes,numDataPoints,startPointIn,endPointIn) RESULT (dataPoints)
! Input nodes x,f(x) for each node, calculate f'(x) and f''(x) from the set of nodes, then spline
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputNodes
    Real(kind=DoubleReal), Dimension(1:size(inputNodes,1),1:4) :: splinePoints
    Integer(kind=StandardInteger) :: i, numDataPoints, nodeKey
    Real(kind=DoubleReal), Dimension(1:numDataPoints,1:4) :: dataPoints
    Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
    Integer(kind=StandardInteger), Optional :: startPointIn, endPointIn
    Integer(kind=StandardInteger) :: startPoint, endPoint
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Init Variables
    dataPoints = 0.0D0
    startPoint = 1
    endPoint = size(inputNodes,1)
    If(Present(startPointIn))Then
      startPoint = startPointIn
    End If
    If(Present(endPointIn))Then
      endPoint = endPointIn
    End If    
! set f'(x) and f''(x)
    Do i=startPoint,endPoint
      x = inputNodes(i,1)
      yArray = PointInterp(inputNodes,x,3,2,startPoint,endPoint)
      splinePoints(i,1) = x
      splinePoints(i,2) = yArray(1)
      splinePoints(i,3) = yArray(2)
      splinePoints(i,4) = yArray(3)
    End Do
! Calculate spline data points
    xStart = splinePoints(startPoint,1)
    xEnd = splinePoints(endPoint,1)
    xIncrement = (xEnd-xStart)/(1.0D0*numDataPoints-1)
! Loop through data points
    nodeKey = startPoint-1
    x = xStart
    Do i=1,numDataPoints
      If((i.eq.1).or.(x.ge.inputNodes(nodeKey+1,1).and.(nodeKey+1).lt.endPoint))Then
        nodeKey = nodeKey + 1
        pointA(1) = inputNodes(nodeKey,1)
        pointA(2) = inputNodes(nodeKey,2)
        pointA(3) = inputNodes(nodeKey,3)
        pointA(4) = inputNodes(nodeKey,4)
        pointB(1) = inputNodes(nodeKey+1,1)
        pointB(2) = inputNodes(nodeKey+1,2)
        pointB(3) = inputNodes(nodeKey+1,3)
        pointB(4) = inputNodes(nodeKey+1,4)
        coefficients = SplineAB(pointA, pointB)
      End If
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomial (coefficients, x, 0)
      dataPoints(i,3) = CalcPolynomial (coefficients, x, 1)
      dataPoints(i,4) = CalcPolynomial (coefficients, x, 2)
! Increment x
      x = x + xIncrement
    End Do
  End Function SplineNodes
  
  Function SplineNodesV(inputNodes,numDataPoints,startPoint,endPoint,dataSize) RESULT (dataPoints)
! Input nodes x,f(x) for each node, calculate f'(x) and f''(x) from the set of nodes, then spline
! Variable length output
    Implicit None  !Force declaration of all variables
! Declare variables - arg
    Real(kind=DoubleReal), Dimension(:,:) :: inputNodes
    Integer(kind=StandardInteger) :: numDataPoints, startPoint, endPoint, dataSize
! Declare variables - priv
    Real(kind=DoubleReal), Dimension(1:size(inputNodes,1),1:4) :: splinePoints
    Integer(kind=StandardInteger) :: i, nodeKey
    Real(kind=DoubleReal), Dimension(1:dataSize,1:4) :: dataPoints
    Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Init Variables
    dataPoints = 0.0D0
    If(startPoint.eq.0)Then
      startPoint = 1
    End If
    If(endPoint.eq.0)Then
      endPoint = size(inputNodes,1)    
    End If
! set f'(x) and f''(x)
    Do i=startPoint,endPoint
      x = inputNodes(i,1)
      yArray = PointInterp(inputNodes,x,3,2,startPoint,endPoint)
      splinePoints(i,1) = x
      splinePoints(i,2) = yArray(1)
      splinePoints(i,3) = yArray(2)
      splinePoints(i,4) = yArray(3)
    End Do
! Calculate spline data points
    xStart = splinePoints(startPoint,1)
    xEnd = splinePoints(endPoint,1)
    xIncrement = (xEnd-xStart)/(1.0D0*numDataPoints-1)
! Loop through data points
    nodeKey = startPoint-1
    x = xStart
    Do i=1,numDataPoints
      If((i.eq.1).or.(x.ge.inputNodes(nodeKey+1,1).and.(nodeKey+1).lt.endPoint))Then
        nodeKey = nodeKey + 1
        pointA(1) = inputNodes(nodeKey,1)
        pointA(2) = inputNodes(nodeKey,2)
        pointA(3) = inputNodes(nodeKey,3)
        pointA(4) = inputNodes(nodeKey,4)
        pointB(1) = inputNodes(nodeKey+1,1)
        pointB(2) = inputNodes(nodeKey+1,2)
        pointB(3) = inputNodes(nodeKey+1,3)
        pointB(4) = inputNodes(nodeKey+1,4)
        coefficients = SplineAB(pointA, pointB)
      End If
      dataPoints(i,1) = x
      dataPoints(i,2) = CalcPolynomial (coefficients, x, 0)
      dataPoints(i,3) = CalcPolynomial (coefficients, x, 1)
      dataPoints(i,4) = CalcPolynomial (coefficients, x, 2)
! Increment x
      x = x + xIncrement
    End Do
  End Function SplineNodesV

  Function SplineComplete(inputPoints,interpSizeIn) RESULT (splinePoints)
! Complete missing points f'(x) f''(x) by interpolation
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputPoints
    Integer(kind=StandardInteger) :: interpSizeIn
    Real(kind=DoubleReal), Dimension(1:size(inputPoints,1),1:4) :: splinePoints
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: x
    Integer(kind=StandardInteger) :: i
! Init variables
    splinePoints = 0.0D0
    Do i=1,size(inputPoints,1)
      x = inputPoints(i,1)
      yArray = PointInterp(inputPoints,x,interpSizeIn,2)
      splinePoints(i,1) = x
      splinePoints(i,2) = yArray(1)  ! f(x)
      splinePoints(i,3) = yArray(2)  ! f'(x)
      splinePoints(i,4) = yArray(3)  ! f''(x)
    End Do
  End Function SplineComplete

  Function Spline(inputNodes,numDataPoints,startPointIn,endPointIn) RESULT (dataPoints)
! NOT WORKING YET - 12092014
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputNodes
    Integer(kind=StandardInteger) :: i, numDataPoints, nodeKey
    Real(kind=DoubleReal), Dimension(1:numDataPoints,1:4) :: dataPoints
    Real(kind=DoubleReal) :: x, xStart, xEnd, xIncrement
    Integer(kind=StandardInteger), Optional :: startPointIn, endPointIn
    Integer(kind=StandardInteger) :: startPoint, endPoint
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Fill in
!
! Init Variables
    dataPoints = 0.0D0
    startPoint = 1
    endPoint = size(inputNodes,1)
    If(Present(startPointIn))Then
      startPoint = startPointIn
    End If
    If(Present(endPointIn))Then
      endPoint = endPointIn
    End If
    xStart = inputNodes(startPoint,1)
    xEnd = inputNodes(endPoint,1)
    xIncrement = (xEnd-xStart)/(1.0D0*numDataPoints-1)
! Loop through data points
    nodeKey = startPoint-1
    x = xStart
    Do i=1,numDataPoints
      If((i.eq.1).or.(x.ge.inputNodes(nodeKey+1,1).and.(nodeKey+1).lt.endPoint))Then
        nodeKey = nodeKey + 1
        pointA(1) = inputNodes(nodeKey,1)
        pointA(2) = inputNodes(nodeKey,2)
        pointA(3) = inputNodes(nodeKey,3)
        pointA(4) = inputNodes(nodeKey,4)
        pointB(1) = inputNodes(nodeKey+1,1)
        pointB(2) = inputNodes(nodeKey+1,2)
        pointB(3) = inputNodes(nodeKey+1,3)
        pointB(4) = inputNodes(nodeKey+1,4)
        coefficients = SplineAB(pointA, pointB)
      End If
      yArray = CalcPolynomial (coefficients, x, 2)
      dataPoints(i,1) = x
      dataPoints(i,2) = yArray(1)
      dataPoints(i,3) = yArray(2)
      dataPoints(i,4) = yArray(3)
! Increment x
      x = x + xIncrement
    End Do
  End Function Spline

  Function VaryNode(nodeValue, varyAmount) RESULT (outputValue)  
! Used by VaryNode
    Implicit None   ! Force declaration of all variables
! Private variables   
    Real(kind=DoubleReal) :: randDouble
    Real(kind=DoubleReal) :: nodeValue, varyAmount, outputValue
! Get rand number
    Call RANDOM_NUMBER(randDouble)
    outputValue = nodeValue + varyAmount*(randDouble-0.5D0)
  End Function VaryNode
  
  Function FillSplineResponse(dataPointsIn, startIn, endIn) RESULT (dataPointsOut)   
! Fill in gaps where there was no response to adjusting the parameter
    Implicit None   ! Force declaration of all variables
! Private variables   
    Real(kind=DoubleReal), Dimension(:,:) :: dataPointsIn
    Real(kind=DoubleReal), Dimension(1:size(dataPointsIn,1),1:size(dataPointsIn,2)) :: dataPointsOut
    Integer(kind=StandardInteger) :: i, j, k, startI, endI
    Real(kind=DoubleReal) :: x, y, xA, xB, yA, yB, grad
    Integer(kind=StandardInteger), optional :: startIn, endIn
! Init
    startI = 1
    endI = size(dataPointsIn,1)
    dataPointsOut = dataPointsIn    
! optional    
    If(Present(startIn))Then
      startI = startIn
    End If
    If(Present(endIn))Then
      endI = endIn
    End If
! Loop 
    Do i=startI+1,endI-1
      If(dataPointsOut(i,2).eq.0.0D0)Then
        xA = dataPointsOut(i-1,1)
        yA = dataPointsOut(i-1,2)
! find next point that isn't 0
        k = 0
        Do j = i+1,endI-1
          k = k + 1
          If(dataPointsOut(j,2).ne.0.0D0)Then
            xB = dataPointsOut(j,1)  
            yB = dataPointsOut(j,2)      
            exit            
          End If          
        End Do
        grad = (yB-yA)/(xB-xA)
        Do j=i,i+k
          x = dataPointsOut(j,1)
          y = yA+(x-xA)*grad
          dataPointsOut(j,2) = y
        End Do
      End If   
    End Do
  End Function FillSplineResponse
  


! ------------------------------------------------------------------------!
! Matrix Functions
! ------------------------------------------------------------------------!

! Standard/Double Precision

  Function InvertMatrix(xMatrix) RESULT (xMatrixInverse)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger) :: row,col,rowb
    Integer(kind=StandardInteger) :: matrixSize
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:2*size(xMatrix,1)) :: xMatrixWorking
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,1)) :: xMatrixInverse
    Real(kind=DoubleReal), Dimension(1:2*size(xMatrix,1)) :: xMatrixRow
! matrix(row,column)
! Initialise variables
    row = 0
    rowb = 0
    col = 0
    matrixSize = size(xMatrix,1)
    xMatrixWorking = 0.0D0
    xMatrixInverse = 0.0D0
    xMatrixRow = 0.0D0
! if a square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
! Fill working array
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
        End Do
      End Do
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.eq.col)Then
            xMatrixWorking(row,col+matrixSize) = 1.0D0
          End If
        End Do
      End Do
! make lower triangle of zeros
      Do row=1,matrixSize-1
        Do rowb=row+1,matrixSize
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the lower triangle
        Do rowb=row+1,matrixSize
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! re-force zeros in the lower triangle
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.gt.col)Then
            xMatrixWorking(row,col) = 0.0D0
          End If
        End Do
      End Do
! make upper triangle of zeros
      Do row=matrixSize,2,-1
        Do rowb=row-1,1,-1
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the upper triangle
        Do rowb=row-1,1,-1
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! Divide rhs by diagonal on lhs and store in inverse
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixInverse(row,col) = 1.0D0*&
          xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
        End Do
      End Do
    End If
  End Function InvertMatrix
  
  Function TransposeMatrix(xMatrix) RESULT (xMatrixTranspose)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,2),1:size(xMatrix,1)) :: xMatrixTranspose
    Integer(kind=StandardInteger) :: row,col
! Transpose
    Do row=1,size(xMatrix,1)
      Do col=1,size(xMatrix,2)
        xMatrixTranspose(col,row) = xMatrix(row,col)
      End Do
    End Do    
  End Function TransposeMatrix
  
  Function IdentityMatrix(iMatrix) RESULT (oMatrix)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Real(kind=DoubleReal), Dimension(:,:) :: iMatrix
    Real(kind=DoubleReal), Dimension(1:size(iMatrix,1),1:size(iMatrix,2)) :: oMatrix
! Transpose
    Do row=1,size(iMatrix,1)
      Do col=1,size(iMatrix,2)
        If(row.eq.col)Then
          oMatrix(row,col) = 1.0D0
        Else
          oMatrix(row,col) = 0.0D0
        End If
      End Do
    End Do    
  End Function IdentityMatrix
  
    
  Function DiagMatrix(iMatrix) RESULT (oMatrix)
! Takes diagonal of a matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: row,col
    Real(kind=DoubleReal), Dimension(:,:) :: iMatrix
    Real(kind=DoubleReal), Dimension(1:size(iMatrix,1),1:size(iMatrix,2)) :: oMatrix
! Transpose
    If(size(iMatrix,1).ne.size(iMatrix,2))Then
      oMatrix = iMatrix
    Else
      Do row=1,size(iMatrix,1)
        Do col=1,size(iMatrix,2)
          If(row.eq.col)Then
            oMatrix(row,col) = 1.0D0*iMatrix(row,col)
          Else
            oMatrix(row,col) = 0.0D0
          End If
        End Do
      End Do  
    End If      
  End Function DiagMatrix
  
  
  Function ScalarMult(scalar,xMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j   
    Real(kind=DoubleReal) :: scalar
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix 
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix    
! Mult 2D
    Do i=1,size(xMatrix,1)
      Do j=1,size(xMatrix,2)
        oMatrix(i,j) = scalar*xMatrix(i,j)
      End Do  
    End Do
  End Function ScalarMult

  Function MatAdd(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix
! Initialise variables
    oMatrix = 0.0D0
! Add matrices
    If(size(xMatrix,1).eq.size(yMatrix,1).and.size(xMatrix,2).eq.size(yMatrix,2))Then
      Do i=1,size(xMatrix,1)
        Do j=1,size(xMatrix,2)
          oMatrix(i,j) = xMatrix(i,j) + yMatrix(i,j)
        End Do
      End Do
    End If
  End Function MatAdd

  Function MatMult(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j,k,xRows,xCols,yRows,yCols
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(yMatrix,2)) :: oMatrix
! matrix(row,column)
! Initialise variables
    oMatrix = 0.0D0
    xRows = size(xMatrix,1)  ! Output rows
    xCols = size(xMatrix,2)
    yRows = size(yMatrix,1)
    yCols = size(yMatrix,2)  ! Output cols
! Multiply matrices
    If(xCols.eq.yRows)Then
      Do i=1,xRows
        Do j=1,yCols
          Do k=1,xCols
            oMatrix(i,j) = oMatrix(i,j) + xMatrix(i,k)*yMatrix(k,j)
          End Do
        End Do
      End Do
    End If
  End Function MatMult

! Quadruple Precision

  Function InvertMatrixQ(xMatrix) RESULT (xMatrixInverse)
! Invert square matrix
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=QuadrupleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger) :: row,col,rowb
    Integer(kind=StandardInteger) :: matrixSize
    Real(kind=QuadrupleReal), Dimension(1:size(xMatrix,1),1:2*size(xMatrix,1)) :: xMatrixWorking
    Real(kind=QuadrupleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,1)) :: xMatrixInverse
    Real(kind=QuadrupleReal), Dimension(1:2*size(xMatrix,1)) :: xMatrixRow
! Initialise variables
    row = 0
    rowb = 0
    col = 0
    matrixSize = size(xMatrix,1)
    xMatrixWorking = 0.0D0
    xMatrixInverse = 0.0D0
    xMatrixRow = 0.0D0
! if a square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
! Fill working array
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
        End Do
      End Do
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.eq.col)Then
            xMatrixWorking(row,col+matrixSize) = 1.0D0
          End If
        End Do
      End Do
! make lower triangle of zeros
      Do row=1,matrixSize-1
        Do rowb=row+1,matrixSize
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the lower triangle
        Do rowb=row+1,matrixSize
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! re-force zeros in the lower triangle
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.gt.col)Then
            xMatrixWorking(row,col) = 0.0D0
          End If
        End Do
      End Do
! make upper triangle of zeros
      Do row=matrixSize,2,-1
        Do rowb=row-1,1,-1
          If(xMatrixWorking(rowb,row).ne.0.0D0)Then !Only do if necessary
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixRow(col) = 1.0D0*&
              ((1.0D0*xMatrixWorking(row,row))/(1.0D0*xMatrixWorking(rowb,row)))*&
              xMatrixWorking(rowb,col)-1.0D0*xMatrixWorking(row,col)
            End Do
! replace row values
            Do col=1,(2*matrixSize) !loop over all columns
              xMatrixWorking(rowb,col) = 1.0D0 * xMatrixRow(col)
            End Do
          End If
        End Do
! force zeros in the upper triangle
        Do rowb=row-1,1,-1
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
! Divide rhs by diagonal on lhs and store in inverse
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixInverse(row,col) = 1.0D0*&
          xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
        End Do
      End Do
    End If
  End Function InvertMatrixQ

  Function MatAddQ(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j
    Real(kind=QuadrupleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=QuadrupleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix
! Initialise variables
    oMatrix = 0.0D0
! Add matrices
    If(size(xMatrix,1).eq.size(yMatrix,1).and.size(xMatrix,2).eq.size(yMatrix,2))Then
      Do i=1,size(xMatrix,1)
        Do j=1,size(xMatrix,2)
          oMatrix(i,j) = xMatrix(i,j) + yMatrix(i,j)
        End Do
      End Do
    End If
  End Function MatAddQ

  Function MatMultQ(xMatrix,yMatrix) RESULT (oMatrix)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j,k,xRows,xCols,yRows,yCols
    Real(kind=QuadrupleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=QuadrupleReal), Dimension(1:size(xMatrix,1),1:size(yMatrix,2)) :: oMatrix
! matrix(row,column)
! Initialise variables
    oMatrix = 0.0D0
    xRows = size(xMatrix,1)  ! Output rows
    xCols = size(xMatrix,2)
    yRows = size(yMatrix,1)
    yCols = size(yMatrix,2)  ! Output cols
! Multiply matrices
    If(xCols.eq.yRows)Then
      Do i=1,xRows
        Do j=1,yCols
          Do k=1,xCols
            oMatrix(i,j) = oMatrix(i,j) + xMatrix(i,k)*yMatrix(k,j)
          End Do
        End Do
      End Do
    End If
  End Function MatMultQ
  
 

! ------------------------------------------------------------------------!
! Unit Vector Functions
! ------------------------------------------------------------------------!

  Function ColToSquare(columnMatrix) RESULT (squareMatrix)
! Converts from PWscf 6 row unit vector to standard 3x3 unit vector
    Implicit None  ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(1:6) :: columnMatrix
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: squareMatrix
    Real(kind=DoubleReal) :: a, b, c, BC, AC, AB, cosBC, &
    cosAC, cosAB, sinBC, sinAC, sinAB
! Initialise variables
    squareMatrix = 0.0D0
! Input variables
    a = columnMatrix(1)
    b = a * columnMatrix(2)
    c = a * columnMatrix(3)
    cosBC = columnMatrix(4)
    cosAC = columnMatrix(5)
    cosAB = columnMatrix(6)
    BC = acos(cosBC)
    AC = acos(cosAC)
    AB = acos(cosAB)
    sinBC = sin(BC)
    sinAC = sin(AC)
    sinAB = sin(AB)
! Fill square matrix
    squareMatrix(1,1) = 1.0D0 * a
    squareMatrix(1,2) = 0.0D0
    squareMatrix(1,3) = 0.0D0
    squareMatrix(2,1) = 1.0D0*b*cosAB
    squareMatrix(2,2) = 1.0D0*b*sinAB
    squareMatrix(2,3) = 0.0D0
    squareMatrix(3,1) = 1.0D0*c*cosAC
    squareMatrix(3,2) = 1.0D0*c*(cosBC-cosAC*cosAB)/sinAB
    squareMatrix(3,3) = 1.0D0*((1.0D0+2.0D0*cosBC*cosAC*cosAB-&
    cosBC**2-cosAC**2-cosAB**2)**0.5)/sinAB
  End Function ColToSquare

  Function SquareToCol(squareMatrix) RESULT (columnMatrix)
! Convert from square matrix to PWscf unit vector
    Implicit None  ! Force declaration of all variables
! declare variables
    Real(kind=DoubleReal), Dimension(1:6) :: columnMatrix
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: squareMatrix
    Real(kind=DoubleReal) :: a, b, c, cosBC, cosAC, cosAB
! Convert
    a = sqrt(squareMatrix(1,1)**2+squareMatrix(1,2)**2+squareMatrix(1,3)**2)
    b = sqrt(squareMatrix(2,1)**2+squareMatrix(2,2)**2+squareMatrix(2,3)**2)/a
    c = sqrt(squareMatrix(3,1)**2+squareMatrix(3,2)**2+squareMatrix(3,3)**2)/a
    cosBC = squareMatrix(2,1)*squareMatrix(3,1)+squareMatrix(2,2)*squareMatrix(3,2)+&
    squareMatrix(2,3)*squareMatrix(3,3)
    cosAC = squareMatrix(1,1)*squareMatrix(3,1)+squareMatrix(1,2)*squareMatrix(3,2)+&
    squareMatrix(1,3)*squareMatrix(3,3)
    cosAB = squareMatrix(1,1)*squareMatrix(2,1)+squareMatrix(1,2)*squareMatrix(2,2)+&
    squareMatrix(1,3)*squareMatrix(2,3)
! Store to PWscf type column vector
    columnMatrix(1) = a
    columnMatrix(2) = b
    columnMatrix(3) = c
    columnMatrix(4) = cosBC
    columnMatrix(5) = cosAC
    columnMatrix(6) = cosAB
  End Function SquareToCol

! ------------------------------------------------------------------------!
! Vector Functions
! ------------------------------------------------------------------------!

  Function CrossProduct(VectorA, VectorB) RESULT (VectorC)
! Calculates cross product of two vectors
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
! Calculate cross product
    VectorC(1) = VectorA(2)*VectorB(3)-VectorA(3)*VectorB(2)
    VectorC(2) = VectorA(3)*VectorB(1)-VectorA(1)*VectorB(3)
    VectorC(3) = VectorA(1)*VectorB(2)-VectorA(2)*VectorB(1)
  End Function CrossProduct

  Function DotProduct(VectorA, VectorB) RESULT (DotProductResult)
! Calculates dot product of two vectors
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB
    Real(kind=DoubleReal) :: DotProductResult
! Calculate dot product
    DotProductResult = VectorA(1)*VectorB(1)+VectorA(2)*VectorB(2)+VectorA(3)*VectorB(3)
  End Function DotProduct

  Function TripleProduct(VectorA, VectorB, VectorC) RESULT (TripleProductResult)
! Calculates (scalar) triple product of three vectors (resulting in volume of 3 vectors)
    Implicit None  ! Force declaration of all variables
! declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
    Real(kind=DoubleReal) :: TripleProductResult
! Calculate cross product
    TripleProductResult = DotProduct(VectorA,CrossProduct(VectorB,VectorC))
  End Function TripleProduct

  Function TripleProductSq(VectorIn) RESULT (TripleProductResult)
! Calculates (scalar) triple product of three vectors input in 3x3 matrix (resulting in volume of 3 vectors)
    Implicit None  ! Force declaration of all variables
! declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: VectorIn
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
    Real(kind=DoubleReal) :: TripleProductResult
! Calculate cross product
    Do i=1,3
      VectorA(i) = VectorIn(1,i)
      VectorB(i) = VectorIn(2,i)
      VectorC(i) = VectorIn(3,i)
    End Do
    TripleProductResult = DotProduct(VectorA,CrossProduct(VectorB,VectorC))
  End Function TripleProductSq

! ------------------------------------------------------------------------!
! Numbers
! ------------------------------------------------------------------------!
  
  Function Odd(input) RESULT (output)
! Returns true if odd, false if even  
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: input
    Real(kind=DoubleReal) :: dpA, dpB
    Logical :: output
    output = .false.
    dpA = 1.0D0*input
    dpB = 2.0D0*ceiling(dpA/2.0D0)
    If(dpB.gt.dpA)Then
      output = .true.
    End If
  End Function Odd
    
  Function Even(input) RESULT (output)
! Returns true if even, false if odd  
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: input
    Real(kind=DoubleReal) :: dpA, dpB
    Logical :: output
    output = .true.
    dpA = 1.0D0*input
    dpB = 2.0D0*ceiling(dpA/2.0D0)
    If(dpB.gt.dpA)Then
      output = .false.
    End If
  End Function Even
  
! ------------------------------------------------------------------------!
! Random Number Related Functions
! ------------------------------------------------------------------------!

  Function RandomInteger(lower,upper) RESULT (randInt)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: lower, upper
    Integer(kind=StandardInteger) :: randInt
    Integer(kind=StandardInteger) :: diff, tempInt
    Real(kind=DoubleReal) :: randDouble
! Make random integer
    If(lower.gt.upper)Then
      tempInt = lower
      lower = upper
      upper = tempInt
    End If
    diff = (upper - lower) + 1
    Call RANDOM_NUMBER(randDouble)
    randInt = lower+floor(1.0D0*diff*randDouble)
  End Function RandomInteger

  Function RandomFloat(lower,upper) RESULT (randFloat)
! force declaration of all variables
    Implicit None
! declare variables
    Real(kind=DoubleReal) :: lower, upper
    Real(kind=DoubleReal) :: randFloat
    Real(kind=DoubleReal) :: diff, tempFloat
    Real(kind=DoubleReal) :: randDouble
! Make random integer
    If(lower.gt.upper)Then
      tempFloat = lower
      lower = upper
      upper = tempFloat
    End If
    diff = (upper - lower) + 1
    Call RANDOM_NUMBER(randDouble)
    randFloat = lower+1.0D0*diff*randDouble
  End Function RandomFloat
  
  Function RandomLCG(seedIn) RESULT (output) 
! Random number - linear congruential generator
    Implicit None ! Force declaration of all variables
! Declare variables  
    Integer(kind=LongInteger) :: m, a, c, clockTime
    Integer(kind=LongInteger) :: seed    
    Integer(kind=StandardInteger), Optional :: seedIn    
    Real(kind=DoubleReal) :: output
! init
    seed = 0
    output = 0.0D0
    m = 4294967296_LongInteger
    a = 1103515245_LongInteger
    c = 12345_LongInteger
! Read input, reset counter       
    If(Present(seedIn))Then
      seed = seedIn
      If(seed.lt.0)Then ! If seed -1 (for example) get random seed
        Call SYSTEM_CLOCK(clockTime) ! "nano seconds" - well, an estimate
        seed = mod(clockTime,m) ! fit in m
      End If
      randomLCG_n = 0
    End If  
! If first iteration
    If(randomLCG_n.eq.0)Then
      If(seed.eq.0)Then
        seed = 12791244 ! Use default seed
      End If
      randomLCG_n = 0
      randomLCG_xn = seed
    End If
! Increment counter    
    randomLCG_n = randomLCG_n + 1
! calculate
    randomLCG_xn = mod((a*randomLCG_xn+c),m)
    output = (1.0D0*randomLCG_xn)/(1.0D0*m)
  End Function RandomLCG
  
  Function RandomDist(distTypeIn,setupDistIn,lowerIn,upperIn,sigmaIn) RESULT (output)
! Random float
! F flat distribution
! S square root (diagonal) 
! G Gaussian - Box Muller - fits mu=2.5 sigma=0.5 multiplied by 0.05
! H Half Gaussian - Box Muller - fits mu=0.0 sigma=1.0 multiplied by 0.1
! P Test distribution
! M Maxwell-Boltzman distribution
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: lower, upper, sigma
    Character(len=1) :: distType, setupDist  !F flat, G Inverse gaussian distribution
    Real(kind=DoubleReal) :: output
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: randNumber, randNumberB, r, theta, x, y
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), optional :: lowerIn, upperIn, sigmaIn
    Character(len=1), optional :: distTypeIn, setupDistIn 
    Real(kind=DoubleReal), Dimension(1:6,1:2) :: fitPoints
    Real(kind=DoubleReal), Dimension(1:21,1:2) :: mbPoints
! Optional arguments  
    distType = "F"
    setupDist = "N"
    sigma = 0.4D0
    lower = 0.0D0
    upper = 1.0D0
    If(Present(distTypeIn))Then
      distType = distTypeIn
    End If
    If(Present(setupDistIn))Then
      setupDist = setupDistIn
    End If
    If(Present(sigmaIn))Then
      sigma = sigmaIn
    End If
    If(Present(lowerIn))Then
      lower = lowerIn
    End If
    If(Present(upperIn))Then
      upper = upperIn
    End If
! Get random number 0<= x <=1
    randNumber = RandomLCG()  ! maths.f90
! If Square Root Type
    If(distTypeIn.eq."S")Then    
      randNumber = Sqrt(randNumber)
    End If
! If Gaussian Type - full curve, centered on 0.5
! Box Muller method - fits mu=2.5 sigma=0.5 multiplied by 0.05
    If(distTypeIn.eq."G")Then
      If(randomDist_randNumberG.gt.-1.0D0)Then
        randNumber = randomDist_randNumberG
        randomDist_randNumberG = -1.0D0
      Else
! Get second random number
        randNumberB = RandomLCG()
! Calculate x and y        
        r = sqrt(-2.0D0*log(randNumber))
        theta = 2.0D0*pi*randNumberB
        x = r*cos(theta)
        y = r*sin(theta)
! Adjust values
        x = (x+5.0D0)/10.0D0
        y = (y+5.0D0)/10.0D0        
! store second random number
        randNumber = x
        randomDist_randNumberG = y
      End If
    End If
! If Gaussian Type - half curve, centered on 0.0
! Box Muller method - fits mu=0.0 sigma=1.0 multiplied by 0.1
    If(distTypeIn.eq."H")Then
      If(randomDist_randNumberH.gt.-1.0D0)Then
        randNumber = randomDist_randNumberH
        randomDist_randNumberH = -1.0D0
      Else
! Get second random number
        randNumberB = RandomLCG()  ! maths.f90
! Calculate x and y        
        r = sqrt(-2.0D0*log(randNumber))
        theta = 2.0D0*pi*randNumberB
        x = r*cos(theta)
        y = r*sin(theta)
! Adjust values
        x = abs(x)/5.0D0
        y = abs(y)/5.0D0        
! store second random number
        randNumber = x
        randomDist_randNumberH = y
      End If
    End If
! Testing dist type - points
    If(distTypeIn.eq."P")Then   
      If(setupDist.eq."Y")Then 
! Example set of points      
        fitPoints(1,1) = 0.0D0  
        fitPoints(2,1) = 0.2D0  
        fitPoints(3,1) = 0.4D0  
        fitPoints(4,1) = 0.6D0  
        fitPoints(5,1) = 0.8D0  
        fitPoints(6,1) = 1.0D0
        fitPoints(1,2) = 0.0D0  
        fitPoints(2,2) = 1.2D0  
        fitPoints(3,2) = 0.8D0  
        fitPoints(4,2) = 1.8D0  
        fitPoints(5,2) = 1.9D0  
        fitPoints(6,2) = 0.5D0
        randomDist_inverseInt = RandomDist_GP(fitPoints,"T")  ! maths.f90
      End If
      yArray = PointInterp(randomDist_inverseInt,randNumber,4)  ! maths.f90
      randNumber = yArray(1)
    End If  
! Maxwell-Boltzman Distribution
! P(x) = srqt(2/pi)*(x^2*exp(...))
    If(distTypeIn.eq."M")Then   
      If(setupDist.eq."Y")Then 
! a = 0.25    
        Do i=1,21
          mbPoints(i,1) = (i-0)/20.0D0
          mbPoints(i,2) = MaxwellBoltzman(mbPoints(i,1),0.25D0)
        End Do  
        randomDist_inverseInt = RandomDist_GP(mbPoints,"T")
      End If
      yArray = PointInterp(randomDist_inverseInt,randNumber,4)  ! maths.f90
      randNumber = yArray(1)
    End If  
! Output (adjust to fall in range) 
    output = lower + randNumber*(upper-lower)
  End Function RandomDist
  
  Function RandomDist_GP(inputPoints, integratorIn) RESULT (outputPoints)
! General purpose distribution function  
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal), Dimension(:,:) :: inputPoints ! 6+ pairs
    Integer(kind=LongInteger) :: i
    Real(kind=DoubleReal), Dimension(0:100,1:2) :: outputPoints
    Real(kind=DoubleReal), Dimension(1:6) :: coefficients
    Real(kind=DoubleReal), Dimension(1:7) :: coefficientsI
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: maxY, dy, y, x, xLast, iVal, iValMax
    Character(len=1), optional :: integratorIn
    Character(len=1) :: integrator
! Init    
    integrator = "P"  ! P Polynomial Fit  T Trapezoidal 
! Optional arguments        
    If(Present(integratorIn))Then
      integrator = integratorIn
    End If
! Polynomial fit + integration
    If(integrator.eq."P")Then
! Fit polynomial to the data points   
      coefficients = PolyFit(inputPoints,5)  ! maths.f90
! Integral coefficients
      coefficientsI(1) = 0.0D0
      Do i=1,6
        coefficientsI(i+1) = (1.0D0/(1.0D0*i))*coefficients(i)
      End Do
! Output points 
      maxY = CalcPolynomial(coefficientsI, 1.0D0)  ! maths.f90 ! f(x) should be always positive, Int(f(x)) always increasing
      Do i=0,100  
        outputPoints(i,1) = CalcPolynomial(coefficientsI, (i/100.0D0)) / maxY  
        outputPoints(i,2) = (i/100.0D0)
      End Do  
    End If      
! Trapezoidal integrate
    If(integrator.eq."T")Then
      dy = 1.0D0/100.0D0
      Do i=0,100
        y = (i/100.0D0)
        outputPoints(i,2) = y
        If(i.eq.0)Then
          iVal = 0.0D0
          outputPoints(i,1) = iVal
          x = inputPoints(1,1)
        Else
          yArray = PointInterp(inputPoints,y,3)  ! maths.f90
          x = yArray(1)
          If(x.lt.0.0D0)Then
            x = 0.0D0
          End If
          iVal = iVal + (0.5D0*dy*(x+xLast))
          outputPoints(i,1) = iVal
        End If
        xLast = x
      End Do
      iValMax = iVal
! Normalize so total integral = 1      
      Do i=0,100
        outputPoints(i,1) = outputPoints(i,1)/iValMax
      End Do
    End If
  End Function RandomDist_GP
  

  Function RandomVaryPoint(pointValue, maxVariation, sigma) RESULT (output)
! Vary a point +/- a max amount using inverse Gaussian distribution
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: pointValue, maxVariation, sigma
    Real(kind=DoubleReal) :: randNumber, variation, output
! Initialise variables
    output = 0.0D0
! Make varied point
    variation = RandomDist("G","N",0.0D0,maxVariation,sigma)
    Call RANDOM_NUMBER(randNumber)
    If(randNumber.gt.0.5D0)Then
      variation = -1.0D0 * variation
    End If
    output = pointValue + variation
  End Function RandomVaryPoint

  Function IntegerList(listStart,listEnd,shuffles) RESULT (list)
! Array filled with integers, possibly shuffled
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: listStart, listEnd, listSize, shuffles, rowA, rowB, shuffleCount
    Integer(kind=StandardInteger), Dimension(1:(listEnd-listStart+1)) :: list
! Initialise variables
    shuffleCount = 0
    listSize = listEnd-listStart+1
! Make list
    Do i=1,listSize
      list(i) = listStart+i-1
    End Do
! Shuffle list
    If(shuffles.gt.0)Then
      Do While(shuffleCount.lt.shuffles)
        rowA = RandomInteger(1,listSize)
        rowB = RandomInteger(1,listSize)
        If(rowA.ne.rowB)Then
          Call swapRows_Integer_1D(list,rowA,rowB)
          shuffleCount = shuffleCount + 1
        End If
      End Do
    End If
  End Function IntegerList

! ------------------------------------------------------------------------!
! Co-ordinates
! ------------------------------------------------------------------------!

  Function TransformCoords (xVect, tVect) RESULT (xPVect)
! Transform coords with transformation matrix
    Implicit None  ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(1:3) :: xVect, xPVect
    Real(kind=DoubleReal), Dimension(1:3, 1:3) :: tVect
    xPVect(1) = xVect(1)*tVect(1,1)+&
    xVect(2)*tVect(1,2)+&
    xVect(3)*tVect(1,3)
    xPVect(2) = xVect(1)*tVect(2,1)+&
    xVect(2)*tVect(2,2)+&
    xVect(3)*tVect(2,3)
    xPVect(3) = xVect(1)*tVect(3,1)+&
    xVect(2)*tVect(3,2)+&
    xVect(3)*tVect(3,3)
  End Function TransformCoords
 
  Function RdCoords (xVect, yVect) RESULT (rd)
! Transform coords with transformation matrix
    Implicit None  ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(1:3) :: xVect, yVect
    Real(kind=DoubleReal) :: rd
    rd = sqrt((xVect(1)-yVect(1))**2+(xVect(2)-yVect(2))**2+(xVect(3)-yVect(3))**2)
  End Function RdCoords
 
  Function HeatCoords (inCoords, maxVar) RESULT (outCoords)
! Transform coords with transformation matrix
    Implicit None  ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(:,:) :: inCoords
    Real(kind=DoubleReal), Dimension(1:size(inCoords,1),1:size(inCoords,2)) :: outCoords
    Real(kind=DoubleReal) :: maxVar
    Integer(kind=StandardInteger) :: i, j
! Vary Coords
    Do i=1,size(inCoords,1)
      Do j=1,size(inCoords,2)      
        outCoords(i,j) = inCoords(i,j) + 2.0D0 * (RandomDist("G") - 0.5D0) *  maxVar      
      End Do
    End Do  
    
  End Function HeatCoords
   
   
   
  

! ------------------------------------------------------------------------!
! Physical/Scientific functions
! ------------------------------------------------------------------------!

  Function Zbl (x, qA, qB) RESULT (y)
! ZBL potential, separation x, charges qA and qB
    Implicit None  ! Force declaration of all variables
! declare variables
    Integer(kind=StandardInteger) :: qA, qB
    Real(kind=DoubleReal) :: xVal, x, y, xa, xs, exa
! Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
    Else
      xVal = x
    End If
! Calculate y
    xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
    xa = 1.0D0*xVal/xs
    exa = 0.1818D0*exp(-3.2D0*xa)+0.5099D0*exp(-0.9423D0*xa)+&
    0.2802D0*exp(-0.4029*xa)+0.02817*exp(-0.2016D0*xa)
    y = ((1.0D0*qA*qB)/xVal)*exa
  End Function Zbl

  Function ZblFull (x, qA, qB) RESULT (yArray)
! y(x), y'(x), y''(x)
    Implicit None ! Force declaration of all variables
! declare variables
    Integer(kind=StandardInteger) :: qA, qB
    Real(kind=DoubleReal) :: xVal, x, y, dy, ddy, xs
    Real(kind=DoubleReal) :: termFa, termFb, termFc, termGa, termGb, termGc
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
    Else
      xVal = x
    End If
    xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
! Calculate y
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    y = termFa * termGa
    yArray(1) = y
! Calculate dy
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                 !f'(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&            !g'(x)
    (-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
    dy = termFa*termGb+termFb*termGa
    yArray(2) = dy
! Calculate ddy
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termFb = (1.0D0*qA*qB)*(xVal)**(-2.0D0)*(-1.0D0)                        !f'(x)
    termFc = (1.0D0*qA*qB)*(xVal)**(-3.0D0)*(-1.0D0)*(-2.0D0)               !f''(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                             !g(x)
    0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    0.02817*exp((-0.2016D0/xs)*xVal)
    termGb = (-3.2D0/xs)*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g'(x)
    (-0.9423D0/xs)*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)*0.02817*exp((-0.2016D0/xs)*xVal)
    termGc = (-3.2D0/xs)**2*0.1818D0*exp((-3.2D0/xs)*xVal)+&                 !g''(x)
    (-0.9423D0/xs)**2*0.5099D0*exp((-0.9423D0/xs)*xVal)+&
    (-0.4029D0/xs)**2*0.2802D0*exp((-0.4029D0/xs)*xVal)+&
    (-0.2016D0/xs)**2*0.02817*exp((-0.2016D0/xs)*xVal)
    ddy = termFa*termGc+2*termFb*termGb+termFc*termGa
    yArray(3) = ddy
  End Function ZblFull

! ------------------------------------------------------------------------!
! Decay Functions
! ------------------------------------------------------------------------!

  Function CalcIsotopeAmount(w,decayDataArray,t,calcOptionIn) RESULT (isotopeChange)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i,j,decaySteps,decayStepCounter, noChanges
    Integer(kind=StandardInteger), optional :: calcOptionIn
    Integer(kind=StandardInteger) :: calcOption
    Real(kind=DoubleReal) :: halfLifeChange, randNumber, w, t
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: decayDataArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChange
    Real(kind=DoubleReal) :: stableLimit
! Quadrupole Reals
    Real(kind=QuadrupoleReal) :: resultQ, resultGS, tQ, tempQ
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: E ! Exp
    Real(kind=QuadrupoleReal), Dimension(1:19) :: B ! Exp
! -------------------------------------------------
! decaySteps really means decay isotopes in chain (steps = decaySteps-1)
! -------------------------------------------------
! Input decay chain array:
! decayDataArray(i,1) !Tally key
! decayDataArray(i,2) !No. Atoms
! decayDataArray(i,3) !Half life
! decayDataArray(i,4) !branching factor
! decayDataArray(i,5) !isotope Z
! decayDataArray(i,6) !isotope A
! -------------------------------------------------
! Output decay chain array:
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate
! isotopeChange(i,11)   !Time
! isotopeChange(i,12)   !GS End
! -------------------------------------------------
! Optional arguments
    calcOption = 1  !(1) 1-4 analytic 5+ SG, (2)1+  SG,  (3) 1-4 analytic+GS 5+ SG
    If(Present(calcOptionIn))Then
      calcOption = calcOptionIn
    End If
! Init variables
    tQ = t
    resultQ = 0.0D0
    resultGS = 0.0D0
! -------------------------------------------------
! Alter decay chain
! -------------------------------------------------
! - If dTime * decay constant lt 1.0D-14 then assume stable for purposes of simulation
    decayStepCounter = 0
    Do i=1,size(decayDataArray,1)
      stableLimit = (log(2.0D0)/decayDataArray(i,3))*t
      decayStepCounter = decayStepCounter + 1
      If(stableLimit.lt.1.0D-14)Then
        decayDataArray(i,3) = -1    !set as stable
        Exit
      End If
    End Do
! Resize array
    decayDataArray = ArraySize2DDouble(decayDataArray,decayStepCounter)
! -------------------------------------------------
! Set stable isotope decay constant very small to avoid infinity error
! -------------------------------------------------
    Do i=1,size(decayDataArray,1)
      If(decayDataArray(i,3).eq.(-1))Then
        decayDataArray(i,3) = 1.0D100
      End If
    End Do
! -------------------------------------------------
! Break same decay constants by ~1E-3% to avoid singularities
! -------------------------------------------------
    noChanges = 0
    Do While(noChanges.eq.0)
      noChanges = 1
      Do i=1,size(decayDataArray,1)
        Do j=1,size(decayDataArray,1)
          If(i.ne.j)Then
            If(decayDataArray(i,3).eq.decayDataArray(j,3))Then
              Call RANDOM_NUMBER(randNumber)
              halfLifeChange = 0.1D0+randNumber*0.9D0
              halfLifeChange = decayDataArray(i,3)*1D-5*halfLifeChange
              decayDataArray(i,3) = decayDataArray(i,3)+halfLifeChange
              decayDataArray(j,3) = decayDataArray(j,3)-halfLifeChange
              noChanges = 0
            End If
          End If
        End Do
      End Do
    End Do
! set decay steps/isotopes
    decaySteps = size(decayDataArray,1)
! allocate isotopeChange array
    Allocate(isotopeChange(1:decaySteps,1:12))
! Fill with starting data
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
      isotopeChange(i,10) = w
      isotopeChange(i,11) = t
      isotopeChange(i,12) = 0.0D0          !default no change
    End Do
! Store lambda starting atom number data
    Do i=1,decaySteps
      If(decayDataArray(i,3).gt.9.9D99)Then
        L(i) = 0.0D0
      Else
        L(i) = lnTwoQ/isotopeChange(i,7)
      End If
      N(i) = isotopeChange(i,3)
      tempQ = -1.0D0*L(i)*tQ
      E(i) = exp(tempQ)
      B(i) = decayDataArray(i,4)
    End Do
!
! nP -> nA -> nB -> nC -> nD ...
!
! Set starting variables
    If(decaySteps.ge.1)Then
! calc nP
      If(calcOption.eq.1)Then
        resultQ = (w/L(1))*(1-E(1))+N(1)*E(1)
        isotopeChange(1,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,1,isotopeChange)
        isotopeChange(1,4) = dble(resultGS)
      End If
    End If
    If(decaySteps.ge.2)Then
! calc nA
      If(calcOption.eq.1)Then ! solve numerically
        resultQ = B(2)*L(1)*w*(1.0D0/(L(1)*L(2))+E(1)/(L(1)*(L(1)-L(2)))-&
        E(2)/(L(2)*(L(1)-L(2))))+&
        B(2)*L(1)*N(1)*(E(1)/(L(2)-L(1))+E(2)/(L(1)-L(2)))+&
        N(2)*E(2)
        isotopeChange(2,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,2,isotopeChange)
        isotopeChange(2,4) = dble(resultGS)
      End If
    End If
    If(decaySteps.ge.3)Then
! child B terms
      If(calcOption.eq.1)Then
        resultQ = &
        w*B(2)*B(3)*L(1)*L(2)*&                   ! Term 1
        (1.0D0/(L(1)*L(2)*L(3))-&
        E(1)/(L(1)*(L(1)-L(2))*(L(1)-L(3)))+&
        E(2)/(L(2)*(L(1)-L(2))*(L(2)-L(3)))+&
        E(3)/(L(3)*(L(1)-L(3))*(L(3)-L(2))))+&
        B(2)*B(3)*L(1)*L(2)*N(1)*&                ! Term 2
        (E(1)/((L(1)-L(2))*(L(1)-L(3)))-&
        E(2)/((L(1)-L(2))*(L(2)-L(3)))-&
        E(3)/((L(1)-L(3))*(L(3)-L(2))))+&
        B(3)*L(2)*N(2)*&                          ! Term 3
        (E(1)/(L(2)-L(1))+E(2)/(L(1)-L(2)))+&
        N(3)*E(3)
        isotopeChange(3,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,3,isotopeChange)
        isotopeChange(3,4) = dble(resultGS)
      End If
    End If
    If(decaySteps.ge.4)Then
! child C terms
      If(calcOption.eq.1)Then
        resultQ = &
        w*B(2)*B(3)*B(4)*L(1)*L(2)*L(3)*&          ! Term 1
        (&
        1.0D0/(L(1)*L(2)*L(3)*L(4))&
        +E(1)/(L(1)*(L(1)-L(2))*(L(1)-L(3))*(L(1)-L(4)))&
        -E(2)/(L(2)*(L(1)-L(2))*(L(1)-L(3))*(L(2)-L(4)))&
        -E(3)/(L(3)*(L(1)-L(3))*(L(3)-L(2))*(L(3)-L(4)))&
        -E(4)/(L(4)*(L(1)-L(4))*(L(4)-L(2))*(L(4)-L(3)))&
        )+&
        B(2)*B(3)*L(1)*L(2)*N(1)*&                  ! Term 2
        (&
        E(2)/((L(1)-L(2))*(L(2)-L(3))*(L(2)-L(4)))&
        -E(1)/((L(1)-L(2))*(L(1)-L(3))*(L(1)-L(4)))&
        +E(3)/((L(1)-L(3))*(L(3)-L(2))*(L(3)-L(4)))&
        +E(4)/((L(1)-L(4))*(L(4)-L(2))*(L(4)-L(3)))&
        )+&
        B(3)*B(4)*L(2)*L(3)*N(2)*&                   ! Term 3
        (&
        E(2)/((L(2)-L(3))*(L(2)-L(4)))&
        -E(3)/((L(2)-L(3))*(L(3)-L(4)))&
        -E(4)/((L(2)-L(4))*(L(4)-L(3)))&
        )+&
        B(4)*L(3)*N(3)*&                   ! Term 4
        (&
        E(3)/(L(4)-L(3))&
        +E(4)/(L(3)-L(4))&
        )+&
        E(4)*N(4)
        isotopeChange(4,4) = dble(resultQ)
      End If
      If(calcOption.eq.2.or.ISNAN(resultQ))Then ! solve numerically
        resultGS = CalcIsotopeAmountGS(tQ,4,isotopeChange)
        isotopeChange(4,4) = dble(resultGS)
      End If
    End If
! Numeric inverse laplace for remainder
    If(decaySteps.ge.5)Then
      Do i=4,decaySteps
        resultGS = CalcIsotopeAmountGS(tQ,i,isotopeChange)
        isotopeChange(i,4) = dble(resultGS)
        isotopeChange(i,12) = dble(resultGS)
      End Do
    End If
! Adjust the isotope values
    Do i=1,decaySteps
      If(isotopeChange(i,4).lt.0.0D0)Then
        isotopeChange(i,4) = 0.0D0
      End If
      If(isotopeChange(i,12).lt.0.0D0)Then
        isotopeChange(i,12) = 0.0D0
      End If
    End Do
! Store changes in isotope amounts
    Do i=1,size(isotopeChange,1)
      isotopeChange(i,2) = isotopeChange(i,4) - isotopeChange(i,3)
    End Do
  End Function CalcIsotopeAmount

  Function CalcIsotopeAmountGS(t,isotopeStep,isotopeChangeIn) RESULT (output)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: i, isotopeStep, M, k
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: isotopeChangeIn
    Real(kind=QuadrupoleReal), Dimension(1:50) :: weightingQ
    Real(kind=QuadrupoleReal), Dimension(1:20) :: L ! Lambda
    Real(kind=QuadrupoleReal), Dimension(1:20) :: N ! Starting number of atoms
    Real(kind=QuadrupoleReal), Dimension(1:20) :: B ! Starting number of atoms
    Real(kind=QuadrupoleReal) :: kQ, w, t, ft, s, FS, output
! -------------------------------------------------
! Output decay chain array:
! isotopeChange(i,1)    !Tally key
! isotopeChange(i,2)    !Change in isotope amount
! isotopeChange(i,3)    !Start amount
! isotopeChange(i,4)    !End amount
! isotopeChange(i,5)    !Isotope Z
! isotopeChange(i,6)    !Isotope A
! isotopeChange(i,7)    !T1/2
! isotopeChange(i,8)    !Decay constant
! isotopeChange(i,9)    !Branching factor
! isotopeChange(i,10)   !Parent production rate
! -------------------------------------------------
! Init variables
    M = 8
    weightingQ = GaverStehfestWeightingQ(M,weightingQ)
    w = isotopeChangeIn(1,10)
    output = 0.0D0
! Adjust the isotope values
    Do i=1,isotopeStep
      If(isotopeChangeIn(i,4).lt.0.0D0)Then
        isotopeChangeIn(i,4) = 0.0D0
      End If
    End Do
! Store lambda starting atom number data
    Do i=1,isotopeStep
      L(i) = lnTwoQ/isotopeChangeIn(i,7)
      N(i) = isotopeChangeIn(i,3)
      If(i.eq.1)Then
        B(i) = 1.0D0
      Else
        B(i) = isotopeChangeIn(i,9)
      End If
    End Do
! Perform calculation
    ft = 0.0D0
    Do k=1,2*M
      kQ = 1.0D0 * k
      s = (kQ*lnTwoQ)/t
! -----------------------
      FS = (1.0D0/(s+L(1)))*(w/s+N(1))
      Do i=2,isotopeStep
        FS = (1.0D0/(s+L(i)))*(B(i)*L(i-1)*FS+N(2))
      End Do
! FS = (1.0D0/(s+L(1)))*(w/s+N(1))
! -----------------------
      ft = ft + weightingQ(k)*FS
    End Do
    ft = (lnTwoQ/t)*ft
    output = Dble(ft)
! isotopeChangeOut(isotopeStep,4) = Dble(ft)
  End Function CalcIsotopeAmountGS

  Function GaverStehfestWeighting(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=DoubleReal) :: factor, wSum
! Real(kind=DoubleReal), Dimension(1:2*N) :: weighting
    Real(kind=DoubleReal), Dimension(:) :: weightingIn
    Real(kind=DoubleReal), Dimension(1:size(weightingIn)) :: weighting
! Init array
    weighting = 0.0D0
! k loop
    Do k=1,2*N
      factor = (-1)**(k+N)/(1.0D0*FactorialDP(N))
      jStart = Floor((k+1)/2.0D0)
      jEnd = min(k,N)
      wSum = 0.0D0
! j loop
      Do j=jStart,jEnd
        wSum = wSum + 1.0D0*(j**(N+1))*BinomialCoefficientDP(N,j)*&
        BinomialCoefficientDP(2*j,j)*BinomialCoefficientDP(j,k-j)
      End Do
      weighting(k) = factor*wSum
    End Do
  End Function GaverStehfestWeighting

  Function GaverStehfestWeightingQ(N, weightingIn) RESULT (weighting)
! Force declaration of all variables
    Implicit None
! Declare variables
    Integer(kind=StandardInteger) :: N
    Integer(kind=StandardInteger) :: j, k, jStart, jEnd
    Real(kind=QuadrupoleReal) :: factor, wSum
    Real(kind=QuadrupoleReal), Dimension(:) :: weightingIn
    Real(kind=QuadrupoleReal), Dimension(1:size(weightingIn)) :: weighting
! Init array
    weighting = 0.0D0
! k loop
    Do k=1,2*N
      factor = (-1)**(k+N)/(1.0D0*FactorialQ(N))
      jStart = Floor((k+1)/2.0D0)
      jEnd = min(k,N)
      wSum = 0.0D0
! j loop
      Do j=jStart,jEnd
        wSum = wSum + 1.0D0*(j**(N+1))*BinomialCoefficientQ(N,j)*&
        BinomialCoefficientQ(2*j,j)*BinomialCoefficientQ(j,k-j)
      End Do
      weighting(k) = factor*wSum
    End Do
  End Function GaverStehfestWeightingQ

  Function MaxTrajDepth(coefficients, maxDepthIn) RESULT (maxDepth)
! Calc max depth (E = 0) for ion trajectory, described by polynomial
    Implicit None  ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, j, nMax
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: c, x, y, maxDepth, maxDepthL
    Real(kind=DoubleReal), Optional :: maxDepthIn
! Set optional argument
    maxDepth = 1.0D10     ! 1m in ang
    If(Present(maxDepthIn))Then
      maxDepth = maxDepthIn
    End If
! Do three refinement loops
    Do i=1,4
      nMax = 20+(i*10)
      c = 10**(log10(maxDepth)/(1.0D0*nMax))
      Do j=1,50
        If(i.lt.3)Then
          x = 1.0D0*c**j
        Else
          x = (maxDepth+1.0D0)-c**(nMax-j)
        End If
        y = CalcPolynomial(coefficients, x)
        If(y.lt.0.0D0)Then
          maxDepth = x
          Exit
        End If
        If(i.eq.4)Then
          maxDepthL = x
        End If
      End Do
    End Do
    maxDepth = SolvePolynomial (coefficients, maxDepthL, maxDepth, 1.0D-6)
  End Function MaxTrajDepth

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! ------------------------------------------------------------------------!
! Swap Matrix Rows
! ------------------------------------------------------------------------!

! Integer, 1D:
  Subroutine swapRows_Integer_1D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:)
    Integer(kind=StandardInteger) :: temp,rowA,rowB
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Integer_1D

! Integer, 2D:
  Subroutine swapRows_Integer_2D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:,:)
    Integer(kind=StandardInteger) :: temp, i, rowA, rowB
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Integer_2D

! Single, 1D:
  Subroutine swapRows_Single_1D(matrix,rowA,rowB)
    Real(kind=SingleReal) :: matrix(:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=SingleReal) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Single_1D

! Single, 2D:
  Subroutine swapRows_Single_2D(matrix,rowA,rowB)
    Real(kind=SingleReal) :: matrix(:,:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=SingleReal) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Single_2D

! Double, 1D:
  Subroutine swapRows_Double_1D(matrix,rowA,rowB)
    Real(kind=DoubleReal) :: matrix(:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=DoubleReal) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Double_1D

! Double, 2D:
  Subroutine swapRows_Double_2D(matrix,rowA,rowB)
    Real(kind=DoubleReal) :: matrix(:,:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=DoubleReal) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Double_2D

! Integer, 1D:
  Subroutine sort_Integer_1D(list)
    Integer(kind=StandardInteger) :: list(:)
    Integer(kind=StandardInteger) :: i, sortComplete
! Sort list
    sortComplete = 0
    Do While(sortComplete.eq.0)
      sortComplete = 1
      Do i=2,size(list,1)
        If(list(i-1).gt.list(i))Then
          Call swapRows_Integer_1D(list,i-1,i)
          sortComplete = 0
        End If
      End Do
    End Do
  End Subroutine sort_Integer_1D

! Integer, 2D:
  Subroutine sort_Integer_2D(list, sortRow)
    Integer(kind=StandardInteger) :: list(:,:)
    Integer(kind=StandardInteger) :: i, sortRow, sortComplete
! Sort list
    sortComplete = 0
    Do While(sortComplete.eq.0)
      sortComplete = 1
      Do i=2,size(list,1)
        If(list(i-1,sortRow).gt.list(i,sortRow))Then
          Call swapRows_Integer_2D(list,i-1,i)
          sortComplete = 0
        End If
      End Do
    End Do
  End Subroutine sort_Integer_2D

! Subroutine SwapRows_Real(matrix,rowA,rowB)
!  Real(kind=SingleReal), Intent(in) :: matrix(:),rowA,rowB

! End Subroutine SwapRows_Real

! Subroutine SwapRows_Double(matrix,rowA,rowB)
!  Real(kind=DoubleReal), Intent(in) :: matrix(:),rowA,rowB

! End Subroutine SwapRows_Double

! ------------------------------------------------------------------------!
! Random Number Related Subroutines
! ------------------------------------------------------------------------!

  Subroutine SetRandomSeedArray()
    Integer(kind=StandardInteger) :: i,j,k
    Integer(kind=StandardInteger) :: clockReturn, gap, seedTemp, multiple
    Integer(kind=StandardInteger), Dimension(1:1000) :: randomSeed
    Integer(kind=StandardInteger), Dimension(0:9) :: randomIntegers
! random number seed from cpu
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
    Do i=1,1000
      k = j + gap
      If(k.gt.9)Then
        k = k - 10
      End If
      seedTemp = i * (randomIntegers(j)-randomIntegers(k))
      seedTemp = abs(seedTemp)
      multiple = floor(1.0D0 * seedTemp / 1.0D0 * clockReturn)
      seedTemp = seedTemp - multiple * clockReturn
      randomSeed(i) = abs(clockReturn - seedTemp)
      If(j.eq.9)Then
        j = 0
      End If
      j = j + 1
    End Do
    Call RANDOM_SEED(put=randomSeed)
  End Subroutine SetRandomSeedArray
  
  Subroutine CompleteNodeData(splineNodes, startIn, endIn)
! Complete y'(x) and y''(x) values for a y(x) spline
! At least 4 data points required, otherwise exits
! At least 4 "columns" in splineNodes x, y(x), y'(x). y''(x)
    Implicit None  ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal), Dimension(:,:) :: splineNodes
    Integer(kind=StandardInteger) :: node, nodes, i, j
    Integer(kind=StandardInteger), optional :: startIn, endIn
    Integer(kind=StandardInteger) :: startNode, endNode, tempNode
    Integer(kind=StandardInteger) :: interpStart, interpEnd
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: interpNodes
! optional arguments
    startNode = 1
    endNode = size(splineNodes,1)
    If(present(startIn))Then
      If(startIn.ge.1)Then
        startNode = startIn
      End If
    End If
    If(present(endIn))Then
      If(endIn.le.endNode)Then
        endNode = endIn
      End If
    End If
! Swap start/end if wrong way around
    If(startNode.gt.endNode)Then
      tempNode = startNode
      startNode = endNode
      endNode = tempNode
    End If
! Exit subroutine if too few points
    nodes = endNode-startNode+1
    If(nodes.lt.4)Then
      return
    End If    
! Interp at each node  
    Do node=startNode,endNode
! Set start-end nodes used for interpolation    
      interpStart = node-2
      interpEnd = interpStart+3
! Adjust these nodes so they are within the bounds of the data
      If(interpStart.lt.startNode)Then
        interpStart = startNode
        interpEnd = interpStart + 3
      End If
      If(interpEnd.gt.endNode)Then
        interpEnd = endNode
        interpStart = interpEnd - 3
      End If
! make interpolation array      
      j = 0
      Do i=interpStart,interpEnd
        j = j + 1
        interpNodes(j,1) = splineNodes(i,1)
        interpNodes(j,2) = splineNodes(i,2)
      End Do
      !print *,node,splineNodes(node,1),splineNodes(node,2),splineNodes(node,3),splineNodes(node,4)
! lagrange interpolation
      splineNodes(node,2) = InterpLagrange(splineNodes(node,1), interpNodes, 0)
      splineNodes(node,3) = InterpLagrange(splineNodes(node,1), interpNodes, 1)
      splineNodes(node,4) = InterpLagrange(splineNodes(node,1), interpNodes, 2)
      !print *,node,splineNodes(node,1),splineNodes(node,2),splineNodes(node,3),splineNodes(node,4)
      !print *,""
    End Do
  End Subroutine CompleteNodeData

! ------------------------------------------------------------------------!
! Miscellaneous Functions (Activity calc)
! ------------------------------------------------------------------------!

  Function ArraySize1DDouble (inputArray,arraySize) RESULT (outputArray)
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: arraySize
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:arraySize))
! transfer data
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
! force declaration of all variables
    Implicit None
! declare variables
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: arraySizeHeight
    Integer(kind=StandardInteger), optional :: arraySizeWidthIn
    Integer(kind=StandardInteger) :: arraySizeWidth
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: outputArray
! catch optional width
    If(present(arraySizeWidthIn))Then
      arraySizeWidth = arraySizeWidthIn
    Else
      arraySizeWidth = size(inputArray,2)
    End If
! Allocate output array
    Allocate(outputArray(1:arraySizeHeight,1:arraySizeWidth))
! transfer data
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
  
  

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE [TEST] FUNCTIONS                                                 !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!
  
  Function expFit_NG(dataPoints) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, iA, m, n, maxLoops
    Real(kind=DoubleReal) :: rss, lastRSS, optRSS, convergence, maxRSSVal
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:2) :: parameters, parametersOpt
    Real(kind=DoubleReal), Dimension(1:3) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: gridA
    Real(kind=DoubleReal) :: a_T, lA_T
    Real(kind=DoubleReal), Dimension(1:100,1:3) :: topBoxes 
    Integer(kind=StandardInteger) :: topBoxCount, maxRSS 
    Logical :: storeFlag
! Newton Gauss
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:2) :: J
    Real(kind=DoubleReal), Dimension(1:2,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:2) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:2) :: P      ! Change
!--------------------------------------------------
! Find Starting Parameters  
!--------------------------------------------------
! Init
    topBoxCount = 4
    topBoxes = 2.0D20
! Set a ranges
    gridA = 12
    Do i=1,gridA
      aRange(i,1) = -1.0D0*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do n=1,1000
        Do m=1,3
          a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
          lA_T = 40.0D0*(RandomLCG()-0.5D0)
! Calc rss
          rss = SingleDecayFitRSS(dataPoints, a_T, lA_T)  
! Check if better than boxed - update if better than already stored for this box         
          storeFlag = .true.
          Do i=1,topBoxCount
            If(topBoxes(i,2).eq.a_T.and.rss.lt.topBoxes(i,1))Then
              topBoxes(i,1) = rss
              topBoxes(i,2) = a_T
              topBoxes(i,3) = lA_T
              storeFlag = .false.
              Exit
            End If  
            If(i.eq.1)Then
              maxRSS = 1
              maxRSSVal = topBoxes(1,1)
            Else
              If(topBoxes(i,1).gt.maxRSSVal)Then
                maxRSS = i
                maxRSSVal = topBoxes(i,1)
              End If            
            End If
          End Do
! If better than any in box
          If(storeFlag)Then
            If(rss.lt.maxRSSVal)Then
              topBoxes(maxRSS,1) = rss
              topBoxes(maxRSS,2) = a_T
              topBoxes(maxRSS,3) = lA_T
            End If
          End If          
        End Do
      End Do
    End Do
!--------------------------------------------------
! Newton Gauss Elimination
!--------------------------------------------------
    Do n=1,topBoxCount
      parameters(1) = topBoxes(n,2)
      parameters(2) = topBoxes(n,3) 
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      maxLoops = 0
      Do While(maxLoops.le.100.and.convergence.gt.1.0D-7)
        maxLoops = maxLoops + 1
        rss = 0.0D0
! Make Jacobian and Residuals matrix
        Do i=1,size(dataPoints,1)
          R(i) = parameters(1)*exp(parameters(2)*dataPoints(i,1))-dataPoints(i,2)   ! f(x)-y
          J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          rss = rss + R(i)**2
        End Do
! calculate change matrix
        !***********     
        ! P = (JTJ)^(-1)(-1*JTR)   
        !***********      
! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)  
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do
        convergence = abs(lastRSS-rss)
        lastRSS = rss
      End Do
      If(n.eq.1)Then
        optRSS = rss
        Do i=1,size(P)
          parametersOpt(i) = parameters(i)
        End Do
      Else
        If(rss.lt.optRSS)Then
          optRSS = rss
          Do i=1,size(P)
            parametersOpt(i) = parameters(i)
          End Do
        End If
      End If  
    End Do  
    output(1) = parametersOpt(1)
    output(2) = parametersOpt(2)
    output(3) = optRSS
  End Function expFit_NG     
  
  Function expFit_LMA(dataPoints) RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, iA, m, n, maxLoops
    Real(kind=DoubleReal) :: rss, lastRSS, optRSS, convergence, maxRSSVal
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:2) :: parameters, parametersOpt
    Real(kind=DoubleReal), Dimension(1:3) :: output
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: aRange
    Integer(kind=StandardInteger) :: gridA
    Real(kind=DoubleReal) :: a_T, lA_T
    Real(kind=DoubleReal), Dimension(1:100,1:3) :: topBoxes 
    Integer(kind=StandardInteger) :: topBoxCount, maxRSS 
    Logical :: storeFlag
! Newton Gauss
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1)) :: R
    Real(kind=DoubleReal), Dimension(1:size(dataPoints,1),1:2) :: J
    Real(kind=DoubleReal), Dimension(1:2,1:size(dataPoints,1)) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:2) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:2) :: P      ! Change
    Real(kind=DoubleReal), Dimension(1:2,1:2) :: D      ! Change
!--------------------------------------------------
! Find Starting Parameters  
!--------------------------------------------------
! Init
    topBoxCount = 4
    topBoxes = 2.0D20
! Set a ranges
    gridA = 12
    Do i=1,gridA
      aRange(i,1) = -1.0D0*10D0**((gridA-4)-i)
      aRange(i,2) = -1.0D0*10D0**((gridA-5)-i)
    End Do    
    Do i=1,gridA
      aRange(i+gridA,1) = 1.0D0*10D0**(i-6)
      aRange(i+gridA,2) = 1.0D0*10D0**(i-5)
    End Do
    aRange(gridA,2) = 0.0D0
    aRange(gridA+1,1) = 0.0D0 
! Reduce search region into 10 smaller "boxes"
    Do iA=1,(gridA+gridA) 
      Do n=1,1000
        Do m=1,3
          a_T = (0.25D0*m)*(aRange(iA,1)+aRange(iA,2))
          lA_T = 40.0D0*(RandomLCG()-0.5D0)
! Calc rss
          rss = SingleDecayFitRSS(dataPoints, a_T, lA_T)  
! Check if better than boxed - update if better than already stored for this box         
          storeFlag = .true.
          Do i=1,topBoxCount
            If(topBoxes(i,2).eq.a_T.and.rss.lt.topBoxes(i,1))Then
              topBoxes(i,1) = rss
              topBoxes(i,2) = a_T
              topBoxes(i,3) = lA_T
              storeFlag = .false.
              Exit
            End If  
            If(i.eq.1)Then
              maxRSS = 1
              maxRSSVal = topBoxes(1,1)
            Else
              If(topBoxes(i,1).gt.maxRSSVal)Then
                maxRSS = i
                maxRSSVal = topBoxes(i,1)
              End If            
            End If
          End Do
! If better than any in box
          If(storeFlag)Then
            If(rss.lt.maxRSSVal)Then
              topBoxes(maxRSS,1) = rss
              topBoxes(maxRSS,2) = a_T
              topBoxes(maxRSS,3) = lA_T
            End If
          End If          
        End Do
      End Do
    End Do
!--------------------------------------------------
! Newton Gauss Elimination
!--------------------------------------------------
    Do n=1,topBoxCount
      parameters(1) = topBoxes(n,2)
      parameters(2) = topBoxes(n,3) 
! NG Opt
      lastRSS = 0
      convergence = 1.0D0
      maxLoops = 0
      Do While(maxLoops.le.100.and.convergence.gt.1.0D-7)
        maxLoops = maxLoops + 1
        rss = 0.0D0
! Make Jacobian and Residuals matrix
        Do i=1,size(dataPoints,1)
          R(i) = parameters(1)*exp(parameters(2)*dataPoints(i,1))-dataPoints(i,2)   ! f(x)-y
          J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
          J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
          rss = rss + R(i)**2
        End Do
! calculate change matrix
        !***********     
        ! P = (JTJ)^(-1)(-1*JTR)   
        !***********      
! Transpose Jacobian
        JT = TransposeMatrix(J)
        JTJ = matmul(JT,J)
        D = DiagMatrix(JTJ)
        JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)
        JTR = matmul(JT,R)
        JTR = -1.0D0*JTR ! Recycle JTR var
        P = matmul(JTJ,JTR)  
        Do i=1,size(P)
          parameters(i) = parameters(i) + P(i)
        End Do
        convergence = abs(lastRSS-rss)
        lastRSS = rss
      End Do
      If(n.eq.1)Then
        optRSS = rss
        Do i=1,size(P)
          parametersOpt(i) = parameters(i)
        End Do
      Else
        If(rss.lt.optRSS)Then
          optRSS = rss
          Do i=1,size(P)
            parametersOpt(i) = parameters(i)
          End Do
        End If
      End If  
    End Do  
    output(1) = parametersOpt(1)
    output(2) = parametersOpt(2)
    output(3) = optRSS
  End Function expFit_LMA   
  
    
  Function expFit_RSS(dataPoints, a, lA) RESULT (rss)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Real(kind=DoubleReal) :: a, lA, rss, x, y
    rss = 0.0D0
    Do i=1,size(dataPoints,1)
      x = dataPoints(i,1)
      y = a*exp(lA*x)
      rss = rss + (dataPoints(i,2)-y)**2
    End Do
  End Function expFit_RSS 
  
  
  Function test_NG() RESULT (output)
! Fits double exponential to data
! f(x) = a exp(lA)
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n
    Real(kind=DoubleReal) :: rss, lastRSS, convergence, output
    Real(kind=DoubleReal) :: x, y, lambda
    Real(kind=DoubleReal), Dimension(1:20,1:2) :: dataPoints   
    Real(kind=DoubleReal), Dimension(1:4) :: parameters, parameters_Last        
! Newton Gauss
    Real(kind=DoubleReal), Dimension(1:20) :: R, R_Last
    Real(kind=DoubleReal), Dimension(1:20,1:4) :: J, J_Last
    Real(kind=DoubleReal), Dimension(1:4,1:20) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ    ! (Jacobian Transpose * Jacobian)
    Real(kind=DoubleReal), Dimension(1:4,1:4) :: JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:4) :: JTR                ! (Jacobian Transpose * Residuals)
    Real(kind=DoubleReal), Dimension(1:4) :: P      ! Change
! Data points
    Do i=1,20
      x = (i-1)/10.D0
      y = 54.0D0*exp(-3.2D0*x)-3.0D0*exp(0.002214D0*x)
      dataPoints(i,1) = x
      dataPoints(i,2) = y
    End Do
    parameters(1) = 20.0D0  ! breaks NG
    parameters(2) = 1.7D0
    parameters(3) = 2.0D0
    parameters(4) = 0.50D0    
!--------------------------------------------------
! LMA
!--------------------------------------------------
    lastRSS = 0
    convergence = 1.0D0
    lambda = 1.0D0
    Do n=1,100
      rss = 0.0D0
! Make Jacobian and Residuals matrix
      Do i=1,20
        R(i) = (parameters(1)*exp(parameters(2)*dataPoints(i,1))+&
        parameters(3)*exp(parameters(4)*dataPoints(i,1)))-dataPoints(i,2)   ! f(x)-y
        J(i,1) = exp(parameters(2)*dataPoints(i,1))  ! d/dx1
        J(i,2) = dataPoints(i,1)*parameters(1)*exp(parameters(2)*dataPoints(i,1))  ! d/dx2
        J(i,3) = exp(parameters(4)*dataPoints(i,1))  ! d/dx3
        J(i,4) = dataPoints(i,1)*parameters(3)*exp(parameters(4)*dataPoints(i,1))  ! d/dx4
        rss = rss + R(i)**2
      End Do
! Choose whether to accept update or increase/decrease lambda      
      If(n.gt.1)Then
        print *,n,rss,lastRSS
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda     
        If(rss.gt.lastRSS)Then  ! If worse...reject, and increase lambda
! Discard changes, increase lambda
          J = J_Last
          R = R_Last
          parameters = parameters_Last
          rss = lastRSS
          lastRSS = -1.0D0
          lambda = lambda * 1.5D0
        End If
        If(rss.lt.lastRSS)Then  ! If better...accept, and decrease lambda
          lambda = lambda * 0.2D0
        End If
      End If
! calculate change matrix
      !***********     
      ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
      !***********      
! Transpose Jacobian
      JT = TransposeMatrix(J)
      JTJ = matmul(JT,J)
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
      JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
      JTR = matmul(JT,R)
      JTR = -1.0D0*JTR ! Recycle JTR var
      P = matmul(JTJ,JTR)  
! convergence of RSS
      convergence = abs(lastRSS-rss)  
      print *,"Loop:         ",n
      print *,"Lambda:       ",lambda
      print *,"RSS:          ",rss
      print *,"Convergence:  ",convergence
      print *,"Parameters:   ",parameters(1),parameters(2),parameters(3),parameters(4)
      print *,"P:            ",P(1),P(2),P(3),P(4)
      print *,""
! Store last loop values
      parameters_Last = parameters
      lastRSS = rss 
      J_Last = J
      R_Last = R
! Update parameters      
      Do i=1,size(P)
        parameters(i) = parameters(i) + P(i)
      End Do   
! Breakout if convergence threshold met      
      If(convergence.lt.1.0D-7)Then
        Exit
      End If
    End Do
    output = rss    
  End Function test_NG    
  

End Module maths
