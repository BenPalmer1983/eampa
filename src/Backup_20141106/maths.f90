Module maths

!--------------------------------------------------------------!
! Maths functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 13th August 2014
!----------------------------------------

!Force declaration of all variables
  Implicit None
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37    
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307    
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer    1 byte
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer    4 bytes -2E31 to 2E31-1
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer 8 bytes -2E63 to 2E63-1
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer

!Make private
  Private  
!--Functions--!
! General Maths Functions
  Public :: Factorial  
  Public :: BinomialCoefficient
  Public :: RoundDP
! Polynomial Related Functions
  Public :: SolvePolynomial  
  Public :: CalcPolynomial
  Public :: DerivativePolynomial
! Fitting, Regression, Interpolation  
  Public :: PolyFit             ! Fit polynomial to large set of data - Vandermonde
  Public :: SP_PolyFit          ! Fit polynomial to large set of data - Superposition
  Public :: MinPolyFit          ! Fit polynomial, find minimum (if in region of data points)
  Public :: InterpLagrange      !find y(x) or y'(x) using lagrange interp from set of x-y
  Public :: PointInterp         !y(x) from large set of data x-y, finds region and uses lagrange interp
  Public :: CalcResidualSquareSum
  Public :: MurnFit
  Public :: MurnCalc
  Public :: MurnRSS
  Public :: MurnFitBP
  Public :: BirchMurnFit
  Public :: BirchMurnCalc
  Public :: BirchMurnRSS
  Public :: BirchMurnFitBP
! Vector Functions  
  Public :: CrossProduct
  Public :: DotProduct
  Public :: TripleProduct
  Public :: TripleProductSq
! Random Number Related Functions 
  Public :: RandomInteger
  Public :: RandomFloat
  Public :: RandomDist
  Public :: RandomVaryPoint
  Public :: IntegerList
! Matrix Functions  
  Public :: InvertMatrix
! Spline Functions
  Public :: SplineAB
  Public :: SplineNodes
  Public :: Spline
! Other Physical Functions
  Public :: Zbl  
  Public :: ZblFull  
   
!--Subroutines--!  
  Public :: setRandomSeedArray
  Public :: swapRows_Double_1D
  Public :: swapRows_Double_2D
  !Public :: sortList_Int
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE INTERFACE                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!    
  
  !Interface swapRows
  !  Module Procedure swapRows_Integer_1D 
    !swapRows_Real, swapRows_Double
  !End Interface swapRows
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
!------------------------------------------------------------------------!
! General Maths Functions
!------------------------------------------------------------------------! 

  Function Factorial(input) RESULT (output)
!force declaration of all variables
  Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,input
    Integer(kind=StandardInteger) :: output
!calculate factorial
    output = 1
    Do i=1,input
      output = i * output
    End Do  
  End Function Factorial
  
  Function BinomialCoefficient(n,k) RESULT (c)
!force declaration of all variables
    Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: c,n,k
!calculate factorial
    c = Factorial(n)/(Factorial(n-k)*Factorial(k))
  End Function BinomialCoefficient
  
  Function RoundDP(dpIn) RESULT (intOut) 
! Round DP to nearest int
    Implicit None ! Force declaration of all variables
    Real(kind=DoubleReal) :: dpIn
    Integer(kind=StandardInteger) :: intOut
    intOut = Floor(dpIn+0.5D0)
  End Function RoundDP
  
      
!------------------------------------------------------------------------!
! Polynomial Related Functions
!------------------------------------------------------------------------!     
  
  Function SolvePolynomial (coefficients, lower, upper) RESULT (output)   
!Solves the polynomial p(x) = 0, in region close to p(x) = 0
    Implicit None  !Force declaration of all variables
!Declare variables
    Real(kind=DoubleReal), Dimension(:) :: coefficients
    Real(kind=DoubleReal) :: upper, lower, output
    Real(kind=DoubleReal) :: x,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget, factor, difference
    Integer(kind=StandardInteger) :: i,maxLoops  
!Set values
    convergenceTarget = 0
    convergence = 1000
    convergenceThreshold = 0.00001
    maxLoops = 0
!set start value for x
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
      do i=2,size(coefficients)
        dydx = dydx + i * x**(i-2) * coefficients(i)      
      End Do      
      convergence = abs(convergenceTarget - y)
      if(convergence.gt.convergenceThreshold)then
        if((dydx.lt.0.and.y.ge.0).or.(dydx.ge.0.and.y.lt.0))then
          x = x + difference  
        else
          x = x - difference
        endif
      endif
    End Do  
    output = x
  End Function solvePolynomial   
  
  Function CalcPolynomial(polyCoefficientsIn, x, derivativeIn) RESULT (y)
!Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
    Implicit None  !Force declaration of all variables
!Declare variables
    Integer(kind=StandardInteger) :: i, j, derivative
    Integer(kind=StandardInteger), Optional :: derivativeIn
    Real(kind=DoubleReal), Dimension( : ) :: polyCoefficientsIn
    Real(kind=DoubleReal), Dimension(1:Size(polyCoefficientsIn,1)) :: polyCoefficients
    Real(kind=DoubleReal) :: x, y
!Handle Optional Arguments
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
!Calculates p(x) by default, p'(x) for derivativeIn = 1 etc
    Implicit None  !Force declaration of all variables
!Declare variables
    Real(kind=DoubleReal), Dimension( : ) :: coefficientsIn
    Real(kind=DoubleReal), Dimension(1:(size(coefficientsIn)-1)) :: coefficientsOut
    Integer(kind=StandardInteger) :: j
    Do j=1,size(coefficientsIn,1)-1    
      coefficientsOut(j) = j * coefficientsIn(j+1)
    End Do
 End Function DerivativePolynomial  
  
  
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
  
!------------------------------------------------------------------------! 
! Fitting, Regression, Interpolation  
!------------------------------------------------------------------------!
  Function PolyFit(points,order,extendedFitIn) RESULT (coefficients)
!Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
!Declare variables
    Integer(kind=StandardInteger) :: k,col,row,exponentValue
    Integer(kind=StandardInteger), Optional :: extendedFitIn
    Integer(kind=StandardInteger) :: extendedFit
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal) :: vRSS, spRSS
    Integer(kind=StandardInteger) :: order
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: coefficients, coefficientsSP
    Real(kind=DoubleReal), Dimension(1:(order+1),1:(order+1)) :: xMatrix
    Real(kind=DoubleReal), Dimension(1:(order+1)) :: yMatrix
!Optional argument
    extendedFit = 0
    If(Present(extendedFitIn))Then
      extendedFit = extendedFitIn
    End If
!Step 1 - Standard fit with Vandermonde matrix
!Build Least Squares Fitting Vandermonde matrix
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
!invert xMatrix
    xMatrix = InvertMatrix(xMatrix)
!multiply inverse by y to get coefficients
    coefficients = matMul(xMatrix,yMatrix)
!Step 2 - extended fit
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
!Fits a polynomial of order to the points input
    Implicit None  !Force declaration of all variables
!Declare variables
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
!Declare variables
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
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, pointToVary
    !Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
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
    !Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
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
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, pointToVary
    !Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
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
    !Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
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
  
  Function InterpLagrange(x, points, derivativeIn) RESULT (output)
! Calculates y(x) or y'(x) using Lagrange interpolation
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
!Calculate y
      y = 0.0D0
      Do n=1,size(points,1)
        y=y+points(n,2)*coefficients(n)
      End Do
      output = y
    End If
!y'(x)
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
!Calculate dy
      dy = 0.0D0
      Do n=1,size(points,1)
        dy = dy + points(n,2)*coefficients(n)
      End Do
      output = dy
    End If
!y''(x)
    If(derivative.eq.2)Then  
!Could use recursive functions for higher orders, but y''(x) high enough for now
!Calculate y'(x) for each x input
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
!Calculate dy
      dyy = 0.0D0
      Do n=1,size(points,1)
        dyy = dyy + pointsTemp(n,2)*coefficients(n)
      End Do
      output = dyy
    End If
  End Function InterpLagrange   
  
  
  Function PointInterp(points,x,subsetSize,derivativeIn,inputSetStartIn,inputSetLengthIn) RESULT (yArray)
! Takes large set of data points, finds region, and interps with lagrange
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
    Else If(subsetSize.lt.inputLength)Then
!Reduce set of data points
      xLower = points(inputStart,1)
      xUpper = points(inputEnd,1)
!Find xPos
      If(x.lt.xLower)Then  !If x lower than data set, use lowest possible points
        xPos = inputStart
      Elseif(x.gt.xUpper)Then  !If x higher than data set, use highest possible points
        xPos = inputEnd
      Else
!Estimate position
        xPos = INT(Floor(((x - xLower) / (xUpper - xLower)) * 1.0D0 * inputLength) + inputStart)
        If(xPos.lt.inputStart)Then
          xPos = inputStart
        End If
        If((xPos+1).gt.inputEnd)Then
          xPos = inputEnd-1
        End If
        xLower = points(xPos,1)
        xUpper = points(xPos+1,1)
!If estimate is incorrect, search for better value
        If (x.lt.xLower) Then
          xPosStart = xPos
          Do xPos=xPosStart,inputStart,-1    !Search down
            xLower = points(xPos,1)
            xUpper = points(xPos+1,1)
            If(x.le.xUpper.and.x.ge.xLower)Then
              Exit  !xPos found
            End If
          End Do
        End If
        If (x.gt.xUpper) Then
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
!Adjust xPos to center of subset  
      xPosOffset = INT(Floor(1.0D0*subsetSize/2))
      xPosOffsetR = subsetSize - xPosOffset
      xPosLower = xPos - xPosOffset
      xPosUpper = xPos + xPosOffsetR - 1
!Adjust xPos start, so it fits in the range of subset of selected data points
      If(xPosLower.lt.inputStart)Then
        xPosLower = inputStart
        xPosUpper = inputStart + subsetSize - 1
      End If 
      If(xPosUpper.gt.inputEnd)Then
        xPosLower = inputEnd - subsetSize + 1
        xPosUpper = inputEnd
      End If    
!Transfer data points to pointsInterp
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
  
  
  
  
  
!------------------------------------------------------------------------!
! Spline Functions
!------------------------------------------------------------------------!   

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
    
 
!------------------------------------------------------------------------!
! Matrix Functions
!------------------------------------------------------------------------! 
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
! Initialise variables   
    row = 0 
    rowb = 0 
    col = 0 
    matrixSize = size(xMatrix,1)
    xMatrixWorking = 0.0D0
    xMatrixInverse = 0.0D0  
    xMatrixRow = 0.0D0  
!if a square matrix
    If(size(xMatrix,1).eq.size(xMatrix,2))Then
!Fill working array
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixWorking(row,col) = 1.0D0*xMatrix(row,col)
        End Do
      End Do
      Do row=1,matrixSize
        Do col=1,matrixSize
          If(row.eq.col)then
            xMatrixWorking(row,col+matrixSize) = 1.0D0
          endif
        End Do
      End Do  
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
            Do col=1,(2*matrixSize) !loop over all columns
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
          If(row.gt.col)then
            xMatrixWorking(row,col) = 0.0D0
          endif
        End Do
      End Do  
!make upper triangle of zeros  
      Do row=matrixSize,2,-1
        Do rowb=row-1,1,-1    
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
        Do rowb=row-1,1,-1
          xMatrixWorking(rowb,row) = 0.0D0
        End Do
      End Do
!Divide rhs by diagonal on lhs and store in inverse
      Do row=1,matrixSize
        Do col=1,matrixSize
          xMatrixInverse(row,col) = 1.0D0*&
          xMatrixWorking(row,col+matrixSize)/xMatrixWorking(row,row)
        End Do
      End Do
    End If 
  End Function InvertMatrix
  
  Function MatAdd(xMatrix,yMatrix) RESULT (oMatrix)   
! Force declaration of all variables
  Implicit None
! Declare variables  
    Integer(kind=StandardInteger) :: i,j
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix, yMatrix
    Real(kind=DoubleReal), Dimension(1:size(xMatrix,1),1:size(xMatrix,2)) :: oMatrix
! Initialise variables   
    oMatrix = 0.0D0
!Add matrices
    If(size(xMatrix,1).eq.size(yMatrix,1).and.size(xMatrix,2).eq.size(yMatrix,2))Then
      Do i=1,size(xMatrix,1)
        Do j=1,size(xMatrix,2)
          oMatrix(i,j) = xMatrix(i,j) + yMatrix(i,j)
        End Do
      End Do    
    End If  
  End Function MatAdd
  
!------------------------------------------------------------------------!
! Unit Vector Functions  
!------------------------------------------------------------------------!   
  
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
!Fill square matrix
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
!declare variables
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
  
  
  
!------------------------------------------------------------------------!
! Vector Functions
!------------------------------------------------------------------------!   
  
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
!declare variables
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
    Real(kind=DoubleReal) :: TripleProductResult
!Calculate cross product
    TripleProductResult = DotProduct(VectorA,CrossProduct(VectorB,VectorC))
  End Function TripleProduct  
  
  Function TripleProductSq(VectorIn) RESULT (TripleProductResult)
! Calculates (scalar) triple product of three vectors input in 3x3 matrix (resulting in volume of 3 vectors)
    Implicit None  ! Force declaration of all variables
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: VectorIn
    Real(kind=DoubleReal), Dimension(1:3) :: VectorA, VectorB, VectorC
    Real(kind=DoubleReal) :: TripleProductResult
!Calculate cross product
    Do i=1,3
      VectorA(i) = VectorIn(1,i)
      VectorB(i) = VectorIn(2,i)
      VectorC(i) = VectorIn(3,i)
    End Do
    TripleProductResult = DotProduct(VectorA,CrossProduct(VectorB,VectorC))
  End Function TripleProductSq 
  
  
!------------------------------------------------------------------------!  
! Random Number Related Functions
!------------------------------------------------------------------------!  
    
  Function RandomInteger(lower,upper) RESULT (randInt)
!force declaration of all variables
  Implicit None
!declare variables
  Integer(kind=StandardInteger) :: lower, upper
  Integer(kind=StandardInteger) :: randInt
  Integer(kind=StandardInteger) :: diff, tempInt
  Real(kind=DoubleReal) :: randDouble
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
  
  Function RandomFloat(lower,upper) RESULT (randFloat)
!force declaration of all variables
  Implicit None
!declare variables
  Real(kind=DoubleReal) :: lower, upper
  Real(kind=DoubleReal) :: randFloat
  Real(kind=DoubleReal) :: diff, tempFloat
  Real(kind=DoubleReal) :: randDouble
!Make random integer
    If(lower.gt.upper)Then
    tempFloat = lower
    lower = upper
    upper = tempFloat
  End If
    diff = (upper - lower) + 1
  Call RANDOM_NUMBER(randDouble)
  randFloat = lower+1.0D0*diff*randDouble  
  End Function RandomFloat
  
  Function RandomDist(lower,upper,distType,sigma) RESULT (output)
! Random float, flat distribution, distribution using inverse gaussian
    Implicit None ! Force declaration of all variables
! Declare variables  
    Real(kind=DoubleReal) :: lower, upper, sigma
    Character(len=1) :: distType   !F flat, G Inverse gaussian distribution
    Real(kind=DoubleReal) :: output
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:21,1:2) :: distributionPoints
    Real(kind=DoubleReal) :: yOne
    Real(kind=DoubleReal) :: mu
    Real(kind=DoubleReal) :: distX, distXInterval
    Real(kind=DoubleReal) :: randNumber
    Real(kind=DoubleReal), Dimension(1:3) :: yArray 
! initialise variables
    output = 0.0D0
    yArray = 0.0D0
!Flat distribution  
    If(distType.eq."F")Then
!Generate random number
      Call RANDOM_NUMBER(randNumber)
      output = lower + (upper-lower) * randNumber
    End If
!Gaussian type distribution
    If(distType.eq."G")Then
      mu = 0.0D0
!Set data points  
      distXInterval = 1.0D0/(21-1)
      yOne = exp(-1*((1.0D0-mu)**2)/(2*sigma**2))
!Make array - inverse Gaussian
      !distributionPoints(1,1) = 0.0D0        !remove when debugging finished
      !distributionPoints(1,2) = 1.0D0
      !distributionPoints(21,1) = 1.0D0
      !distributionPoints(21,2) = 0.0D0
      !Do i=2,20
      !  distX = (i-1)*1.0D0*distXInterval
      !  distributionPoints(i,1) = (exp(-1*((distX-mu)**2)/(2*sigma**2))-yOne)/(1-yOne)
      !  distributionPoints(i,2) = distX
      !End Do
      distributionPoints(1,1) = 0.0D0
      distributionPoints(1,2) = 1.0D0
      distributionPoints(21,1) = 1.0D0
      distributionPoints(21,2) = 0.0D0
      Do i=2,20
        distX = 1.0D0-(i-1)*distXInterval
        distributionPoints(i,1) = (exp(-1*((distX-mu)**2)/(2*sigma**2))-yOne)/(1-yOne)
        distributionPoints(i,2) = distX
      End Do      
!Generate random number
      Call RANDOM_NUMBER(randNumber)      
      yArray = PointInterp(distributionPoints,randNumber,3)
      output = lower + (upper-lower) * yArray(1)      
      !Do i=1,1000      !remove when debugging finished
      !  yArray = PointInterp(distributionPoints,(i/1000.0D0),3)
      !  print *,(i/1000.0D0),yArray(1)
      !End Do
    End If  
  End Function RandomDist
  
  Function RandomVaryPoint(pointValue, maxVariation, sigma) RESULT (output)
! Vary a point +/- a max amount using inverse Gaussian distribution
    Implicit None ! Force declaration of all variables
! Declare variables   
    Real(kind=DoubleReal) :: pointValue, maxVariation, sigma
    Real(kind=DoubleReal) :: randNumber, variation, output
! Initialise variables
    output = 0.0D0    
! Make varied point
    variation = RandomDist(0.0D0,maxVariation,"G",sigma)
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
  
  
!------------------------------------------------------------------------!
! Physical/Scientific functions
!------------------------------------------------------------------------! 
    
  Function Zbl (x, qA, qB) RESULT (y)
! ZBL potential, separation x, charges qA and qB
    Implicit None  ! Force declaration of all variables
!declare variables
    Integer(kind=StandardInteger) :: qA, qB
    Real(kind=DoubleReal) :: xVal, x, y, xa, xs, exa
!Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
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
  
  
  Function ZblFull (x, qA, qB) RESULT (yArray)
! y(x), y'(x), y''(x)
    Implicit None ! Force declaration of all variables
!declare variables
    Integer(kind=StandardInteger) :: qA, qB
    Real(kind=DoubleReal) :: xVal, x, y, dy, ddy, xs
    Real(kind=DoubleReal) :: termFa, termFb, termFc, termGa, termGb, termGc
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
!Force none infinite result for 0
    If(x.eq.0.0D0)Then
      xVal = 0.00001D0
    Else 
      xVal = x
    End If
    xs = 0.4683766 * (qA**(2.0D0/3.0D0)+qB**(2.0D0/3.0D0))**0.5
!Calculate y
    termFa = (1.0D0*qA*qB)*(xVal)**(-1.0D0)                          !f(x)
    termGa = 0.1818D0*exp((-3.2D0/xs)*xVal)+&                        !g(x)
             0.5099D0*exp((-0.9423D0/xs)*xVal)+&
             0.2802D0*exp((-0.4029D0/xs)*xVal)+&
             0.02817*exp((-0.2016D0/xs)*xVal)
    y = termFa * termGa
    yArray(1) = y
!Calculate dy
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
!Calculate ddy
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
  
  
  
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
!------------------------------------------------------------------------!  
! Swap Matrix Rows
!------------------------------------------------------------------------!    

! Integer, 1D:
  Subroutine swapRows_Integer_1D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:)
    Integer(kind=StandardInteger) :: temp,rowA,rowB
!Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp        
  End Subroutine swapRows_Integer_1D  
  
! Integer, 2D:  
  Subroutine swapRows_Integer_2D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:,:)
    Integer(kind=StandardInteger) :: temp, i, rowA, rowB
!Swap rows
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
!Swap rows
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
!Swap rows
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
!Swap rows
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
!Swap rows
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
!Sort list
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
!Sort list
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
  
  !Subroutine SwapRows_Real(matrix,rowA,rowB)
  !  Real(kind=SingleReal), Intent(in) :: matrix(:),rowA,rowB
        
        
  !End Subroutine SwapRows_Real   
  
  !Subroutine SwapRows_Double(matrix,rowA,rowB)
  !  Real(kind=DoubleReal), Intent(in) :: matrix(:),rowA,rowB
        
        
  !End Subroutine SwapRows_Double  
  
  
!------------------------------------------------------------------------!  
! Random Number Related Subroutines
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
  
End Module maths  


