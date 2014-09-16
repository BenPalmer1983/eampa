Program eampa

! University of Birmingham
! Ben Palmer
!

!Setup Modules
  Use kinds        ! data kinds
  Use mpif          ! physical constants module
  Use constants      ! physical constants module
  Use units        ! unit conversion and normalisation 
  Use general        ! string functions
  Use maths        ! maths functions
  Use initialise    ! initialise program
  Use loadData        ! load important data
  Use globals      ! declare all globals
  Use readinput      ! read input
  
!Include MPI header
  Include 'mpif.h'
  
!Variables
  Integer(kind=StandardInteger) :: error  
  
  !Integer(kind=StandardInteger) :: i
  Real(kind=DoubleReal), Dimension(1:4) :: polyCoefficients
  Real(kind=DoubleReal), Dimension(1:8) :: polyCoefficientsB
  Real(kind=DoubleReal) :: x, y, dy, ddy
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: xMat, yMat
  Real(kind=DoubleReal), Dimension(1:12,1:2) :: points
  Real(kind=DoubleReal), Dimension(1:4,1:2) :: pointsB
  Integer(kind=StandardInteger), Dimension(1:4,1:2) :: pointsC
  Integer(kind=StandardInteger), Dimension(1:4) :: pointsD
  Real(kind=DoubleReal), Dimension(1:3) :: yArray
  
!Init MPI
  Call MPI_Init(error)  

! Run initialisation module:
! Make and store output/temp directories
! Create output files
! Init a file cleanup list
  Call runInitialise()
  
! Run load data module:
! Loads isotope data into 4 arrays
! Any other data useful should be loaded here
  Call loadIsotopeData()

! Read user input file:
  Call readUserInput()
  
  polyCoefficients(1) = 3.2D0  
  polyCoefficients(2) = -4.1D0  
  polyCoefficients(3) = 7.5D0  
  polyCoefficients(4) = 6.3D0
  x = 1.2D0
  y = CalcPolynomial(polyCoefficients,x,1)
  print *,y
  
  xMat(1,1) = 0.1D0
  xMat(1,2) = 0.2D0
  xMat(1,3) = 1.1D0
  xMat(2,1) = 2.4D0
  xMat(2,2) = 0.2D0
  xMat(2,3) = 0.5D0
  xMat(3,1) = 2.3D0
  xMat(3,2) = 1.1D0
  xMat(3,3) = 4.1D0
  yMat = InvertMatrix(xMat)
  print *,yMat(1,1),yMat(1,2),yMat(1,3)
  print *,yMat(2,1),yMat(2,2),yMat(2,3)
  print *,yMat(3,1),yMat(3,2),yMat(3,3)
  
  
  
  points(1,1) = 0.1
  points(2,1) = 0.2
  points(3,1) = 0.4
  points(4,1) = 0.5
  points(5,1) = 0.6
  points(6,1) = 0.8
  points(7,1) = 1.1
  points(8,1) = 1.2
  points(9,1) = 1.5
  points(10,1) = 1.5
  points(11,1) = 1.6
  points(12,1) = 1.7
  points(1,2) = -2.0
  points(2,2) = 3.1
  points(3,2) = 1.8
  points(4,2) = 4.9
  points(5,2) = 5.7
  points(6,2) = 6.2
  points(7,2) = 6.4
  points(8,2) = 5.9
  points(9,2) = 5.6
  points(10,2) = 5.5
  points(11,2) = 5.4
  points(12,2) = 5.3
  
  print *,""
  
  !Do i=1,12
  !print *,i,points(i,1),points(i,2)
  !End Do
  !print *,""
  
  x = 0.55
  yArray = PointInterp(points,x,5,1,2,8)
  
  !x = 0.01
  !yArray = PointInterp(points,x,5,2,8)
  !x = 0.21
  !yArray = PointInterp(points,x,5,2,8)
  !x = 1.65
  !yArray = PointInterp(points,x,5,2,8)
  !x = 5.21
  !yArray = PointInterp(points,x,5,2,8)
  
  !y = InterpLagrange(x, points, 1)
  !print *,y
  x = 0.3
  pointsB(1,1) = 0.1
  pointsB(2,1) = 0.2
  pointsB(3,1) = 0.4
  pointsB(4,1) = 0.5
  pointsB(1,2) = -2.0
  pointsB(2,2) = 3.1
  pointsB(3,2) = 1.8
  pointsB(4,2) = 4.9
  
  y = InterpLagrange(x, pointsB)
  dy = InterpLagrange(x, pointsB,1)
  ddy = InterpLagrange(x, pointsB,2)
  print *,x,y,dy,ddy
  
  
  
  pointsC(1,1) = 1
  pointsC(2,1) = 2
  pointsC(3,1) = 4
  pointsC(4,1) = 5
  pointsC(1,2) = -2
  pointsC(2,2) = 3
  pointsC(3,2) = 1
  pointsC(4,2) = 4
  
  
  pointsD(1) = 1
  pointsD(2) = 2
  pointsD(3) = 4
  pointsD(4) = 5
  
  !print *,pointsB(2,2),pointsB(3,2)
  
  !Call swapRows_Double_2D(pointsB,2,3)
  
  !print *,pointsB(2,2),pointsB(3,2)
  

  
  points(1,1) = 0.1
  points(2,1) = 0.2
  points(3,1) = 0.4
  points(4,1) = 0.5
  points(5,1) = 0.6
  points(6,1) = 0.8
  points(7,1) = 1.1
  points(8,1) = 1.2
  points(9,1) = 1.5
  points(10,1) = 1.5
  points(11,1) = 1.6
  points(12,1) = 1.7
  points(1,2) = -2.0
  points(2,2) = 3.1
  points(3,2) = 1.8
  points(4,2) = 4.9
  points(5,2) = 5.7
  points(6,2) = 6.2
  points(7,2) = 6.4
  points(8,2) = 5.9
  points(9,2) = 5.6
  points(10,2) = 5.5
  points(11,2) = 5.4
  points(12,2) = 5.3
  
  
  polyCoefficientsB = PolyFit(points,7)
  print *,""
  print *,"rss: ",CalcResidualSquareSum(points,polyCoefficientsB)
  
  
  points(1,1) = 0.1
  points(2,1) = 0.2
  points(3,1) = 0.4
  points(4,1) = 0.5
  points(5,1) = 0.6
  points(6,1) = 0.8
  points(7,1) = 1.1
  points(8,1) = 1.2
  points(9,1) = 1.5
  points(10,1) = 1.55
  points(11,1) = 1.6
  points(12,1) = 1.7
  points(1,2) = -2.0
  points(2,2) = 3.1
  points(3,2) = 1.8
  points(4,2) = 4.9
  points(5,2) = 5.7
  points(6,2) = 6.2
  points(7,2) = 6.4
  points(8,2) = 5.9
  points(9,2) = 5.6
  points(10,2) = 5.5
  points(11,2) = 5.4
  points(12,2) = 5.3
  
  polyCoefficientsB = SP_PolyFit(points,7)
  print *,""
  print *,"rss: ",CalcResidualSquareSum(points,polyCoefficientsB)
  
  
  !pointsD = IntegerList(3,6,10)
  
  
!Finalise MPI
  Call MPI_Finalize(error)
  
  
  
End