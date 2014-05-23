Program test

! University of Birmingham
! Ben Palmer
!

!Setup Modules
  Use kinds				!data kinds
  Use constants			!physical constants module
  Use strings		        !string functions
  Use maths				!maths functions
  
!Variables
  Integer(kind=StandardInteger) :: i
  Real(kind=DoubleReal) :: x,y
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splineXY
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsA
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsB
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsC
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsD
  Real(kind=DoubleReal) :: valueA
  Real(kind=SingleReal) :: startTime, endTime
  
  Call cpu_time(startTime)

  Allocate(dataPointsA(1:3,1:2))
  Allocate(dataPointsB(1:4,1:2))
  Allocate(dataPointsC(1:5,1:2))
  Allocate(dataPointsD(1:8,1:2))
  !Allocate(coefficients(0:4))
  
  dataPointsA(1,1) = -0.5
  dataPointsA(1,2) = 10.69
  dataPointsA(2,1) = 0
  dataPointsA(2,2) = 12.4
  dataPointsA(3,1) = 1
  dataPointsA(3,2) = 13.96
  
  dataPointsB(1,1) = -0.5
  dataPointsB(1,2) = 2.445625
  dataPointsB(2,1) = 0
  dataPointsB(2,2) = 3
  dataPointsB(3,1) = 1
  dataPointsB(3,2) = 3.535
  dataPointsB(4,1) = 1.2
  dataPointsB(4,2) = 3.56448
  
  dataPointsC(1,1) = -0.5D0
  dataPointsC(1,2) = 2.877188D0
  dataPointsC(2,1) = 0.0D0
  dataPointsC(2,2) = 4.1D0
  dataPointsC(3,1) = 1.0D0
  dataPointsC(3,2) = 7.635D0
  dataPointsC(4,1) = 1.2D0
  dataPointsC(4,2) = 8.377376D0
  dataPointsC(5,1) = 1.4D0
  dataPointsC(5,2) = 9.100856D0
  
  dataPointsD(1,1) = -0.5D0
  dataPointsD(1,2) = 2.877188D0
  dataPointsD(2,1) = 0.0D0
  dataPointsD(2,2) = 4.1D0
  dataPointsD(3,1) = 1.0D0
  dataPointsD(3,2) = 7.635D0
  dataPointsD(4,1) = 1.2D0
  dataPointsD(4,2) = 8.377376D0
  dataPointsD(5,1) = 1.4D0
  dataPointsD(5,2) = 9.100856D0
  dataPointsD(6,1) = 1.5D0
  dataPointsD(6,2) = 9.03D0
  dataPointsD(7,1) = 1.6D0
  dataPointsD(7,2) = 8.76D0
  dataPointsD(8,1) = 1.8D0
  dataPointsD(8,2) = 9.321D0
    
  !coefficients = FivePointInterpolationValue(dataPointsC,1.2D0) 
  !print *,coefficients(1),coefficients(2),coefficients(3)
  
  !yArray = MatrixInterpolation(dataPointsC,1.2D0) 
  !print *,yArray(1),yArray(2),yArray(3),yArray(4),yArray(5)
  
  !print *,coefficients(1),coefficients(2),coefficients(3)
  
  !splineXY = Spline(dataPointsD, 1001, 3, 4)
  
  !print *,splineXY(1,1),splineXY(1,2)
  !print *,splineXY(500,1),splineXY(500,2)
  !print *,splineXY(1000,1),splineXY(1000,2)
  !print *,splineXY(1001,1),splineXY(1001,2)
  
  
  !dataPoints(1,1) = -0.5
  !dataPoints(1,2) = 2.445625
  !dataPoints(2,1) = 0
  !dataPoints(2,2) = 3
  !dataPoints(3,1) = 1
  !dataPoints(3,2) = 3.535
  !dataPoints(4,1) = 1.2
  !dataPoints(4,2) = 3.56448
  
  !coefficients = PointInterpolationFull(dataPointsC,1.1D0)
  !print *,coefficients(1)
  !print *,coefficients(2)
  !valueA = PointInterpolation(dataPointsC,1.1D0)
  !print *,valueA
  
  !print *,BinomialCoefficient(1,3)
  !print *,BinomialCoefficient(1,4)
  !print *,BinomialCoefficient(2,4)
  !print *,BinomialCoefficient(3,4)
  !print *,BinomialCoefficient(4,4)
  
  !valueA = PermutationCoefficients(dataPointsC,1,0)
  
  !x = 100
  !Do i=1,1000
  !  y = VaryPoint(x,0.01D0,1)
!	print *,i,x,y
  !End Do
  
  x = 1.25
  Do i=1,1000000
    yArray = PointInterpolationArr(dataPointsD,x)
  End Do
  print *,yArray(1),yArray(2)
  
  
  Call cpu_time(endTime)
  print *,(endTime-startTime)
  
  
End