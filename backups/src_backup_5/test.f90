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
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsA
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsB
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPointsC
  Real(kind=DoubleReal) :: valueA

  Allocate(dataPointsA(1:3,1:2))
  Allocate(dataPointsB(1:4,1:2))
  Allocate(dataPointsC(1:5,1:2))
  
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
  
  dataPointsC(1,1) = -0.5
  dataPointsC(1,2) = 2.877188
  dataPointsC(2,1) = 0
  dataPointsC(2,2) = 4.1
  dataPointsC(3,1) = 1
  dataPointsC(3,2) = 7.635
  dataPointsC(4,1) = 1.2
  dataPointsC(4,2) = 8.377376
  dataPointsC(5,1) = 1.4
  dataPointsC(5,2) = 9.100856
  
  
  !dataPoints(1,1) = -0.5
  !dataPoints(1,2) = 2.445625
  !dataPoints(2,1) = 0
  !dataPoints(2,2) = 3
  !dataPoints(3,1) = 1
  !dataPoints(3,2) = 3.535
  !dataPoints(4,1) = 1.2
  !dataPoints(4,2) = 3.56448
  
  coefficients = PointInterpolationFull(dataPointsC,1.1D0)
  print *,coefficients(1)
  print *,coefficients(2)
  !valueA = PointInterpolation(dataPointsC,1.1D0)
  !print *,valueA
  
  !print *,BinomialCoefficient(1,3)
  !print *,BinomialCoefficient(1,4)
  !print *,BinomialCoefficient(2,4)
  !print *,BinomialCoefficient(3,4)
  !print *,BinomialCoefficient(4,4)
  
  !valueA = PermutationCoefficients(dataPointsC,1,0)
  
  
End