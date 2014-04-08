Program polyTest

! University of Birmingham
! Ben Palmer
!


!Setup Modules
Use kinds				!data kinds
Use constants			!physical constants module
Use units				!unit conversion and normalisation 
Use strings		        !string functions
Use maths				!maths functions


!force declaration of all variables
  Implicit None

!declare variables
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPoints
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: dataPointsS
  Real(kind=SingleReal), Dimension( : ), Allocatable :: coefficientsS
  Integer(kind=StandardInteger) :: i,j,k,polyOrder
  
  polyOrder=3
  Allocate(dataPoints(1:21,1:2))
  Allocate(coefficients(0:polyOrder))
  Allocate(dataPointsS(1:21,1:2))
  Allocate(coefficientsS(0:polyOrder))
  
  dataPoints(1,1) = 7970.0132
  dataPoints(1,2) = -1680.3405
  dataPoints(2,1) = 7994.188
  dataPoints(2,2) = -1680.6499
  dataPoints(3,1) = 8018.4116
  dataPoints(3,2) = -1680.9008
  dataPoints(4,1) = 8042.6865
  dataPoints(4,2) = -1681.0769
  dataPoints(5,1) = 8067.0107
  dataPoints(5,2) = -1681.01
  dataPoints(6,1) = 8091.3813
  dataPoints(6,2) = -1680.9685
  dataPoints(7,1) = 8115.8008
  dataPoints(7,2) = -1680.9734
  dataPoints(8,1) = 8140.272
  dataPoints(8,2) = -1681.0387
  dataPoints(9,1) = 8164.7896
  dataPoints(9,2) = -1680.9819
  dataPoints(10,1) = 8189.3589
  dataPoints(10,2) = -1680.8303
  dataPoints(11,1) = 8213.9756
  dataPoints(11,2) = -1680.7107
  dataPoints(12,1) = 8238.6436
  dataPoints(12,2) = -1680.7112
  dataPoints(13,1) = 8263.3584
  dataPoints(13,2) = -1680.4535
  dataPoints(14,1) = 8288.1221
  dataPoints(14,2) = -1680.1946
  dataPoints(15,1) = 8312.9365
  dataPoints(15,2) = -1679.8531
  dataPoints(16,1) = 8337.8018
  dataPoints(16,2) = -1679.5502
  dataPoints(17,1) = 8362.7168
  dataPoints(17,2) = -1679.3109
  dataPoints(18,1) = 8387.6797
  dataPoints(18,2) = -1678.856
  dataPoints(19,1) = 8412.6934
  dataPoints(19,2) = -1678.2948
  dataPoints(20,1) = 8437.7529
  dataPoints(20,2) = -1678.0349
  dataPoints(21,1) = 8462.8672
  dataPoints(21,2) = -1677.4668

  coefficients = polynomialFit(dataPoints,polyOrder)
	
  print *,""
  print *,""
  do i=0,(size(coefficients)-1)
    print *,i,coefficients(i)
  enddo
  
  print *,""
  print *,""
  

End