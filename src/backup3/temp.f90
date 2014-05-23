Program temp

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
  
  Call SetRandomSeedArray()
  
  x = 100.0D0
  Do i=1,100
  y = VaryPointRand(x,2.0D1,0.05D0,2.0D0,0.05D0)
  !y = VaryPointRand(x,2.0D1,0.05D0)
  print *,x,y
  End Do
  
  
End