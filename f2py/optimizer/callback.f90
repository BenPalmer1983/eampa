! Testing for python

SUBROUTINE foo2(f, x, y)
!############################################################
IMPLICIT NONE
!############################################################
EXTERNAL f
REAL(kind=DoubleReal) :: f
REAL(kind=DoubleReal), INTENT(IN) :: x
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
DO n=-5,5
  y = y + f(x)
END DO
!############################################################
END SUBROUTINE foo2