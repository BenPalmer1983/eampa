SUBROUTINE interpn(xi, x, y, yi)
! INTERP N SIZES DATA ARRAY
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(OUT) :: yi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, n
REAL(kind=DoubleReal) :: li
!############################################################
n = SIZE(x,1)
yi = 0.0D0
IF (SIZE(x,1) .EQ. SIZE(y,1)) THEN
  DO i = 1, n
    li = 1.0D0
    DO j = 1, n
      IF(i .NE. j) THEN
        li = li * (xi - x(j)) / (x(i) - x(j))
      END IF
    END DO
    yi = yi + li * y(i)
  END DO
END IF
!############################################################
END SUBROUTINE interpn