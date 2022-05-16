!#
!#  Interpolates y array at x and returns single value for y
!#

SUBROUTINE interpolate(xi, x, y, n_interp, d_interp, yi)
!############################################################
IMPLICIT NONE
!############################################################
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: n_interp
INTEGER(kind=StandardInteger), INTENT(IN) :: d_interp
REAL(kind=DoubleReal), INTENT(OUT) :: yi
!############################################################
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, n
REAL(kind=DoubleReal) :: li
REAL(kind=DoubleReal) :: x_min, x_max, x_range
INTEGER(kind=StandardInteger) :: point_count
LOGICAL :: loop
!############################################################

yi = 0.0D0
x_min = x(1)
x_max = x(SIZE(x,1))
x_range = x_max - x_min

IF(SIZE(x,1) .LT. n_interp)THEN
  RETURN
END IF

IF(SIZE(x,1) .EQ. n_interp)THEN
  IF(d_interp .EQ. 0)THEN
    CALL interpn(xi, x(:), y(:), yi)
  ELSE
    CALL interpndydxn(xi, x(:), y(:), d_interp, yi)
  END IF
  RETURN
END IF

!# Estimate Location
n = INT((xi / x_range) * SIZE(x,1))
n = MAX(1, MIN(n, SIZE(x,1)-1))

DO WHILE(n .GE. 1 .AND. n .LE. SIZE(x,1)-1)
  IF (xi .EQ. x(n) .AND. d_interp .EQ. 0) THEN
    yi = y(n)
    RETURN
  ELSE IF (xi .EQ. x(n + 1) .AND. d_interp .EQ. 0) THEN
    yi = y(n + 1)
    RETURN
  ELSE IF ((xi .GT. x(n)) .AND. (xi .LT. x(n+1))) THEN
    EXIT
  ELSE IF (xi .LT. x(n)) THEN
    n = n - 1
  ELSE IF (xi .GT. x(n+1)) THEN
    n = n + 1
  END IF
END DO

n = MAX(1, MIN(n - n_interp / 2, SIZE(x,1)-n_interp + 1))

IF(d_interp .EQ. 0)THEN
  CALL interpn(xi, x(n:n+n_interp-1), y(n:n+n_interp-1), yi)
ELSE
  CALL interpndydxn(xi, x(n:n+n_interp-1), y(n:n+n_interp-1), d_interp, yi)
END IF

END SUBROUTINE interpolate 