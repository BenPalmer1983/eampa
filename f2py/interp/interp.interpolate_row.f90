!#
!#  Interpolates array
!#  Format of array columns   x, y, f(x), f'(x), f''(x)...
!# and returns single value for y
!#

SUBROUTINE interpolate_row(xi, arr, n_interp, row)
!############################################################
IMPLICIT NONE
!############################################################
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: arr(:,:)
INTEGER(kind=StandardInteger), INTENT(IN) :: n_interp
REAL(kind=DoubleReal), INTENT(OUT) :: row(1:SIZE(arr, 2))
!############################################################
! PRIVATE
INTEGER(kind=StandardInteger) :: n, k
INTEGER(kind=StandardInteger) :: width, length
REAL(kind=DoubleReal) :: x_min, x_max, x_range
INTEGER(kind=StandardInteger) :: point_count
LOGICAL :: loop
!############################################################

row = 0.0D0
x_min = arr(1, 1)
x_max = arr(SIZE(arr,1), 1)
x_range = x_max - x_min
length = SIZE(arr,1) 
width = SIZE(arr,2) 

IF(length .LT. n_interp .OR. width .LT. 2)THEN
  RETURN
END IF

! Store x 
row(1) = xi

IF(length .EQ. n_interp)THEN  
  DO k = 2, width
    CALL interpn(xi, arr(:, 1), arr(:, k), row(k))
  END DO
  RETURN
END IF

!# Estimate Location
n = INT((xi / x_range) * length)
n = MAX(1, MIN(n, length-1))

DO WHILE(n .GE. 1 .AND. n .LE. length-1)
  IF (xi .EQ. arr(n,1)) THEN
    row(2:width) = arr(n, 2:width)
    RETURN
  ELSE IF (xi .EQ. arr(n + 1, 1)) THEN
    row(2:width) = arr(n+1, 2:width)
    RETURN
  ELSE IF ((xi .GT. arr(n, 1)) .AND. (xi .LT. arr(n+1, 1))) THEN
    EXIT
  ELSE IF (xi .LT. arr(n, 1)) THEN
    n = n - 1
  ELSE IF (xi .GT. arr(n+1, 1)) THEN
    n = n + 1
  END IF
END DO

! Compute row values
n = MAX(1, MIN(n - n_interp / 2, length-n_interp + 1))
DO k = 2, width
  CALL interpn(xi, arr(n:n+n_interp-1, 1), arr(n:n+n_interp-1, k), row(k))
END DO

END SUBROUTINE interpolate_row 