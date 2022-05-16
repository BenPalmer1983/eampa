!#
!#  Fills array with new x,y points, also dy/dx and higher if needed
!#

SUBROUTINE fill(x, y, n_interp, out_arr)
!############################################################
IMPLICIT NONE
!############################################################
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: n_interp
REAL(kind=DoubleReal), INTENT(INOUT) :: out_arr(:,:)
!############################################################
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, n, nn, arr_l, arr_w
REAL(kind=DoubleReal) :: xmult
REAL(kind=DoubleReal) :: x_min, x_max, xi
INTEGER(kind=StandardInteger) :: point_count
LOGICAL :: loop
!############################################################

! INTERP POINTS
arr_l = SIZE(out_arr, 1)
arr_w = SIZE(out_arr, 2)

x_min = minval(x)
x_max = maxval(x)
xmult = ((x_max - x_min) / (arr_l - 1.0D0))

n = 1
DO i = 1, arr_l
  ! X VAL
  xi = x_min + (i - 1.0D0) * xmult
  out_arr(i,1) = xi

  ! FIND n
  loop = .TRUE.
  IF((x(n) .LE. xi .AND. x(n+1) .GE. xi))THEN
    loop = .FALSE.
  END IF
  DO WHILE(loop)
    n = n + 1
    IF(n .GT. SIZE(x, 1))THEN
      n = SIZE(x, 1)
      loop = .FALSE.
    END IF
    IF((x(n) .LE. xi .AND. x(n+1) .GE. xi))THEN
      loop = .FALSE.
    END IF
  END DO
  
  nn = MAX(1, MIN(n - n_interp / 2, SIZE(x,1)-n_interp + 1))
  
  CALL interpn(xi, x(nn:nn+n_interp-1), y(nn:nn+n_interp-1), out_arr(i,2))
END DO


DO j = 3, arr_w
  DO i = 1, arr_l  
    nn = MAX(1, MIN(i - n_interp / 2, SIZE(x,1)-n_interp + 1))
    CALL interpndydx(out_arr(i,1), out_arr(nn:nn+n_interp-1,1), out_arr(nn:nn+n_interp-1,j-1), out_arr(i,j))    
  END DO
END DO


END SUBROUTINE fill