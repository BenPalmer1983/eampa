
SUBROUTINE interpn(xi, x, y, yi)
! Identity for square matrix
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(OUT) :: yi
!############################################################
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


SUBROUTINE interpndydxn(xi, x, y, interp_n, y_out)
! Interpolate and return derivative at xi
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: interp_n
REAL(kind=DoubleReal), INTENT(OUT) :: y_out
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, k, n, m, dn
REAL(kind=DoubleReal) :: fx, gx, psum, li, yi, ypi
REAL(kind=DoubleReal) :: ta(1:SIZE(y,1))
REAL(kind=DoubleReal) :: tb(1:SIZE(y,1))
!############################################################
IF(interp_n .EQ. 0)THEN
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
  y_out = yi
ELSE IF(interp_n .EQ. 1)THEN
  n = SIZE(x,1)
  IF (SIZE(x,1) .EQ. SIZE(y,1)) THEN
    ypi = 0.0D0
    Do i=1,SIZE(x,1)
      fx = 1.0D0
      gx = 0.0D0
        Do j=1,SIZE(x,1)
        If(i .NE. j) Then
          fx = fx / (x(i) - x(j))
          psum = 1.0D0
          Do k=1,SIZE(x,1)
            If((i .NE. k) .AND. (j .NE. k))Then
              psum = psum * (xi - x(k))
            End If
          End Do
          gx = gx + psum
        End If
      End Do
      ypi = ypi + fx * gx * y(i)
    End Do
  END IF
  y_out = ypi
ELSE
  n = SIZE(x,1)
  IF (SIZE(x,1) .EQ. SIZE(y,1)) THEN
    ta(:) = y(:)
    DO dn=2, interp_n
      DO m =1,n 
        ypi = 0.0D0
        Do i=1,SIZE(x,1)
          fx = 1.0D0
          gx = 0.0D0
          Do j=1,SIZE(x,1)
            If(i .NE. j) Then
              fx = fx / (x(i) - x(j))
              psum = 1.0D0
              Do k=1,SIZE(x,1)
                If((i .NE. k) .AND. (j .NE. k))Then
                  psum = psum * (x(m) - x(k))
                End If
              End Do
              gx = gx + psum
            End If
          End Do
          ypi = ypi + fx * gx * ta(i)
        End Do
        tb(m) = ypi
      END DO
      ta(:) = tb(:)
    END DO
    ! INTERP
    ypi = 0.0D0
    Do i=1,SIZE(x,1)
      fx = 1.0D0
      gx = 0.0D0
        Do j=1,SIZE(x,1)
        If(i .NE. j) Then
          fx = fx / (x(i) - x(j))
          psum = 1.0D0
          Do k=1,SIZE(x,1)
            If((i .NE. k) .AND. (j .NE. k))Then
              psum = psum * (xi - x(k))
            End If
          End Do
          gx = gx + psum
        End If
      End Do
      ypi = ypi + fx * gx * tb(i)
    End Do
    y_out = ypi
  END IF
END IF
!############################################################
END SUBROUTINE interpndydxn


SUBROUTINE interp_fill_derivatives(nodes, interp_size)
!############################################################
REAL(kind=DoubleReal) :: nodes(:,:)
INTEGER(kind=StandardInteger) :: interp_size
!############################################################
INTEGER(kind=StandardInteger) :: n, j
!############################################################
DO n = 1, SIZE(nodes, 1)
  j = MAX(1,MIN(n,SIZE(nodes, 1)-(interp_size-1))) 
  CALL interpndydxn(nodes(n, 1), nodes(j:j+(interp_size-1), 1), &
     nodes(j:j+(interp_size-1), 2), 1, nodes(n, 3))  
END DO
END SUBROUTINE interp_fill_derivatives


SUBROUTINE interp_fill_derivatives_notends(nodes, interp_size)
!############################################################
REAL(kind=DoubleReal) :: nodes(:,:)
INTEGER(kind=StandardInteger) :: interp_size
!############################################################
INTEGER(kind=StandardInteger) :: n, j
REAL(kind=DoubleReal) :: dr_start, dr_end
!############################################################
dr_start = nodes(1,3)
dr_end = nodes(SIZE(nodes, 1),3)
DO n = 1, SIZE(nodes, 1)
  j = MAX(1,MIN(n,SIZE(nodes, 1)-(interp_size-1))) 
  CALL interpndydxn(nodes(n, 1), nodes(j:j+(interp_size-1), 1), &
     nodes(j:j+(interp_size-1), 2), 1, nodes(n, 3))  
END DO
nodes(1,3) = dr_start
nodes(SIZE(nodes, 1),3) = dr_end
END SUBROUTINE interp_fill_derivatives_notends


SUBROUTINE interp_fill_d(nodes_x, nodes_y, interp_size, nodes_dydx)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: nodes_x(:)
REAL(kind=DoubleReal), INTENT(IN) :: nodes_y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: interp_size
REAL(kind=DoubleReal), INTENT(INOUT) :: nodes_dydx(:)
!############################################################
INTEGER(kind=StandardInteger) :: n, j
!############################################################
nodes_dydx(:) = 0.0D0
IF(SIZE(nodes_x,1) .EQ. SIZE(nodes_dydx,1) .AND. SIZE(nodes_y,1) .EQ. SIZE(nodes_dydx,1))THEN
  DO n = 1, SIZE(nodes_x, 1)
    j = MAX(1,MIN(n,SIZE(nodes_x, 1)-(interp_size-1))) 
    CALL interpndydxn(nodes_x(n), nodes_x(j:j+(interp_size-1)), &
      nodes_y(j:j+(interp_size-1)), 1, nodes_dydx(n))  
  END DO
END IF
END SUBROUTINE interp_fill_d


SUBROUTINE fill_d(nodes_x, nodes_y, interp_size, nodes_dydx)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: nodes_x(:)
REAL(kind=DoubleReal), INTENT(IN) :: nodes_y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: interp_size
REAL(kind=DoubleReal), INTENT(INOUT) :: nodes_dydx(:)
!############################################################
INTEGER(kind=StandardInteger) :: n, j, nend
!############################################################
nodes_dydx(:) = 0.0D0
nend = SIZE(nodes_x, 1)
IF(SIZE(nodes_x,1) .EQ. SIZE(nodes_dydx,1) .AND. SIZE(nodes_y,1) .EQ. SIZE(nodes_dydx,1))THEN
  DO n = 1, nend
    IF(n .EQ. 1)THEN
      nodes_dydx(n) = (nodes_y(2) - nodes_y(1)) / (nodes_x(2) - nodes_x(1))
    ELSE IF(n .EQ. nend)THEN
      nodes_dydx(n) = (nodes_y(nend) - nodes_y(nend-1)) / (nodes_x(nend) - nodes_x(nend-1))
    ELSE
      nodes_dydx(n) = (nodes_y(n+1) - nodes_y(n-1)) / (nodes_x(n+1) - nodes_x(n-1))
    END IF
  END DO
END IF
END SUBROUTINE fill_d


SUBROUTINE interp_fill_xy(nodes_in_x, nodes_in_y, interp_size, nodes_out_x, nodes_out_y)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: nodes_in_x(:)
REAL(kind=DoubleReal), INTENT(IN) :: nodes_in_y(:)
INTEGER(kind=StandardInteger), INTENT(IN) :: interp_size
REAL(kind=DoubleReal), INTENT(INOUT) :: nodes_out_x(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: nodes_out_y(:)
!############################################################
INTEGER(kind=StandardInteger) :: s_in, s_out
INTEGER(kind=StandardInteger) :: m, n, j
REAL(kind=DoubleReal) :: xa, xb, xr
!############################################################
nodes_out_x = 0.0D0
nodes_out_y = 0.0D0
IF(SIZE(nodes_in_x,1) .NE. SIZE(nodes_in_y,1))THEN
  RETURN
END IF
IF(SIZE(nodes_out_x,1) .NE. SIZE(nodes_out_y,1))THEN
  RETURN
END IF
s_in = SIZE(nodes_in_x,1)
s_out = SIZE(nodes_out_x,1)
xa = nodes_in_x(1)
xb = nodes_in_x(s_in)
xr = xb - xa
DO n =1, s_out
  nodes_out_x(n) = xa + ((n - 1.0D0) / (s_out - 1.0D0)) * xr
END DO
m = 1
DO n=1,s_out
  DO WHILE((m .LT. s_in) .AND. &
  (.NOT.(nodes_in_x(m) .LE. nodes_out_x(n) .AND. nodes_out_x(n) .LE. nodes_in_x(m + 1))))
    m = m + 1
  END DO
  j = MAX(1,MIN(m,s_in-(interp_size-1))) 
  CALL interpn(nodes_out_x(n), nodes_in_x(j:j+(interp_size-1)), &
      nodes_in_y(j:j+(interp_size-1)), nodes_out_y(n))    
END DO
END SUBROUTINE interp_fill_xy


SUBROUTINE interp_fill(nodes_in_xy, interp_size, nodes_out_xy)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: nodes_in_xy(:,:)
INTEGER(kind=StandardInteger), INTENT(IN) :: interp_size
REAL(kind=DoubleReal), INTENT(INOUT) :: nodes_out_xy(:,:)
!############################################################
INTEGER(kind=StandardInteger) :: s_in, s_out
INTEGER(kind=StandardInteger) :: m, n, j, k, w
REAL(kind=DoubleReal) :: xa, xb, xr
!############################################################
nodes_out_xy = 0.0D0
s_in = SIZE(nodes_in_xy,1)
s_out = SIZE(nodes_out_xy,1)
w = SIZE(nodes_out_xy,2)
IF(w .LT. 2)THEN
  RETURN
END IF
! FILL X
xa = nodes_in_xy(1,1)
xb = nodes_in_xy(s_in,1)
xr = xb - xa
DO n =1, s_out
  nodes_out_xy(n, 1) = xa + ((n - 1.0D0) / (s_out - 1.0D0)) * xr
END DO
! INTERP Y
m = 1
DO n=1,s_out
  DO WHILE((m .LT. s_in) .AND. &
  (.NOT.(nodes_in_xy(m, 1) .LE. nodes_out_xy(n,1) .AND. nodes_out_xy(n,1) .LE. nodes_in_xy(m + 1, 1))))
    m = m + 1
  END DO
  j = MAX(1,MIN(m,s_in-(interp_size-1))) 
  CALL interpn(nodes_out_xy(n,1), nodes_in_xy(j:j+(interp_size-1),1), &
      nodes_in_xy(j:j+(interp_size-1),2), nodes_out_xy(n,2))    
END DO
! INTERP DnY/DXn
DO k=3, w  
  DO n=1,s_out
    j = MAX(1,MIN(n,s_out-(interp_size-1))) 
    CALL interpndydxn(nodes_out_xy(n,1), nodes_out_xy(j:j+(interp_size-1),1), &
        nodes_out_xy(j:j+(interp_size-1),k-1), 1, nodes_out_xy(n,k))    
  END DO
END DO
END SUBROUTINE interp_fill






