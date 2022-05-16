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
