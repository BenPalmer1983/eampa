SUBROUTINE interpndydx(xi, x, y, ypi)
! Interpolate and return derivative at xi
IMPLICIT NONE
! IN/OUT
REAL(kind=DoubleReal), INTENT(IN) :: xi
REAL(kind=DoubleReal), INTENT(IN) :: x(:)
REAL(kind=DoubleReal), INTENT(IN) :: y(:)
REAL(kind=DoubleReal), INTENT(OUT) :: ypi
! PRIVATE
INTEGER(kind=StandardInteger) :: i, j, k, n
REAL(kind=DoubleReal) :: fx, gx, psum
!############################################################
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
!############################################################
END SUBROUTINE interpndydx