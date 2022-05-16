! Just a zero function
! f(x) = 0
! Can be used to turn off certain potential functions

SUBROUTINE zero(r, p, p_fixed, y)
!############################################################
! f(x) = 0
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
y = 0.0D0
END SUBROUTINE zero



SUBROUTINE zero_v(r, p, p_fixed, y)
!############################################################
! f(x) = 0
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL zero(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE zero_v