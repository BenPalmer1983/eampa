! ### GRADIENT ###
! Just a zero function
! f(x) = 0
! Can be used to turn off certain potential functions

SUBROUTINE zero_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = 0
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
!############################################################
dydr = 0.0D0
END SUBROUTINE zero_grad



SUBROUTINE zero_grad_v(r, p, p_fixed, dydr)
!############################################################
! f(x) = 0
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL zero_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE zero_grad_v