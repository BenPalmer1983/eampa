! ### GRADIENT ###
! Slater 4S density


SUBROUTINE slater_4s_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: h = 1.0D-7
REAL(kind=DoubleReal) :: a, b
!############################################################
CALL slater_4s(r-h, p, p_fixed, a)
CALL slater_4s(r+h, p, p_fixed, b)
dydr = (b - a) / (2.0D0 * h)
END SUBROUTINE slater_4s_grad


SUBROUTINE slater_4s_grad_v(r, p, p_fixed, dydr)
!############################################################
! BUCKINGHAM POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL slater_4s_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE slater_4s_grad_v




