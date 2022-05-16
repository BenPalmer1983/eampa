

SUBROUTINE mishin_density_grad(r, p, p_fixed, dydr)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
!############################################################
IF(r .LT. 0.0D0)THEN
  dydr = 0.0D0
ELSE
  dydr = -0.5D0 * r**(-0.5D0) + 2 * p(1) * r
END IF
END SUBROUTINE mishin_density_grad



SUBROUTINE mishin_density_grad_v(r, p, p_fixed, dydr)
!############################################################
! BUCKINGHAM POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL mishin_density_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE mishin_density_grad_v


