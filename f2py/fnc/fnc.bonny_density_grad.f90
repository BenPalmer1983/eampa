! ### GRADIENT ###
! Bonny Density


SUBROUTINE bonny_density_grad(r, p, p_fixed, dydr)
!############################################################
! 
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: h = 1.0D-7
REAL(kind=DoubleReal) :: a, b
!############################################################
CALL bonny_density(r-h, p, p_fixed, a)
CALL bonny_density(r+h, p, p_fixed, b)
dydr = (b - a) / (2.0D0 * h)
END SUBROUTINE bonny_density_grad


SUBROUTINE bonny_density_grad_v(r, p, p_fixed, dydr)
!############################################################
! Bonny Density
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL bonny_density_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE bonny_density_grad_v