!# Bonny Pesianot Terentyev Malerba Density

SUBROUTINE bonny_density(r, p, p_fixed, y)
!############################################################
! 
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: rc
REAL(kind=DoubleReal) :: rho0, beta
REAL(kind=DoubleReal) :: fcut
!############################################################
rc = p_fixed(1)
rho0 = p(1)
beta = p(2)

IF(r .GT. rc)THEN
  y = 0.0D0
ELSE 
  IF(r .LE. 1)THEN
    fcut = 1.0D0
  ELSE
    fcut = 1 - ((r - 1.0D0)**3 / (rc - 1)**3)
  END IF
  y = fcut * rho0 * (exp(-1.0D0 * beta * r) / r - exp(-1.0D0 * beta * rc) / rc)
END IF
END SUBROUTINE bonny_density


SUBROUTINE bonny_density_v(r, p, p_fixed, y)
!############################################################
! 
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL bonny_density(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE bonny_density_v