!# Mishin

SUBROUTINE mishin_density(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:7)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: rc, r0, h
REAL(kind=DoubleReal) :: A, B, C, gamma, yb
REAL(kind=DoubleReal) :: psi, e, x, z
!############################################################
rc = p_fixed(1)
r0 = p(1)
h = p(2)
A = p(3)
B = p(4)
C = p(5)
gamma = p(6)
yb = p(7)

IF(r .GT. rc)THEN
  y = 0.0D0
ELSE
  x = ((r - rc) / h)**4
  psi = x / (1 + x)
  z = r - r0
  e = exp(-1.0D0 * gamma * z)
  y = psi * (A * z**yb * e * (1 + B * e) + C)
END IF
END SUBROUTINE mishin_density


SUBROUTINE mishin_density_v(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL mishin_density(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE mishin_density_v
