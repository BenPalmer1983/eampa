!# Polynomial embedding (GRAD)
!# General function for embedding functions
!# like F(r) = sqrt(r)
!#      F(r) = A sqrt(r) + B r**2 + C r**3


SUBROUTINE polynomial_embedding_grad(r, p, p_fixed, y)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)           !# Coeff
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)     !# Power
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0

IF(SIZE(p,1) .NE. SIZE(p_fixed,1))THEN
  RETURN
END IF

DO n = 1, SIZE(p,1)
  y = y + (p(n) * p_fixed(n)) * r**(p_fixed(n) - 1.0D0)
END DO

END SUBROUTINE polynomial_embedding_grad



SUBROUTINE polynomial_embedding_grad_v(r, p, p_fixed, y)
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
  CALL polynomial_embedding_grad(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE polynomial_embedding_grad_v

