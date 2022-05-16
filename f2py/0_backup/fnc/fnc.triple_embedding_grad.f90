

SUBROUTINE triple_embedding_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = A * sqrt(r) + B * r**2 + C * r**4
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
IF(r .LT. 0.0D0)THEN
  dydr = 0.0D0
ELSE
  dydr = 0.5D0 * p(1) * r**(-0.5D0) + 2.0D0 * p(2) * (r) + 4.0D0 * p(3) * (r)**3
END IF
END SUBROUTINE triple_embedding_grad


SUBROUTINE triple_embedding_grad_v(r, p, p_fixed, dydr)
!############################################################
! f(x) = A * sqrt(r) + B * r**2 + C * r**4
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL triple_embedding_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE triple_embedding_grad_v

