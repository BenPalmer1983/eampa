!# Olsson/Walenius
!# Ackland Mendelev 2004

SUBROUTINE triple_embedding(r, p, p_fixed, y)
!############################################################
! f(x) = -sqrt(r) + A * r**2 + B * r**4
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
!############################################################
IF(r .LT. 0.0D0)THEN
  y = 0.0D0
ELSE
  y =  p(1) * sqrt(r) + p(2) * (r)**2 + p(3) * (r)**4
END IF
END SUBROUTINE triple_embedding



SUBROUTINE triple_embedding_v(r, p, p_fixed, y)
!############################################################
! f(x) = A * sqrt(r) + B * r**2 + C * r**4
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL triple_embedding(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE triple_embedding_v
