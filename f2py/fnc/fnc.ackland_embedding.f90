!# Olsson/Walenius
!# Ackland Mendelev 2004

SUBROUTINE ackland_embedding(r, p, p_fixed, y)
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
  y = -1.0D0 * sqrt(r) + p(1) * (r)**2 + p(2) * (r)**4
END IF
END SUBROUTINE ackland_embedding



SUBROUTINE ackland_embedding_v(r, p, p_fixed, y)
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
  CALL ackland_embedding(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE ackland_embedding_v
