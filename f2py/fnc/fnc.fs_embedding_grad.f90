SUBROUTINE fs_embedding_grad(r, p, p_fixed, dydr)
!############################################################
! f(x) = (N r^3 e^(-z r))^2
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
  dydr = -0.5D0 * p(1) * r**(-0.5D0)
END IF
END SUBROUTINE fs_embedding_grad



SUBROUTINE fs_embedding_grad_v(r, p, p_fixed, dydr)
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
  CALL fs_embedding_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE fs_embedding_grad_v

