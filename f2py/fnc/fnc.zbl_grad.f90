

SUBROUTINE zbl_grad(r, p, p_fixed, dydr)
!############################################################
! ZBL
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
CALL zbl(r-h, p, p_fixed, a)
CALL zbl(r+h, p, p_fixed, b)
dydr = (b - a) / (2.0D0 * h)
END SUBROUTINE zbl_grad



SUBROUTINE zbl_grad_v(r, p, p_fixed, dydr)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL zbl_grad(r(n), p,  p_fixed, dydr(n))
END DO
END SUBROUTINE zbl_grad_v


