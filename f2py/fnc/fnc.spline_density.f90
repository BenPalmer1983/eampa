SUBROUTINE spline_density(r, p, p_fixed, y)
!############################################################
! Single "spline" by Bonny et al with additional cutoff approaching 0
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: c0
REAL(kind=DoubleReal) :: rc1
REAL(kind=DoubleReal) :: rc2
INTEGER(kind=StandardInteger) :: rc1_power
REAL(kind=DoubleReal) :: rc1mult
!############################################################

y = 0.0D0

IF(SIZE(p,1) .NE. 1 .OR. SIZE(p_fixed,1) .NE. 3)THEN
  RETURN
END IF

c0 = p(1)
rc1 = p_fixed(1)
rc2 = p_fixed(2)
rc1_power = INT(p_fixed(3))


IF(r .LT. 0.0D0 .OR. r .GT. rc2)THEN
  RETURN
END IF

rc1mult = 1.0D0
IF(r .LT. rc1 .AND. rc1 .GT. 0.0D0)THEN
  rc1mult = 1.0D0 - (rc1 - r)**rc1_power / (rc1)**rc1_power
END IF
y = c0 * rc1mult * (rc2 - r)**3


!############################################################
END SUBROUTINE spline_density

SUBROUTINE spline_density_v(r, p, p_fixed, y)
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
  CALL spline_density(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE spline_density_v