

SUBROUTINE heaviside(r, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
IF(r < 0.0D0)THEN
  y = 0.0D0
ELSE
  y = 1.0D0
END IF
END SUBROUTINE heaviside



SUBROUTINE heaviside_v(r, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  IF(r(n) < 0.0D0)THEN
    y(n) = 0.0D0
  ELSE
    y(n) = 1.0D0
  END IF
END DO
END SUBROUTINE heaviside_v

