SUBROUTINE step(r, h, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: h
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
IF(r < h)THEN
  y = 0.0D0
ELSE
  y = 1.0D0
END IF
END SUBROUTINE step


SUBROUTINE step_v(r, h, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: h
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  IF(r(n) < h)THEN
    y(n) = 0.0D0
  ELSE
    y(n) = 1.0D0
  END IF
END DO
END SUBROUTINE step_v
