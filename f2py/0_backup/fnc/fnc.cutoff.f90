! Cutoff

SUBROUTINE cutoff(r, rcut, n, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: rcut
INTEGER(kind=StandardInteger), INTENT(IN) :: n
REAL(kind=DoubleReal), INTENT(INOUT) :: y
!############################################################
IF(r .GT. rcut)THEN
  y = 0.0D0
ELSE
  y = y * ((rcut - r) / rcut)**n
END IF
END SUBROUTINE cutoff
