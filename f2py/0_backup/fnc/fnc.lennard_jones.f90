

SUBROUTINE lennard_jones(r, p, p_fixed, y)
!############################################################
! LENNARD JONES FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
y = p(1) * ((p(2) / r)**12 - 2 * (p(2)/r)**6)
END SUBROUTINE lennard_jones

SUBROUTINE lennard_jones_v(r, p, p_fixed, y)
!############################################################
! LENNARD JONES VECTOR FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL lennard_jones(r(n), p, p_fixed,  y(n))
END DO
END SUBROUTINE lennard_jones_v





