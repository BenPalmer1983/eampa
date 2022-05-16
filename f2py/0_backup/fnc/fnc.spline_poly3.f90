

!############################################################
! Fix the end nodes with p_fix (x, y, y', y'' for start and end)


SUBROUTINE spline_one_node(r, p_var, p_fix, y)
!############################################################
! MORSE FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p_var(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: p_fix(1:8)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
INTEGER(kind=StandardInteger) :: c
REAL(kind=DoubleReal) :: coeffs(1:6)
!############################################################
IF(r .LE. p_fix(1))THEN
  y = p_fix(2)
ELSE IF(r .EQ. p_var(1))THEN
  y = p_var(2)
ELSE IF(r .GE. p_fix(5))THEN
  y = p_fix(6)
ELSE IF(r .GT. p_fix(1) .AND. r .LT. p_var(1))THEN
  CALL spline_ab(p_fix(1:4), p_var(1:4), coeffs(1:6))
  y = 0.0D0
  DO c = 1, 6
    y = y + coeffs(c) * r**(c-1)
  END DO
ELSE IF(r .GT. p_var(1) .AND. r .LT. p_fix(5))THEN
  CALL spline_ab(p_var(1:4), p_fix(5:8), coeffs(1:6))
  y = 0.0D0
  DO c = 1, 6
    y = y + coeffs(c) * r**(c-1)
  END DO
END IF
END SUBROUTINE spline_poly3


SUBROUTINE spline_poly3_v(r, p_var, p_fix, y)
!############################################################
! MORSE FUNCTION
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_var(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fix(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL spline_poly3(r(n), p_var, p_fix, y(n))
END DO
END SUBROUTINE spline_poly3_v
