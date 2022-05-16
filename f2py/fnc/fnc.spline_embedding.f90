! POWER SPLINE
! sum (ai (r - ri)^3 H(ri - r))     - cubic example
! sum (ai (r - ri)^5 H(ri - r))     - quintic example

! SCALAR SUBROUTINE
SUBROUTINE spline_embedding(r, p, p_fixed, y)
!############################################################
! p coefficients
! pf r cutoffs
! they must be the same size
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
INTEGER(kind=StandardInteger) :: ps, pfs, ss
REAL(kind=DoubleReal) :: rk
INTEGER(kind=StandardInteger) :: n
!############################################################

y = 0.0D0
ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)

IF(ps .NE. pfs + 1)THEN
  RETURN
END IF

y = p(1) * sqrt(r)

DO n = 2, SIZE(p,1)
  rk = p_fixed(n-1)
  IF(r .GT. rk)THEN
    y = y + p(n) * (r - rk)**4
  END IF
END DO


END SUBROUTINE spline_embedding

! VECTOR SUBROUTINE
SUBROUTINE spline_embedding_v(r, p, p_fixed, y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
!############################################################
! Loop through all the values in r(:), calculate and store in y(:)
DO n = 1, SIZE(r,1)
  CALL spline_embedding(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE spline_embedding_v