! Two Band Modelling Fe-Cr Olsson, Wallenius
! CUBIC SPLINE
! sum (ai (r - ri)^3 H(ri - r))

! SCALAR SUBROUTINE
SUBROUTINE ackland_spline(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: settings(1:10)
REAL(kind=DoubleReal) :: nodes(1:SIZE(p,1))
REAL(kind=DoubleReal) :: rk, xi
INTEGER(kind=StandardInteger) :: n
! Settings
REAL(kind=DoubleReal) :: a, qa, qb, r1, r2, b0, b1, b2, b3
!############################################################

y = 0.0D0
ss = 10
ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)

IF(pfs .NE. (ps + ss))THEN
  RETURN
END IF

settings(1:ss) = p_fixed(1:ss)
nodes(1:ps) = p_fixed(ss+1:pfs)
r1 = settings(1)
r2 = settings(2)
a = settings(3)
qa = settings(4)
qb = settings(5)
b0 = settings(6)
b1 = settings(7)
b2 = settings(8)
b3 = settings(9)

IF(r .LE. r1)THEN
  CALL fzbl(r, qa, qb, xi)
  y = (a * xi) / r
ELSE IF(r .GT. r1 .AND. r .LT. r2)THEN
  y = EXP(b0 + b1 * r + b2 * r**2 + b3 * r**3)
ELSE
  y = 0.0D0
  DO n = 1, SIZE(p,1)
    rk = nodes(n)
    IF((rk - r) .GE. 0.0D0)THEN
      y = y + p(n) * (rk - r)**3
    END IF
  END DO
END IF



END SUBROUTINE ackland_spline

! VECTOR SUBROUTINE
SUBROUTINE ackland_spline_v(r, p, p_fixed, y)
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
  CALL ackland_spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE ackland_spline_v
