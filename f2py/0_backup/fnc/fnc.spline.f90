! POWER SPLINE
! sum (ai (r - ri)^3 H(ri - r))     - cubic example
! sum (ai (r - ri)^5 H(ri - r))     - quintic example

! SCALAR SUBROUTINE
SUBROUTINE spline(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: rk
INTEGER(kind=StandardInteger) :: n, k
! Settings
INTEGER(kind=StandardInteger) :: power
INTEGER(kind=StandardInteger) :: mult
REAL(kind=DoubleReal) :: qa, qb
INTEGER(kind=StandardInteger) :: power_upper
REAL(kind=DoubleReal) :: rc_lower_h = 0.0D0
REAL(kind=DoubleReal) :: psi = 1.0D0
REAL(kind=DoubleReal) :: ypoly
!REAL(kind=DoubleReal) :: cutoff(1:5)
!############################################################
psi = 1.0D0
y = 0.0D0
ss = 10
ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)

IF(pfs .NE. (ps + ss))THEN
  RETURN
END IF

settings(1:ss) = p_fixed(1:ss)
nodes(1:ps) = p_fixed(ss+1:pfs)

power = INT(settings(1))    ! 3, 5 etc
mult = INT(settings(2))   ! 1 = (r - rk)^n   -1 = (rk - r)^n 
qa = settings(3)
qb = settings(4)
power_upper = INT(settings(5)) 
rc_lower_h = settings(6)
!cutoff(1:5) = settings(6:10)

IF(qa .GT. 0.0D0 .AND. qb .GT. 0.0D0)THEN
  IF(rc_lower_h .GT. 0.0D0)THEN
    psi = (r/rc_lower_h)**4
  END IF
  CALL fzbl(r, qa, qb, y)  
END IF

ypoly = 0.0D0
IF(power .GT. 0 .AND. power_upper .LE. power)THEN
  DO n = 1, SIZE(p,1)
    rk = nodes(n)
    IF((rk - r) .GE. 0.0D0)THEN
      ypoly = ypoly + p(n) * (mult * (r - rk))**power
    END IF
  END DO
ELSE IF(power .GT. 0 .AND. power_upper .GT. power)THEN
  k = power
  DO n = 1, SIZE(p,1)
    rk = nodes(n)
    IF((rk - r) .GE. 0.0D0)THEN
      ypoly = ypoly + p(n) * (mult * (r - rk))**k
    END IF
    IF(k .LT. power_upper)THEN
      k = k + 1
    ELSE
      k = power
    END IF
  END DO
END IF
y = y + psi * ypoly
END SUBROUTINE spline

! VECTOR SUBROUTINE
SUBROUTINE spline_v(r, p, p_fixed, y)
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
  CALL spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE spline_v
