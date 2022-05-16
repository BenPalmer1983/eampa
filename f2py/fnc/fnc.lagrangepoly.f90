SUBROUTINE lagrangepoly(r, p, p_fixed, y)
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
REAL(kind=DoubleReal) :: settings(1:3)
REAL(kind=DoubleReal) :: nodes(1:SIZE(p,1))
REAL(kind=DoubleReal) :: fcut
INTEGER(kind=StandardInteger) :: n


REAL(kind=DoubleReal) :: rcut1, rcut2
INTEGER(kind=StandardInteger) :: n_interp
!############################################################

y = 0.0D0

ss = 3
ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)
settings(1:ss) = p_fixed(1:ss)

rcut1 = settings(1)
rcut2 = settings(2)
n_interp = INT(settings(3))

IF(pfs .NE. (ps + ss))THEN
  RETURN
END IF

! Compute Cutoff
IF(r .LE. rcut1)THEN
  fcut = 1.0D0
ELSE IF(r .GT. rcut2)THEN
  fcut = 0.0D0
ELSE
	fcut = MAX(0.0D0, 1.0D0 - (r - rcut1)**3 / (rcut2 - rcut1)**3)
END IF


! Set Nodes
nodes(1:ps) = p_fixed(ss+1:pfs)

IF(r .LT. nodes(1))THEN
  n = 1
ELSE IF(r .GT. nodes(ps))THEN
  n = ps
ELSE
  DO n = 1, ps-1
    IF(r .EQ. nodes(n))THEN
      y = fcut * p(n)
      RETURN
		ELSE IF(r .EQ. nodes(n+1))THEN
		  y = fcut * p(n+1)
	    RETURN
    ELSE IF(r .GE. nodes(n) .AND. r .LE. nodes(n+1))THEN
      EXIT
    END IF
  END DO
END IF

n = MAX(1, MIN(n - n_interp / 2, ps - n_interp + 1))
CALL interpn(r, nodes(n:n+n_interp-1), p(n:n+n_interp-1), y)
y = y * fcut



!############################################################
END SUBROUTINE lagrangepoly


SUBROUTINE lagrangepoly_v(r, p, p_fixed, y)
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
  CALL lagrangepoly(r(n), p,  p_fixed, y(n))
END DO
END SUBROUTINE lagrangepoly_v





