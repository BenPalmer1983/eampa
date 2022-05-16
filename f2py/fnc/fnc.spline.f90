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
INTEGER(kind=StandardInteger) :: power_lower
INTEGER(kind=StandardInteger) :: power_upper
INTEGER(kind=StandardInteger) :: rkoption
REAL(kind=DoubleReal) :: zbl_qa
REAL(kind=DoubleReal) :: zbl_qb
REAL(kind=DoubleReal) :: h_node_r = 0.0D0
REAL(kind=DoubleReal) :: h_node_c = 0.0D0
INTEGER(kind=StandardInteger) :: h_node_p
REAL(kind=DoubleReal) :: yzbl
!############################################################
y = 0.0D0
ss = 10
ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)

IF(pfs .LT. ss)THEN
  RETURN
END IF

!Read Settings
settings(1:ss) = p_fixed(1:ss)
power_lower = INT(settings(1))    ! 3, 5 etc
power_upper = INT(settings(2))    !
rkoption = INT(settings(3)) ! 0 = (rk - r)^n    1 = (r - rk)^n   
zbl_qa = settings(4)
zbl_qb = settings(5)
h_node_r = settings(6)    ! Hardening node rcutoff
h_node_c = settings(7)    ! Hardening node coefficient
h_node_p = INT(settings(8))    ! Hardening node power


IF(pfs .NE. (ps + ss))THEN
  RETURN
END IF
nodes(1:ps) = p_fixed(ss+1:pfs)

IF(r .GT. nodes(ps))THEN
  RETURN
END IF

! Hardening node (See G. Bonny 2012 code)
IF(h_node_r .GT. 0.0D0 .AND. r .LE. h_node_r)Then
  y = y + h_node_c * (h_node_r - r)**h_node_p
END IF

! Power spline
k = power_lower
DO n = 1, ps
  rk = nodes(n)
  IF((rk - r) .GE. 0.0D0)THEN
    IF(rkoption .EQ. 0)THEN
      y = y + p(n) * (rk - r)**k
    ELSE
      y = y + p(n) * (r - rk)**k
    END IF
  END IF
  IF(k .LT. power_upper)THEN
    k = k + 1
  ELSE
    k = power_lower
  END IF
END DO

! Add ZBL
IF(zbl_qa .GT. 0.0D0 .AND. zbl_qb .GT. 0.0D0)THEN
  CALL fzbl(r, zbl_qa, zbl_qb, yzbl)
  y = y + yzbl
END IF






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
