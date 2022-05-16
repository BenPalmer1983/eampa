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
INTEGER(kind=StandardInteger) :: rkoption
REAL(kind=DoubleReal) :: rkmult
REAL(kind=DoubleReal) :: qa, qb
INTEGER(kind=StandardInteger) :: power_upper
REAL(kind=DoubleReal) :: rc_lower_h = 0.0D0
REAL(kind=DoubleReal) :: rc_lower_rcut = 0.0D0
REAL(kind=DoubleReal) :: knot_spline_rc1 = 0.0D0
REAL(kind=DoubleReal) :: knot_spline_rc2 = 0.0D0
INTEGER(kind=StandardInteger) :: nodes_in_p
REAL(kind=DoubleReal) :: x
REAL(kind=DoubleReal) :: psi = 1.0D0
REAL(kind=DoubleReal) :: ypoly
REAL(kind=DoubleReal) :: node_a(1:4)
REAL(kind=DoubleReal) :: node_b(1:4)
REAL(kind=DoubleReal) :: spline_arr(1:4)
!REAL(kind=DoubleReal) :: cutoff(1:5)
!############################################################
psi = 1.0D0
y = 0.0D0
ss = 10
ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)

IF(pfs .LT. ss)THEN
  RETURN
END IF

!Read Settings
settings(1:ss) = p_fixed(1:ss)
power = INT(settings(1))    ! 3, 5 etc
rkoption = INT(settings(2)) ! 0 = (rk - r)^n    1 = (r - rk)^n   
qa = settings(3)
qb = settings(4)
power_upper = INT(settings(5)) 
rc_lower_rcut = settings(6)
rc_lower_h = settings(7)
knot_spline_rc1 = settings(8) 
knot_spline_rc2 = settings(9)
nodes_in_p = INT(settings(10))   ! 0 = nodes are in p_fixed, 1 = nodes are variable and in p

! If the node x positions are fixed, and contained in the
IF(nodes_in_p .EQ. 0)THEN
  IF(pfs .NE. (ps + ss))THEN
    RETURN
  END IF
  nodes(1:ps) = p_fixed(ss+1:pfs)
ELSE
  ! change ps size to half
  ps = ps / 2
  nodes(1:ps) = p(ps+1:2*ps)
END IF

IF(rkoption .EQ. 0)THEN
  rkmult = 1.0D0
ELSE
  rkmult = -1.0D0
END IF


! Knot from spline to ZBL - must have Qa Qb and the two cutoffs
IF(knot_spline_rc1 .GT. 0.0D0 .AND. knot_spline_rc2 .GT. knot_spline_rc1 &
    .AND. qa .GT. 0.0D0 .AND. qb .GT. 0.0D0 )THEN
  IF(r .GT. 0.0 .AND. r .LE. knot_spline_rc1)THEN
    CALL fzbl(r, qa, qb, y)  
  ELSE IF(r .GT. knot_spline_rc1 .AND. r .LT. knot_spline_rc2)THEN
    !# Node A
    node_a(1) = knot_spline_rc1
    CALL fzbl(knot_spline_rc1, qa, qb, node_a(2))
    CALL fzbl_grad(knot_spline_rc1, qa, qb, node_a(3))
    CALL fzbl_ggrad(knot_spline_rc1, qa, qb, node_a(4))

    !# Node B
    node_b = 0.0D0
    node_b(1) = knot_spline_rc2
    !######################### Spline Start
    IF(power .GT. 0 .AND. power_upper .LE. power)THEN
      DO n = 1, ps
        rk = nodes(n)
        IF((rk - node_b(1)) .GE. 0.0D0)THEN
          node_b(2) = node_b(2) + p(n) * (rkmult * (rk - node_b(1)))**power
          node_b(3) = node_b(3) + power * p(n) * (rkmult * (rk - node_b(1)))**(power-1)
          node_b(4) = node_b(4) + power * (power - 1) * p(n) * (rkmult * (rk - node_b(1)))**(power-2)
        END IF
      END DO
    ELSE IF(power .GT. 0 .AND. power_upper .GT. power)THEN
      k = power
      DO n = 1, ps
        rk = nodes(n)
        IF((rk - node_b(1)) .GE. 0.0D0)THEN
          node_b(2) = node_b(2) + p(n) * (rkmult * (rk - node_b(1)))**k
          node_b(3) = node_b(3) + k * p(n) * (rkmult * (rk - node_b(1)))**(k-1)
          node_b(4) = node_b(4) + k * (k-1) * p(n) * (rkmult * (rk - node_b(1)))**(k-2)
        END IF
        IF(k .LT. power_upper)THEN
          k = k + 1
        ELSE
          k = power
        END IF
      END DO
    END IF

    !# Knot Spline
    CALL spline_nodes(r, 5, node_a, node_b, spline_arr)
    y = spline_arr(2)
    !print *, node_a(:)
    !print *, node_b(:)
    !######################### Spline End  
  ELSE IF(r .GE. knot_spline_rc2 .AND. r .LT. nodes(ps))THEN
    !######################### Spline Start
    IF(power .GT. 0 .AND. power_upper .LE. power)THEN
      DO n = 1, ps
        rk = nodes(n)
        IF((rk - r) .GE. 0.0D0)THEN
          y = y + p(n) * (rkmult * (rk - r))**power
        END IF
      END DO
    ELSE IF(power .GT. 0 .AND. power_upper .GT. power)THEN
      k = power
      DO n = 1, ps
        rk = nodes(n)
        IF((rk - r) .GE. 0.0D0)THEN
          y = y + p(n) * (rkmult * (rk - r))**k
        END IF
        IF(k .LT. power_upper)THEN
          k = k + 1
        ELSE
          k = power
        END IF
      END DO
    END IF
    !######################### Spline End
  ELSE
    y = 0.0D0
  END IF

! No knot from spline to ZBL
ELSE
  ! ZBL
  IF(qa .GT. 0.0D0 .AND. qb .GT. 0.0D0)THEN
    CALL fzbl(r, qa, qb, y)  
  END IF

  ! Force spline to zero at 0
  IF(rc_lower_h .GT. 0.0D0 .AND. rc_lower_rcut .GE. 0.0D0)THEN
    IF(r .LT. rc_lower_rcut)THEN
      psi = 0.0D0
    ELSE
      x = ((r - rc_lower_rcut) / rc_lower_h)**4
      psi = x / (1 + x)
    END IF
  END IF

  ypoly = 0.0D0
  !######################### Spline Start
  IF(power .GT. 0 .AND. power_upper .LE. power)THEN
    DO n = 1, ps
      rk = nodes(n)
      IF((rk - r) .GE. 0.0D0)THEN
        ypoly = ypoly + p(n) * (rkmult * (rk - r))**power
      END IF
    END DO
  ELSE IF(power .GT. 0 .AND. power_upper .GT. power)THEN
    k = power
    DO n = 1, ps
      rk = nodes(n)
      IF((rk - r) .GE. 0.0D0)THEN
        ypoly = ypoly + p(n) * (rkmult * (rk - r))**k
      END IF
      IF(k .LT. power_upper)THEN
        k = k + 1
      ELSE
        k = power
      END IF
    END DO
  END IF
  !######################### Spline End
  y = y + psi * ypoly
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
