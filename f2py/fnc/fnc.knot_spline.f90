SUBROUTINE knot_spline(r, p, p_fixed, y)
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
! Sizes
INTEGER(kind=StandardInteger) :: ps, pns, pfs, ns, spline_option
! Settings
INTEGER(kind=StandardInteger) :: ss = 20
REAL(kind=DoubleReal) :: settings(1:20)
INTEGER(kind=StandardInteger) :: poly = 3
INTEGER(kind=StandardInteger) :: interp_size = 5
INTEGER(kind=StandardInteger) :: node_size = 0
INTEGER(kind=StandardInteger) :: fix_start_node = 0
REAL(kind=DoubleReal) :: start_x = 0
REAL(kind=DoubleReal) :: start_y = 0
REAL(kind=DoubleReal) :: start_dydx = 0
REAL(kind=DoubleReal) :: start_d2ydx2 = 0
INTEGER(kind=StandardInteger) :: fix_end_node = 0
REAL(kind=DoubleReal) :: end_x = 0
REAL(kind=DoubleReal) :: end_y = 0
REAL(kind=DoubleReal) :: end_dydx = 0
REAL(kind=DoubleReal) :: end_d2ydx2 = 0
INTEGER(kind=StandardInteger) :: zbl_on = 0
REAL(kind=DoubleReal) :: ra = 0.0D0
REAL(kind=DoubleReal) :: rb = 0.0D0
REAL(kind=DoubleReal) :: za = 0.0D0
REAL(kind=DoubleReal) :: zb = 0.0D0
!
INTEGER(kind=StandardInteger) :: offset
INTEGER(kind=StandardInteger) :: fill_d = 0
REAL(kind=DoubleReal) :: nodes(1:SIZE(p,1)+10,1:4)
REAL(kind=DoubleReal) :: temp(1:SIZE(p,1)+10)
REAL(kind=DoubleReal) :: node_min_x
REAL(kind=DoubleReal) :: rspline
LOGICAL :: make_rb_a_node = .FALSE.
LOGICAL :: compute_spline = .FALSE.
LOGICAL :: zbl_spline = .FALSE.

REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb
INTEGER(kind=StandardInteger) :: i, j, k, n
REAL(kind=DoubleReal) :: xyarr(1:4)
REAL(kind=DoubleReal) :: na(1:10) = 0.0D0
REAL(kind=DoubleReal) :: nb(1:10) = 0.0D0
REAL(kind=DoubleReal) :: xmat(1:10,1:10)
REAL(kind=DoubleReal) :: ymat(1:10)
REAL(kind=DoubleReal) :: c(1:10)
REAL(kind=DoubleReal) :: rb_node(1:4) = 0.0D0
!############################################################

ps = SIZE(p,1)
pfs = SIZE(p_fixed,1)

! zero y default and load settings
y = 0.0D0
settings(1:ss) = p_fixed(1:ss)

IF(pfs .LE. ss)THEN
  y = 0.0D0  
  RETURN
END IF

IF(pfs .EQ. ps + ss)THEN
  pns = ps
  spline_option = 1
ELSE IF(pfs .EQ. (ps / 2) + ss)THEN
  pns = ps / 2
  spline_option = 2
ELSE IF(pfs .EQ. (ps / 3) + ss)THEN
  pns = ps / 3
  spline_option = 3
ELSE
	y = 0.0D0  
	RETURN
END IF


! settings(1)      3 = cubic
! settings(2)      2 = derivative interpolation points - default is 5
! settings(3)      0 = None, n nodes using interpolation
! settings(4)      1 = fix start node
! settings(5)      0 = start node x
! settings(6)      0 = start node y
! settings(7)      0 = start node dydx
! settings(8)      0 = start node dy2dy2
! settings(9)      1 = fix end node
! settings(10)     0 = end node x
! settings(11)     0 = end node y
! settings(12)     0 = end node dydx
! settings(13)     0 = end node dy2dy2
! settings(14)     0 = zbl off, zbl on
! settings(15)     0 = ra
! settings(16)     0 = rb
! settings(17)     0 = za
! settings(18)     0 = zb


IF(settings(1) .EQ. 3 .OR. settings(1) .EQ. 5) poly = settings(1)
IF(settings(2) .GE. 3 .AND. settings(2) .LE. 6) interp_size = settings(2)
IF(settings(3) .GE. 0 .AND. settings(3) .LE. 30) node_size = settings(3)

fix_start_node = INT(settings(4))
start_x = settings(5)
start_y = settings(6)
start_dydx = settings(7)
start_d2ydx2 = settings(8)

fix_end_node = INT(settings(9))
end_x = settings(10)
end_y = settings(11)
end_dydx = settings(12)
end_d2ydx2 = settings(13)

zbl_on = INT(settings(14))
ra = settings(15)
rb = settings(16)
za = settings(17)
zb = settings(18)

! IF ZBL AND IN ZBL PART
IF((zbl_on .eq. 1) .AND. (r .LE. ra))THEN
  IF(r .LT. 0.0D0)THEN
    y = 0.0D0
    RETURN 
  END IF
  CALL fzbl(r, za, zb, y)
  RETURN 
END IF

node_min_x = p_fixed(ss+1)
rspline = r

! Is ZBL on and is this this spline between ZBL and the remaining nodes?
zbl_spline = .FALSE.
IF((zbl_on .eq. 1) .AND. (rspline .GT. ra) .AND. (rspline .LT. rb))THEN
  rspline = rb
  zbl_spline = .TRUE.
END IF

!
!   TO DO - MAKE RB from ZBL a NODE
!


! READ NODES AND COMPUTE DERIVATIVES WHERE NOT PROVIDE

! Zero array
nodes = 0.0D0
! If the zbl b spline node is less than the existing nodes, add it as a node
make_rb_a_node = .FALSE.
IF((zbl_on .eq. 1) .AND. (rb .LT. node_min_x) .AND. &
   (fix_start_node .EQ. 0 .OR. ((fix_start_node .EQ. 1) .AND. (rb .LT. start_x))))THEN
  make_rb_a_node = .TRUE.
END IF

! User provides x and f(x) ONLY
! Calculate first and second derivatives
IF (spline_option .EQ. 1) THEN  
  offset = 0
  ns = pns
  ! If zbl on and rb is less than the first node then add rb as first node
  IF(make_rb_a_node)THEN
    nodes(1+offset, 1) = rb
    nodes(1+offset, 2) = 0.0D0
    nodes(1+offset, 3) = 0.0D0
    nodes(1+offset, 4) = 0.0D0
    offset = offset + 1
    ns = ns + 1
  END IF
  IF(fix_start_node .EQ. 1)THEN
    nodes(1+offset, 1) = start_x
    nodes(1+offset, 2) = start_y
    offset = offset + 1
    ns = ns + 1
  END IF
  nodes(1+offset:pns+offset, 1) = p_fixed(ss+1:pns+ss)
  nodes(1+offset:pns+offset, 2) = p(1:pns)
  IF(fix_end_node .EQ. 1)THEN
    nodes(pns+1+offset, 1) = end_x
    nodes(pns+1+offset, 2) = end_y
    ns = ns + 1
  END IF
  ! Fill in f(rb) by interpolating the existing nodes
  IF(make_rb_a_node)THEN
    CALL interpn(nodes(1, 1), nodes(2:6, 1), nodes(2:6, 2), nodes(1, 2))
  END IF
  ! Fill in the first and second derivatives of nodes
  CALL interp_fill_d(nodes(1:ns, 1), nodes(1:ns, 2), 5, nodes(1:ns, 3))
  CALL interp_fill_d(nodes(1:ns, 1), nodes(1:ns, 3), 5, nodes(1:ns, 4))
  ! Overwrite start/end node if fixed
  offset = 0  
  IF(make_rb_a_node)THEN
    offset = offset + 1    
  END IF
  IF(fix_start_node .EQ. 1)THEN
    nodes(1+offset, 3) = start_dydx
    nodes(1+offset, 4) = start_d2ydx2
    offset = 1
  END IF
  IF(fix_end_node .EQ. 1)THEN
    nodes(pns+1+offset, 3) = end_dydx
    nodes(pns+1+offset, 4) = end_d2ydx2
  END IF
! User provides x,f(x) and f'(x)
! Calculate second derivatives
ELSE IF (spline_option .EQ. 2) THEN
  offset = 0
  ns = pns
  IF(make_rb_a_node)THEN
    nodes(1+offset, 1) = rb
    nodes(1+offset, 2) = 0.0D0
    nodes(1+offset, 3) = 0.0D0
    nodes(1+offset, 4) = 0.0D0
    offset = offset + 1
    ns = ns + 1
  END IF
  IF(fix_start_node .EQ. 1)THEN
    nodes(1+offset, 1) = start_x
    nodes(1+offset, 2) = start_y
    nodes(1+offset, 3) = start_dydx
    offset = 1
    ns = ns + 1
  END IF
  nodes(1+offset:pns+offset, 1) = p_fixed(ss+1:pns+ss)
  nodes(1+offset:pns+offset, 2) = p(1:pns)
  nodes(1+offset:pns+offset, 3) = p(pns+1:2*pns)
  IF(fix_end_node .EQ. 1)THEN 
    nodes(pns+1+offset, 1) = end_x
    nodes(pns+1+offset, 2) = end_y
    nodes(pns+1+offset, 3) = end_dydx
    ns = ns + 1
  END IF
  ! Fill in f(rb) and f'(rb) by interpolating the existing nodes
  IF(make_rb_a_node)THEN
    CALL interpn(nodes(1, 1), nodes(2:6, 1), nodes(2:6, 2), nodes(1, 2))
    CALL interpndydxn(nodes(1, 1), nodes(1:5, 1), nodes(1:5, 2), 1, nodes(1, 3))
  END IF
  ! Fill in the second derivatives of nodes
  CALL interp_fill_d(nodes(1:ns, 1), nodes(1:ns, 3), 5, nodes(1:ns, 4))
  offset = 0
  IF(make_rb_a_node)THEN
    offset = offset + 1    
  END IF
  IF(fix_start_node .EQ. 1)THEN
    nodes(1, 4) = start_d2ydx2
    offset = 1
  END IF
  IF(fix_end_node .EQ. 1)THEN
    nodes(pns+offset+1, 4) = end_d2ydx2
  END IF
! User provides x,f(x), f'(x) and f''(x)
ELSE IF (spline_option .EQ. 3) THEN
  offset = 0
  ns = pns
  IF(make_rb_a_node)THEN
    nodes(1+offset, 1) = rb
    nodes(1+offset, 2) = 0.0D0
    nodes(1+offset, 3) = 0.0D0
    nodes(1+offset, 4) = 0.0D0
    offset = offset + 1
    ns = ns + 1
  END IF
  IF(fix_start_node .EQ. 1)THEN
    nodes(1+offset, 1) = start_x
    nodes(1+offset, 2) = start_y
    nodes(1+offset, 3) = start_dydx
    nodes(1+offset, 4) = start_d2ydx2
    offset = 1
    ns = ns + 1
  END IF
  nodes(1+offset:pns+offset, 1) = p_fixed(ss+1:pns+ss)
  nodes(1+offset:pns+offset, 2) = p(1:pns)
  nodes(1+offset:pns+offset, 3) = p(pns+1:2*pns)
  nodes(1+offset:pns+offset, 4) = p(2*pns+1:3*pns)
  IF(fix_end_node .EQ. 1)THEN
    nodes(pns+1+offset, 1) = end_x
    nodes(pns+1+offset, 2) = end_y
    nodes(pns+1+offset, 3) = end_dydx
    nodes(pns+1+offset, 4) = end_d2ydx2
    ns = ns + 1
  END IF
  ! Fill in f(rb), f'(rb) and f''(rb) by interpolating the existing nodes
  IF(make_rb_a_node)THEN
    CALL interpn(nodes(1, 1), nodes(2:6, 1), nodes(2:6, 2), nodes(1, 2))
    CALL interpndydxn(nodes(1, 1), nodes(1:5, 1), nodes(1:5, 2), 1, nodes(1, 3))
    CALL interpndydxn(nodes(1, 1), nodes(1:5, 1), nodes(1:5, 3), 1, nodes(1, 4))
  END IF
ELSE
  RETURN
END IF

compute_spline = .TRUE.

! If out of range (or an end node)
IF(rspline .LE. nodes(1, 1))THEN
  y = nodes(1, 2)
  IF(.NOT. zbl_spline)THEN
    RETURN
  END IF
  compute_spline = .FALSE.
  xyarr(:) =  nodes(1, :)
ELSE IF(rspline .GE. nodes(ns, 1))THEN
  y = nodes(ns, 2)
  IF(.NOT. zbl_spline)THEN
    RETURN
  END IF
  compute_spline = .FALSE.
  xyarr(:) =  nodes(1, :)
END IF

! If exactly a node
DO n = 2, ns-1
  IF(rspline .EQ. nodes(n, 1))THEN
    y = nodes(n, 2)
    IF(.NOT. zbl_spline)THEN
      RETURN
    END IF
    compute_spline = .FALSE.
    xyarr(:) =  nodes(n, :)
  END IF
END DO


IF(compute_spline)THEN
  ! Find which nodes r is between
  n = 1
  DO WHILE(n .LT. ns-1)
    IF((nodes(n, 1) .LT. rspline) .AND. (rspline .LT. nodes(n+1, 1)))THEN
      EXIT
    END IF
    n = n + 1
  END DO

  na(1) = nodes(n, 1)
  na(2) = nodes(n, 2)
  na(3) = nodes(n, 3)
  na(4) = nodes(n, 4)

  nb(1) = nodes(n+1, 1)
  nb(2) = nodes(n+1, 2)
  nb(3) = nodes(n+1, 3)
  nb(4) = nodes(n+1, 4)

  CALL spline_nodes(rspline, poly, na(1:4), nb(1:4), xyarr)
  y = xyarr(2)
END IF


IF(zbl_spline)THEN
  rspline = r

  na(1) = ra
  CALL fzbl(ra, za, zb, na(2))
  CALL fzbl_grad(ra, za, zb, na(3))
  CALL fzbl_ggrad(ra, za, zb, na(4))

  nb(1) = xyarr(1)
  nb(2) = xyarr(2)
  nb(3) = xyarr(3)
  nb(4) = xyarr(4)

  CALL spline_nodes(rspline, poly, na, nb, xyarr)
  y = xyarr(2)
END IF

END SUBROUTINE knot_spline






! VECTOR SUBROUTINE
SUBROUTINE knot_spline_v(r, p, p_fixed, y)
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
  CALL knot_spline(r(n), p, p_fixed, y(n))
END DO
END SUBROUTINE knot_spline_v















SUBROUTINE spline_nodes1(rspline, poly, na, nb, xyarr)
!############################################################
! p coefficients
! pf r cutoffs
! they must be the same size
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: rspline
INTEGER(kind=StandardInteger), INTENT(IN) :: poly
REAL(kind=DoubleReal), INTENT(IN) :: na(:)
REAL(kind=DoubleReal), INTENT(IN) :: nb(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: xyarr(1:4)
!############################################################
REAL(kind=DoubleReal) :: xmat(1:10,1:10)
REAL(kind=DoubleReal) :: ymat(1:10)
REAL(kind=DoubleReal) :: c(1:10)
!############################################################
xyarr(:) = 0.0D0
IF(poly .EQ. 3)THEN
  ! X VALUES IN p_fixed
  ! Y VALUES IN p
  ! DERIV VALUES IN dydr  
  xmat(1,1) = 1.0
  xmat(1,2) = na(1)
  xmat(1,3) = na(1)**2
  xmat(1,4) = na(1)**3
  xmat(2,1) = 1.0
  xmat(2,2) = nb(1)
  xmat(2,3) = nb(1)**2
  xmat(2,4) = nb(1)**3
  xmat(3,1) = 0.0
  xmat(3,2) = 1.0D0
  xmat(3,3) = 2.0D0 * na(1)
  xmat(3,4) = 3.0D0 * na(1)**2
  xmat(4,1) = 0.0
  xmat(4,2) = 1.0D0
  xmat(4,3) = 2.0D0 * nb(1)
  xmat(4,4) = 3.0D0 * nb(1)**2
  
  ymat(1) = na(2)
  ymat(2) = nb(2)
  ymat(3) = na(3)
  ymat(4) = nb(3)
    
  CALL sls_solve(xmat(1:4,1:4), ymat(1:4), c(1:4))

  xyarr(1) = rspline
  xyarr(2) = c(1) + c(2) * rspline + c(3) * rspline**2 + c(4) * rspline**3
  xyarr(3) = c(2) + 2 * c(3) * rspline + 3 * c(4) * rspline**2
  xyarr(4) = 2 * c(3) + 6 * c(4) * rspline

ELSE IF(poly .EQ. 5)THEN
  ! X VALUES IN p_fixed
  ! Y VALUES IN p
  ! DERIV VALUES IN dydr  
  xmat(1,1) = 1.0D0
  xmat(1,2) = na(1)
  xmat(1,3) = na(1)**2
  xmat(1,4) = na(1)**3
  xmat(1,5) = na(1)**4
  xmat(1,6) = na(1)**5
  xmat(2,1) = 1.0D0
  xmat(2,2) = nb(1)
  xmat(2,3) = nb(1)**2
  xmat(2,4) = nb(1)**3
  xmat(2,5) = nb(1)**4
  xmat(2,6) = nb(1)**5
  xmat(3,1) = 0.0D0
  xmat(3,2) = 1.0D0
  xmat(3,3) = 2.0D0 * na(1)
  xmat(3,4) = 3.0D0 * na(1)**2
  xmat(3,5) = 4.0D0 * na(1)**3
  xmat(3,6) = 5.0D0 * na(1)**4
  xmat(4,1) = 0.0D0
  xmat(4,2) = 1.0D0
  xmat(4,3) = 2.0D0 * nb(1)
  xmat(4,4) = 3.0D0 * nb(1)**2  
  xmat(4,5) = 4.0D0 * nb(1)**3  
  xmat(4,6) = 5.0D0 * nb(1)**4  
  xmat(5,1) = 0.0D0
  xmat(5,2) = 0.0D0
  xmat(5,3) = 2.0D0 
  xmat(5,4) = 6.0D0 * na(1)
  xmat(5,5) = 12.0D0 * na(1)**2
  xmat(5,6) = 20.0D0 * na(1)**3
  xmat(6,1) = 0.0D0
  xmat(6,2) = 0.0D0
  xmat(6,3) = 2.0D0 
  xmat(6,4) = 6.0D0 * nb(1)
  xmat(6,5) = 12.0D0 * nb(1)**2
  xmat(6,6) = 20.0D0 * nb(1)**3
  
  ymat(1) = na(2)
  ymat(2) = nb(2)
  ymat(3) = na(3)
  ymat(4) = nb(3)
  ymat(5) = na(4)
  ymat(6) = nb(4)
    
  CALL sls_solve(xmat(1:6,1:6), ymat(1:6), c(1:6))

  xyarr(1) = rspline
  xyarr(2) = c(1) + c(2) * rspline + c(3) * rspline**2 + c(4) * rspline**3 + &
                 c(5) * rspline**4 + c(6) * rspline**5
  xyarr(3) = c(2) + 2 * c(3) * rspline + 3 * c(4) * rspline**2 + &
                 4 * c(5) * rspline**3 + 5 * c(6) * rspline**4
  xyarr(4) = 2 * c(3) + 6 * c(4) * rspline + &
                 12 * c(5) * rspline**2 + 20 * c(6) * rspline**3
END IF


END SUBROUTINE spline_nodes1



