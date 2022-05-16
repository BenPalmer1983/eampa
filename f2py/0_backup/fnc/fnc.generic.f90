


FUNCTION f_zbl(r, q1, q2) RESULT (y)
!############################################################
! ZBL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: q1, q2, r
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: rs, x, e
!############################################################
rs = 0.4683766D0 / (q1**(2.0D0/3.0D0) + q2**(2.0D0/3.0D0))
x = (r / rs)
e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
y = ((q1 * q2) / r) * e
END FUNCTION f_zbl


FUNCTION f_zbl_ackland_mendelev(r, q1, q2) RESULT (y)
!############################################################
! ZBL Development of an interatomic potential for phosphorus impurities in α-iron 2004
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: q1, q2, r
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: rs, x, e
!############################################################
rs = 0.4683766D0 / (q1**(2.0D0/3.0D0) + q2**(2.0D0/3.0D0))
x = (r / rs)
e = 0.1818D0 * exp(-3.2 * x) + 0.5099D0 * exp(-0.9423 * x) + 0.2802D0 * exp(-0.4029 * x) + 0.02817D0 * exp(-0.2016 * x)
y = (q1 * q2 * e) / r
END FUNCTION f_zbl_ackland_mendelev


FUNCTION f_spline(r, p, pf, k) RESULT (y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: p(:)
REAL(kind=DoubleReal) :: pf(:)
INTEGER(kind=StandardInteger) :: k
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: H, rk
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
IF(SIZE(p,1) == SIZE(pf,1) .AND. SIZE(pf,1) > 0)THEN 
  DO n = 1, SIZE(p,1)
    rk = pf(n)
    H = f_heaviside(rk - r)
    y = y + H * p(n) * (rk-r)**3
  END DO
END IF
END FUNCTION f_spline


FUNCTION f_cubic_spline(r, p, pf) RESULT (y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: p(:)
REAL(kind=DoubleReal) :: pf(:)
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: H, rk
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
IF(SIZE(p,1) == SIZE(pf,1) .AND. SIZE(pf,1) > 0)THEN 
  DO n = 1, SIZE(p,1)
    rk = pf(n)
    H = f_heaviside(rk - r)
    y = y + H * p(n) * (rk-r)**3
  END DO
END IF
END FUNCTION f_cubic_spline


FUNCTION f_quintic_spline(r, p, pf) RESULT (y)
!############################################################
! CUBIC SPLINE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: p(:)
REAL(kind=DoubleReal) :: pf(:)
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: H, rk
INTEGER(kind=StandardInteger) :: n
!############################################################
y = 0.0D0
IF(SIZE(p,1) == SIZE(pf,1) .AND. SIZE(pf,1) > 0)THEN 
  DO n = 1, SIZE(p,1)
    rk = pf(n)
    H = f_heaviside(rk - r)
    y = y + H * p(n) * (rk-r)**5
  END DO
END IF
END FUNCTION f_quintic_spline



FUNCTION f_heaviside(r) RESULT (y)
!############################################################
! HEAVISIDE
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
INTEGER(kind=StandardInteger) :: y
!############################################################
IF(r < 0.0D0)THEN
  y = 0
ELSE
  y = 1
END IF
END FUNCTION f_heaviside



FUNCTION f_smootherstep_01(r) RESULT (y)
!############################################################
! https://en.wikipedia.org/wiki/Smoothstep
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: y
!############################################################
IF(r <= 0.0D0)THEN
  y = 0.0D0
ELSE IF(r >= 1.0D0)THEN
  y = 1.0D0
ELSE
  y = 6.0D0*r**5 - 15.0D0*r**4 + 10.0D0*r**3
END IF
END FUNCTION f_smootherstep_01



FUNCTION f_smootherstep_10(r) RESULT (y)
!############################################################
! https://en.wikipedia.org/wiki/Smoothstep
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: y
!############################################################
IF(r <= 0.0D0)THEN
  y = 1.0D0
ELSE IF(r >= 1.0D0)THEN
  y = 0.0D0
ELSE
  y = 1.0D0 - (6.0D0*r**5 - 15.0D0*r**4 + 10.0D0*r**3)
END IF
END FUNCTION f_smootherstep_10




FUNCTION f_cubic_knot_ab(r, xa, ya, ypa, xb, yb, ypb) RESULT (y)
! CUBIC KNOT
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: xa, ya, ypa, xb, yb, ypb
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: xmat(1:4,1:4)
REAL(kind=DoubleReal) :: xmat_inv(1:4,1:4)
REAL(kind=DoubleReal) :: ymat(1:4)
REAL(kind=DoubleReal) :: c(1:4)
!############################################################
xmat(1,1) = 1.0
xmat(1,2) = xa
xmat(1,3) = xa**2
xmat(1,4) = xa**3
xmat(2,1) = 1.0
xmat(2,2) = xb
xmat(2,3) = xb**2
xmat(2,4) = xb**3
xmat(3,1) = 0.0
xmat(3,2) = 1.0D0
xmat(3,3) = 2.0D0 * xa
xmat(3,4) = 3.0D0 * xa**2
xmat(4,1) = 0.0
xmat(4,2) = 1.0D0
xmat(4,3) = 2.0D0 * xb
xmat(4,4) = 3.0D0 * xb**2
    
ymat(1) = ya
ymat(2) = yb
ymat(3) = ypa
ymat(4) = ypb
    
!CALL sls_solve(xmat, ymat, c)
CALL minverse(xmat, xmat_inv, 4)
c = MATMUL(xmat_inv, ymat)
y = c(1) + c(2) * r + c(3) * r**2 + c(4) * r**3

END FUNCTION f_cubic_knot_ab





FUNCTION f_quintic_knot_ab(r, xa, ya, ypa, yppa, xb, yb, ypb, yppb) RESULT (y)
! CUBIC KNOT
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: xa, ya, ypa, yppa, xb, yb, ypb, yppb
REAL(kind=DoubleReal) :: y
!############################################################
REAL(kind=DoubleReal) :: xmat(1:6,1:6)
REAL(kind=DoubleReal) :: ymat(1:6)
REAL(kind=DoubleReal) :: c(1:6)
!############################################################
xmat(1,1) = 1.0
xmat(1,2) = xa
xmat(1,3) = xa**2
xmat(1,4) = xa**3
xmat(1,5) = xa**4
xmat(1,6) = xa**5
xmat(2,1) = 1.0
xmat(2,2) = xb
xmat(2,3) = xb**2
xmat(2,4) = xb**3
xmat(2,5) = xb**4
xmat(2,6) = xb**5

xmat(3,1) = 0.0D0
xmat(3,2) = 1.0D0
xmat(3,3) = 2.0D0 * xa
xmat(3,4) = 3.0D0 * xa**2
xmat(3,5) = 4.0D0 * xa**3
xmat(3,6) = 5.0D0 * xa**4
xmat(4,1) = 0.0D0
xmat(4,2) = 1.0D0
xmat(4,3) = 2.0D0 * xb
xmat(4,4) = 3.0D0 * xb**2
xmat(4,5) = 4.0D0 * xb**3
xmat(4,6) = 5.0D0 * xb**4

xmat(3,1) = 0.0D0
xmat(3,2) = 0.0D0
xmat(3,3) = 2.0D0
xmat(3,4) = 6.0D0 * xa
xmat(3,5) = 12.0D0 * xa**2
xmat(3,6) = 20.0D0 * xa**3

xmat(4,1) = 0.0D0
xmat(4,2) = 0.0D0
xmat(4,3) = 2.0D0
xmat(4,4) = 6.0D0 * xb
xmat(4,5) = 12.0D0 * xb**2
xmat(4,6) = 20.0D0 * xb**3
    
ymat(1) = ya
ymat(2) = yb
ymat(3) = ypa
ymat(4) = ypb
ymat(5) = yppa
ymat(6) = yppb
    
CALL sls_solve(xmat, ymat, c)
y = c(1) + c(2) * r + c(3) * r**2 + c(4) * r**3 + c(5) * r**4 + c(6) * r**5

END FUNCTION f_quintic_knot_ab



FUNCTION f_cubic_knot_spline(r, nodes) RESULT (y)
!############################################################
REAL(kind=DoubleReal) :: r
REAL(kind=DoubleReal) :: nodes(:,:)
REAL(kind=DoubleReal) :: y
!############################################################
INTEGER(kind=StandardInteger) :: n, j
!############################################################

! FIND n, n+1 TO SPLINE BETWEEN
IF(r .LT. nodes(1,1))THEN
  n = 1    
ELSE IF(r .GT. nodes(SIZE(nodes,1), 1))THEN
  n = SIZE(nodes,1)  - 1  
ELSE
  DO n = 1, SIZE(nodes,1) - 1
    IF(r .GE. nodes(n, 1) .AND. r .LT. nodes(n+1, 1))THEN
      EXIT
    END IF
  END DO
END IF  
y = f_cubic_knot_ab(r, nodes(n, 1), nodes(n, 2), nodes(n, 3), nodes(n+1, 1), nodes(n+1, 2), nodes(n+1, 3))
END FUNCTION f_cubic_knot_spline

