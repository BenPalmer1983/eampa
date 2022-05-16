SUBROUTINE nl_change_cell(a0_new, uv_new, a0_old, uv_old, r, rvec)
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a0_new
REAL(kind=DoubleReal), INTENT(IN) :: uv_new(1:3,1:3)
REAL(kind=DoubleReal), INTENT(IN) :: a0_old
REAL(kind=DoubleReal), INTENT(IN) :: uv_old(1:3,1:3)
REAL(kind=DoubleReal), INTENT(INOUT) :: r(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: rvec(:, :)
!############################################################
REAL(kind=DoubleReal) :: uv_old_inv(1:3,1:3)
REAL(kind=DoubleReal) :: m(1:3,1:3)
REAL(kind=DoubleReal) :: amult, rnew, ra(1:3)
INTEGER(kind=StandardInteger) :: n
!############################################################
CALL minverse3(uv_old, uv_old_inv)
amult = a0_new / a0_old
m = amult * MATMUL(uv_new, uv_old_inv)

DO n=1, SIZE(r)
  ra(1) = m(1,1) * rvec(n, 1) + m(1,2) * rvec(n, 2) + m(1,3) * rvec(n, 3)
  ra(2) = m(2,1) * rvec(n, 1) + m(2,2) * rvec(n, 2) + m(2,3) * rvec(n, 3)
  ra(3) = m(3,1) * rvec(n, 1) + m(3,2) * rvec(n, 2) + m(3,3) * rvec(n, 3)
  ra(1:3) = r(n) * ra(1:3)  
  r(n) = SQRT(SUM(ra(1:3)**2))
  rvec(n, 1:3) = ra(1:3) / r(n)
END DO

!############################################################
END SUBROUTINE nl_change_cell
