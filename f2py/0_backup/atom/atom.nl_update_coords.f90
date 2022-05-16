SUBROUTINE nl_update_coords(cin, a0, uv, nl_labels, nlhalo, r, rvec)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: cin(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: a0
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3,1:3)
INTEGER(kind=StandardInteger), INTENT(IN) :: nl_labels(:,:)
INTEGER(kind=StandardInteger), INTENT(IN) :: nlhalo(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: r(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: rvec(:, :)
!############################################################
INTEGER(kind=StandardInteger) :: n, i, j, k
REAL(kind=DoubleReal) :: ijk(-13:13, 1:3)
REAL(kind=DoubleReal) :: a(1:3)
REAL(kind=DoubleReal) :: b(1:3)
INTEGER(kind=StandardInteger) :: nl_size
!############################################################
DO k=-1,1
  DO j=-1,1
    DO i=-1,1
      n = i+3*j+9*k
      ijk(n,1) = i
      ijk(n,2) = j
      ijk(n,3) = k
    END DO
  END DO
END DO



nl_size = SIZE(r, 1)

DO n =1, nl_size
  a(:) = cin(nl_labels(n, 3),:)
  b(:) = cin(nl_labels(n, 4),:) + ijk(nlhalo(n),:)
  rvec(n, :) = MATMUL(a0 * uv(1:3,1:3), b(:) - a(:))
	r(n) = SQRT(SUM(rvec(n, :)**2))
  rvec(n, :) = rvec(n, :) / r(n)
END DO




!############################################################
END SUBROUTINE nl_update_coords