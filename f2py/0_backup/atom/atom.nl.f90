SUBROUTINE nl(clabels, cin, a0, uv, rcut, nl_labels, r, rvec, nlhalo, nlsize)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: clabels(:)
REAL(kind=DoubleReal), INTENT(IN) :: cin(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: a0
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3,1:3)
REAL(kind=DoubleReal), INTENT(IN) :: rcut
INTEGER(kind=StandardInteger), INTENT(INOUT) :: nl_labels(:,:)
REAL(kind=DoubleReal), INTENT(INOUT) :: r(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: rvec(:, :)
INTEGER(kind=StandardInteger), INTENT(INOUT) :: nlhalo(:)
INTEGER(kind=StandardInteger), INTENT(OUT) :: nlsize
!############################################################
REAL(kind=DoubleReal) :: ccoords(1:SIZE(cin,1),1:3)
REAL(kind=DoubleReal) :: gcoords(1:27*SIZE(cin,1),1:3)
INTEGER(kind=StandardInteger) :: glabels(1:27*SIZE(cin,1))
INTEGER(kind=StandardInteger) :: gids(1:27*SIZE(cin,1))
INTEGER(kind=StandardInteger) :: halo(1:27*SIZE(cin,1))
REAL(kind=DoubleReal) :: rcutsq
REAL(kind=DoubleReal) :: rsq
REAL(kind=DoubleReal) :: shift(1:3)
REAL(kind=DoubleReal) :: rvect(1:3)
INTEGER(kind=StandardInteger) :: m, n, i, j, k
!############################################################
DO n=1,SIZE(cin, 1)
  ccoords(n, :) = a0 * MATMUL(uv, cin(n, :))
END DO

m = 0
DO i=-1,1
  DO j=-1,1 
    DO k=-1,1
      DO n=1, SIZE(cin, 1)
        m = m + 1
        shift(1) = i
        shift(2) = j
        shift(3) = k
        gcoords(m, :) = a0 * MATMUL(uv, cin(n, :) + shift(:))
        glabels(m) = clabels(n)
        gids(m) = n
        halo(m) = i + 3 * j + 9 * k
      END DO
    END DO
  END DO
END DO

rcutsq = rcut * rcut

k = 0
DO n=1, SIZE(ccoords, 1)
  DO m=1, SIZE(gcoords, 1)
    IF(gids(m) .GT. n)THEN
      rvect(:) = gcoords(m, :) - ccoords(n, :)
      rsq = SUM(rvect(:)**2)
      IF(rsq .LE. rcutsq)THEN
        k = k + 1
        nl_labels(k, 1) = clabels(n)
        nl_labels(k, 2) = glabels(m)
        nl_labels(k, 3) = n
        nl_labels(k, 4) = gids(m)
        r(k) = sqrt(rsq)
        rvec(k, :) = rvect(:) / r(k)
        nlhalo(k) = halo(m)
      END IF
    END IF
  END DO
END DO
nlsize = k
!CALL nl_ncount(nl_labels, SIZE(ccoords, 1))

!############################################################
END SUBROUTINE nl


SUBROUTINE nl_ncount(nlabels, csize)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: nlabels(:,:)
INTEGER(kind=StandardInteger), INTENT(IN) :: csize
!############################################################
INTEGER(kind=StandardInteger) :: n, ns
INTEGER(kind=StandardInteger) :: ncount(1:csize)
!############################################################
ns = SIZE(nlabels, 1)
ncount = 0
print *, csize
DO n =1, ns
  ncount(nlabels(n,3)) = ncount(nlabels(n,3)) + 1
  ncount(nlabels(n,4)) = ncount(nlabels(n,4)) + 1  
END DO

print *, ncount(:)
print *, SUM(ncount(:))





!############################################################
END SUBROUTINE nl_ncount



