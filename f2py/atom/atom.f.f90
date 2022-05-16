SUBROUTINE f(nlabels, r, rvec, clabels, forces)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::      nlabels(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              r(:)
REAL(kind=DoubleReal), INTENT(IN) ::              rvec(:, :)
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(INOUT) ::           forces(:, :)
!############################################################
INTEGER(KIND=StandardInteger) :: n, cc, nc, fn, fgroup
REAL(kind=DoubleReal), ALLOCATABLE :: density(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: dgab(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: dgba(:,:)
REAL(kind=DoubleReal), ALLOCATABLE :: epA(:)
REAL(kind=DoubleReal), ALLOCATABLE :: epB(:)
REAL(kind=DoubleReal) :: y(1:2)
REAL(kind=DoubleReal) :: fv(1:3)
!############################################################
! Neighbourlist Count
nc = SIZE(nlabels, 1)
cc = SIZE(clabels, 1)

! Allocate density array
ALLOCATE(density(1:cc, 0:group_max-1))
ALLOCATE(dgab(1:nc, 0:group_max-1))
ALLOCATE(dgba(1:nc, 0:group_max-1))
ALLOCATE(epA(0:group_max-1))
ALLOCATE(epB(0:group_max-1))

! Zero
density = 0.0D0
dgab = 0.D0
dgba = 0.D0
forces = 0.0D0

DO n = 1, nc
  ! PAIR ENERGY
  CALL getparr(1, nlabels(n,1), nlabels(n,2), 0, r(n), y)
  energy(1) = energy(1) + y(1)
  ! PAIR FORCE
  fv(:) = y(2) * rvec(n, 1:3)
  forces(nlabels(n, 3), 1:3) = forces(nlabels(n, 3), 1:3) + fv(1:3)
  forces(nlabels(n, 4), 1:3) = forces(nlabels(n, 4), 1:3) - fv(1:3)

  ! DENSITY
  IF(nlabels(n,1) .EQ. nlabels(n,2))THEN
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getparr(2, nlabels(n,1), -1, fgroup, r(n), y(:))
      density(nlabels(n,4), fgroup) = density(nlabels(n,4), fgroup) + y(1)
      density(nlabels(n,3), fgroup) = density(nlabels(n,3), fgroup) + y(1)
      dgab(n, fgroup) = y(2)
      dgba(n, fgroup) = y(2)
      fn = fn + 1
    END DO
  ELSE
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n,1), -1, fgroup, r(n), y(1))
      density(nlabels(n,4), fgroup) = density(nlabels(n,4), fgroup) + y(1)
      fn = fn + 1
    END DO
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,2), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n, 2), -1, fgroup, r(n), y(1))
      density(nlabels(n,3), fgroup) = density(nlabels(n,3), fgroup) + y(1)
      fn = fn + 1
    END DO
  END IF  
END DO

! Embedding Force
DO n = 1, nc
  epA = 0.0D0
  epB = 0.0D0
	fn = 0
	DO WHILE(embe_index_a(nlabels(n, 1), fn) .GT. -1)
		fgroup = embe_index_a(nlabels(n, 1), fn)
    CALL getpgrad(3, nlabels(n, 1), -1, fgroup, density(nlabels(n, 3), fgroup), epA(fgroup))
    CALL getpgrad(3, nlabels(n, 2), -1, fgroup, density(nlabels(n, 4), fgroup), epB(fgroup))
    !fv(:) = (dgab(n, fgroup) * epA(fgroup) + dgba(n, fgroup) * epB(fgroup)) * rvec(n, 1:3)
    !forces(nlabels(n, 3), 1:3) = forces(nlabels(n, 3), 1:3) + fv(1:3)
    !forces(nlabels(n, 4), 1:3) = forces(nlabels(n, 4), 1:3) - fv(1:3)
		fn = fn + 1
  END DO
  fv(1:3) = SUM(dgab(n, :) * epA(:) + dgba(n, :) * epB(:)) * rvec(n, 1:3)
  forces(nlabels(n, 3), 1:3) = forces(nlabels(n, 3), 1:3) + fv(1:3)
  forces(nlabels(n, 4), 1:3) = forces(nlabels(n, 4), 1:3) - fv(1:3)
END DO


DEALLOCATE(density)
DEALLOCATE(dgab)
DEALLOCATE(dgba)
DEALLOCATE(epA)
DEALLOCATE(epB)

END SUBROUTINE f