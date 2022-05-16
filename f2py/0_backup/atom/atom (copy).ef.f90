SUBROUTINE ef(nlabels, r, rvec, clabels, energy, forces)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::      nlabels(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              r(:)
REAL(kind=DoubleReal), INTENT(IN) ::              rvec(:, :)
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(INOUT) ::           energy(1:3)
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
REAL(kind=DoubleReal) :: aaa, bbb, ccc, ddd, eee
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
energy = 0.0D0
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
  !forces(nlabels(n, 3), 1:3) = forces(nlabels(n, 3), 1:3) + fv(1:3)
  !forces(nlabels(n, 4), 1:3) = forces(nlabels(n, 4), 1:3) - fv(1:3)

  ! DENSITY
  IF(nlabels(n,1) .EQ. nlabels(n,2))THEN
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n,1), -1, fgroup, r(n), y(1))
      density(nlabels(n,4), fgroup) = density(nlabels(n,4), fgroup) + y(1)
      density(nlabels(n,3), fgroup) = density(nlabels(n,3), fgroup) + y(1)
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

! Embedding Energy
DO n = 1, cc
	fn = 0
	DO WHILE(embe_index_a(clabels(n), fn) .GT. -1)
		fgroup = embe_index_a(clabels(n), fn)
	  CALL getp(3, clabels(n), -1, fgroup, density(n, fgroup), y(1))
	  energy(2) = energy(2) + y(1)
		fn = fn + 1
  END DO
END DO


! Embedding Force
DO n = 1, nc
  !aa, bb, cc, dd, ee
  CALL getpgrad(2, nlabels(n, 1), -1, 0, r(n), aaa)
  CALL getpgrad(2, nlabels(n, 2), -1, 0, r(n), bbb)
  dgab(n, 0) = aaa
  dgba(n, 0) = bbb
!  CALL getp(3, nlabels(n, 1), -1, 0, aaa, ddd)
!  CALL getp(3, nlabels(n, 2), -1, 0, bbb, ddd)

!  fv(:) = (aaa*ddd + bbb*ccc) * rvec(n, 1:3)
!  forces(nlabels(n, 3), 1:3) = forces(nlabels(n, 3), 1:3) + fv(1:3)
!  forces(nlabels(n, 4), 1:3) = forces(nlabels(n, 4), 1:3) - fv(1:3)
END DO
DO n = 1, nc
  CALL getpgrad(3, nlabels(n, 1), -1, 0, density(nlabels(n, 3), 0), ccc)
  CALL getpgrad(3, nlabels(n, 2), -1, 0, density(nlabels(n, 4), 0), ddd)
  fv(:) = (dgab(n, 0) * ccc + dgba(n, 0) * ddd) * rvec(n, 1:3)
  forces(nlabels(n, 3), 1:3) = forces(nlabels(n, 3), 1:3) + fv(1:3)
  forces(nlabels(n, 4), 1:3) = forces(nlabels(n, 4), 1:3) - fv(1:3)
END DO




energy(3) = SUM(energy(1:2))


DEALLOCATE(density)
DEALLOCATE(dgab)
DEALLOCATE(dgba)
DEALLOCATE(epA)
DEALLOCATE(epB)

END SUBROUTINE ef