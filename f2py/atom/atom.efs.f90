SUBROUTINE efs(nlabels, r, rvec, inhalo, clabels, volume, energy, forces, stress)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::      nlabels(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              r(:)
REAL(kind=DoubleReal), INTENT(IN) ::              rvec(:, :)
INTEGER(KIND=StandardInteger), INTENT(IN) ::      inhalo(:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(IN) ::              volume
REAL(kind=DoubleReal), INTENT(INOUT) ::           energy(1:3)
REAL(kind=DoubleReal), INTENT(INOUT) ::           forces(:, :)
REAL(kind=DoubleReal), INTENT(INOUT) ::           stress(:, :)
!############################################################
INTEGER(KIND=StandardInteger) :: n, cc, nc, fn, fgroup
INTEGER(KIND=StandardInteger) :: i, j
REAL(kind=DoubleReal), ALLOCATABLE ::             density(:,:)
REAL(kind=DoubleReal), ALLOCATABLE ::             nl_force(:,:)
REAL(kind=DoubleReal), ALLOCATABLE ::             fgradA(:)
REAL(kind=DoubleReal), ALLOCATABLE ::             fgradB(:)
REAL(kind=DoubleReal), ALLOCATABLE ::             dgradA(:)
REAL(kind=DoubleReal), ALLOCATABLE ::             dgradB(:)
REAL(kind=DoubleReal), ALLOCATABLE ::             density_thread(:,:)
REAL(kind=DoubleReal), ALLOCATABLE ::             nl_force_thread(:,:)
REAL(kind=DoubleReal) ::                          energy_thread(1:3)
REAL(kind=DoubleReal) ::                          forces_thread(1:SIZE(forces,1), 1:SIZE(forces,2))
REAL(kind=DoubleReal) ::                          y(1:2)
REAL(kind=DoubleReal) ::                          fv(1:3)
!############################################################
! Neighbourlist Count
nc = SIZE(nlabels, 1)
cc = SIZE(clabels, 1)

! Allocate density array
ALLOCATE(density(1:cc, 0:group_max-1))
ALLOCATE(density_thread(1:cc, 0:group_max-1))
ALLOCATE(fgradA(0:group_max-1))
ALLOCATE(fgradB(0:group_max-1))
ALLOCATE(dgradA(0:group_max-1))
ALLOCATE(dgradB(0:group_max-1))
ALLOCATE(nl_force(1:nc, 1:3))
ALLOCATE(nl_force_thread(1:nc, 1:3))


! Zero
energy = 0.0D0
density = 0.0D0
forces = 0.0D0
stress = 0.0D0
nl_force = 0.0D0

!$OMP PARALLEL &
!$OMP PRIVATE(n, y, fgroup, fn, energy_thread, density_thread, &
!$OMP forces_thread, fv, nl_force_thread) &
!$OMP SHARED(nc, nlabels, r, rvec, energy, density, forces, nl_force)
energy_thread(1) = 0.0D0
density_thread(:,:) = 0.0D0
forces_thread(:,:) = 0.0D0
nl_force_thread(:,:) = 0.0D0
!$OMP DO
DO n = 1, nc
  ! PAIR ENERGY
  CALL getparr(1, nlabels(n,1), nlabels(n,2), 0, r(n), y)
  energy_thread(1) = energy_thread(1) + y(1)
  ! PAIR FORCE
  fv(:) = y(2) * rvec(n, 1:3)
  forces_thread(nlabels(n, 3), 1:3) = forces_thread(nlabels(n, 3), 1:3) + fv(1:3)
  forces_thread(nlabels(n, 4), 1:3) = forces_thread(nlabels(n, 4), 1:3) - fv(1:3)
  nl_force_thread(n, 1:3) = -fv(1:3)

  ! DENSITY
  IF(nlabels(n,1) .EQ. nlabels(n,2))THEN
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n,1), -1, fgroup, r(n), y(1))
      density_thread(nlabels(n,4), fgroup) = density_thread(nlabels(n,4), fgroup) + y(1)
      density_thread(nlabels(n,3), fgroup) = density_thread(nlabels(n,3), fgroup) + y(1)
      fn = fn + 1
    END DO
  ELSE
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n,1), -1, fgroup, r(n), y(1))
      density_thread(nlabels(n,4), fgroup) = density_thread(nlabels(n,4), fgroup) + y(1)
      fn = fn + 1
    END DO
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,2), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,2), fn)
      CALL getp(2, nlabels(n, 2), -1, fgroup, r(n), y(1))
      density_thread(nlabels(n,3), fgroup) = density_thread(nlabels(n,3), fgroup) + y(1)
      fn = fn + 1
    END DO
  END IF  
END DO
!$OMP END DO
!$OMP CRITICAL
energy(1) = energy(1) + energy_thread(1)
density(:,:) = density(:,:) + density_thread(:,:)
forces(:,:) = forces(:,:) + forces_thread(:,:)
nl_force(:,:) = nl_force(:,:) + nl_force_thread(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

! Embedding Energy
! If a small number of atoms, just use one thread
IF(cc .LE. 100)THEN
  DO n = 1, cc
	  fn = 0
	  DO WHILE(embe_index_a(clabels(n), fn) .GT. -1)
	  	fgroup = embe_index_a(clabels(n), fn)
	    CALL getp(3, clabels(n), -1, fgroup, density(n, fgroup), y(1))
	    energy(2) = energy(2) + y(1)
	  	fn = fn + 1
    END DO
  END DO
ELSE
!$OMP PARALLEL &
!$OMP PRIVATE(n, y, fgroup, fn, energy_thread) &
!$OMP SHARED(cc, nlabels, clabels, energy, density)
  energy_thread(2) = 0.0D0
!$OMP DO
  DO n = 1, cc
	  fn = 0
  	DO WHILE(embe_index_a(clabels(n), fn) .GT. -1)
	  	fgroup = embe_index_a(clabels(n), fn)
	    CALL getp(3, clabels(n), -1, fgroup, density(n, fgroup), y(1))
	    energy_thread(2) = energy_thread(2) + y(1)
		  fn = fn + 1
    END DO
  END DO
!$OMP END DO
!$OMP CRITICAL
  energy(2) = energy(2) + energy_thread(2)
!$OMP END CRITICAL
!$OMP END PARALLEL
END IF

! Embedding Force
!$OMP PARALLEL &
!$OMP PRIVATE(n, fgroup, fn, forces_thread, dgradA, dgradB, fgradA, fgradB, fv) &
!$OMP SHARED(nc, nlabels, rvec, density, forces)
forces_thread(:,:) = 0.0D0
nl_force_thread(:,:) = 0.0D0
!$OMP DO
DO n = 1, nc
  dgradA(:) = 0.0D0
  dgradB(:) = 0.0D0
  fgradA(:) = 0.0D0
  fgradB(:) = 0.0D0

  ! Loop through groups and compute density and embedding gradients
	fn = 0
	DO WHILE(embe_index_a(nlabels(n, 1), fn) .GT. -1)
		fgroup = embe_index_a(nlabels(n, 1), fn)
    CALL getpgrad(2, nlabels(n,1), -1, fgroup, r(n), dgradA(fgroup))
    CALL getpgrad(2, nlabels(n,2), -1, fgroup, r(n), dgradB(fgroup))
    CALL getpgrad(3, nlabels(n,1), -1, fgroup, density(nlabels(n, 3), fgroup), fgradA(fgroup))
    CALL getpgrad(3, nlabels(n,2), -1, fgroup, density(nlabels(n, 4), fgroup), fgradB(fgroup))
		fn = fn + 1
  END DO

  ! Compute forces
  fv(1:3) = SUM(fgradA(:) * dgradB(:) + fgradB(:) * dgradA(:)) * rvec(n, 1:3)
  forces_thread(nlabels(n, 3), 1:3) = forces_thread(nlabels(n, 3), 1:3) + fv(1:3)
  forces_thread(nlabels(n, 4), 1:3) = forces_thread(nlabels(n, 4), 1:3) - fv(1:3)
  nl_force_thread(n, 1:3) = -fv(1:3)
END DO
!$OMP END DO
!$OMP CRITICAL
forces(:,:) = forces(:,:) + forces_thread(:,:)
nl_force(:,:) = nl_force(:,:) + nl_force_thread(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

! Energy (pair + embed)
energy(3) = SUM(energy(1:2))

! Stress
DO n = 1, nc
  IF(inhalo(n) .NE. 0)THEN  
    DO i = 1,3
      DO j = 1,3
        stress(i,j) = stress(i,j) + (rvec(n, i) * nl_force(n,j))
      END DO
    END DO
  END IF
END DO
stress(1:3, 1:3) = stress(1:3, 1:3) / (2.0D0 * volume)



DEALLOCATE(density)
DEALLOCATE(density_thread)
DEALLOCATE(fgradA)
DEALLOCATE(fgradB)
DEALLOCATE(dgradA)
DEALLOCATE(dgradB)

END SUBROUTINE efs