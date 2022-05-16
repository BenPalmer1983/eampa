SUBROUTINE e(nlabels, r, clabels, energy)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::      nlabels(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              r(:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(INOUT) ::           energy(1:3)
!############################################################
INTEGER(KIND=StandardInteger) ::                  n, cc, nc, fn, fgroup
REAL(kind=DoubleReal), ALLOCATABLE ::             density(:,:)
REAL(kind=DoubleReal), ALLOCATABLE ::             density_thread(:,:)
REAL(kind=DoubleReal) ::                          energy_thread(1:3)
REAL(kind=DoubleReal) ::                          y
!############################################################
! Neighbourlist Count
nc = SIZE(nlabels, 1)
cc = SIZE(clabels, 1)

! Allocate density array
ALLOCATE(density(1:cc, 0:group_max-1))
ALLOCATE(density_thread(1:cc, 0:group_max-1))

! Zero
energy = 0.0D0
density = 0.0D0

!$OMP PARALLEL &
!$OMP PRIVATE(n, y, fgroup, fn, energy_thread, &
!$OMP density_thread) &
!$OMP SHARED(nc, nlabels, r, energy, density)
energy_thread(1) = 0.0D0
density_thread(:,:) = 0.0D0
!$OMP DO
DO n = 1, nc
  ! PAIR ENERGY
  CALL getp(1, nlabels(n,1), nlabels(n,2), 0, r(n), y)
  energy_thread(1) = energy_thread(1) + y

  ! DENSITY
  IF(nlabels(n,1) .EQ. nlabels(n,2))THEN
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n,1), -1, fgroup, r(n), y)
      density_thread(nlabels(n,4), fgroup) = density_thread(nlabels(n,4), fgroup) + y
      density_thread(nlabels(n,3), fgroup) = density_thread(nlabels(n,3), fgroup) + y
      fn = fn + 1
    END DO
  ELSE
    fn = 0
    DO WHILE(dens_index_a(nlabels(n,1), fn) .GT. -1)
      fgroup = dens_index_a(nlabels(n,1), fn)
      CALL getp(2, nlabels(n,1), -1, fgroup, r(n), y)
      density_thread(nlabels(n,4), fgroup) = density_thread(nlabels(n,4), fgroup) + y
      CALL getp(2, nlabels(n, 2), -1, fgroup, r(n), y)
      density_thread(nlabels(n,3), fgroup) = density_thread(nlabels(n,3), fgroup) + y
      fn = fn + 1
    END DO
  END IF  
END DO
!$OMP END DO
!$OMP CRITICAL
energy(1) = energy(1) + energy_thread(1)
density(:,:) = density(:,:) + density_thread(:,:)
!$OMP END CRITICAL
!$OMP END PARALLEL

! Embedding Energy
! If a small number of atoms, just use one thread
IF(cc .LE. 100)THEN
  DO n = 1, cc
	  fn = 0
	  DO WHILE(embe_index_a(clabels(n), fn) .GT. -1)
	  	fgroup = embe_index_a(clabels(n), fn)
	    CALL getp(3, clabels(n), -1, fgroup, density(n, fgroup), y)
	    energy(2) = energy(2) + y
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
	    CALL getp(3, clabels(n), -1, fgroup, density(n, fgroup), y)
	    energy_thread(2) = energy_thread(2) + y
		  fn = fn + 1
    END DO
  END DO
!$OMP END DO
!$OMP CRITICAL
  energy(2) = energy(2) + energy_thread(2)
!$OMP END CRITICAL
!$OMP END PARALLEL
END IF

! Energy (pair + embed)
energy(3) = SUM(energy(1:2))


DEALLOCATE(density)
DEALLOCATE(density_thread)


END SUBROUTINE e



