


SUBROUTINE d(nlabels, r, clabels, density)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::      nlabels(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              r(:)
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(INOUT) ::           density(:, :)
!############################################################
INTEGER(KIND=StandardInteger) ::                  n, cc, nc, fn, fgroup
REAL(kind=DoubleReal), ALLOCATABLE ::             density_thread(:,:)
REAL(kind=DoubleReal) ::                          y
!############################################################
! Neighbourlist Count
nc = SIZE(nlabels, 1)
cc = SIZE(clabels, 1)

! Allocate density array
ALLOCATE(density_thread(1:cc, 0:group_max-1))

! Zero
density = 0.0D0

!$OMP PARALLEL &
!$OMP PRIVATE(n, y, fgroup, fn, &
!$OMP density_thread) &
!$OMP SHARED(nc, nlabels, r, density)
density_thread(:,:) = 0.0D0
!$OMP DO
DO n = 1, nc
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
density(1:cc,1:group_max) = density(1:cc,1:group_max) + density_thread(1:cc,0:group_max-1)
!$OMP END CRITICAL
!$OMP END PARALLEL

DEALLOCATE(density_thread)


END SUBROUTINE d



