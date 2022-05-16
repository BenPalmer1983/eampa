SUBROUTINE verlet_step(clabels, cin, a0, uv, rcut, vinit, vdamp, steps, rebuild, dt, &
                       keep_history, cout, vout)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(IN) ::              cin(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              a0
REAL(kind=DoubleReal), INTENT(IN) ::              uv(1:3,1:3)
REAL(kind=DoubleReal), INTENT(IN) ::              rcut
REAL(kind=DoubleReal), INTENT(IN) ::              vinit(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              vdamp
INTEGER(kind=StandardInteger), INTENT(IN) ::      steps
INTEGER(kind=StandardInteger), INTENT(IN) ::      rebuild
REAL(kind=DoubleReal), INTENT(IN) ::              dt
LOGICAL, INTENT(IN) ::                            keep_history
REAL(kind=DoubleReal), INTENT(INOUT) ::           cout(:,:)
REAL(kind=DoubleReal), INTENT(INOUT) ::           vout(:,:)
!############################################################
INTEGER(kind=StandardInteger) ::                  n, k
INTEGER(kind=StandardInteger) ::                  nl_array_size
INTEGER(kind=StandardInteger) ::                  atom_count
INTEGER(kind=StandardInteger), ALLOCATABLE ::     nlabels(:,:)
REAL(kind=DoubleReal), ALLOCATABLE  ::            r(:)
REAL(kind=DoubleReal), ALLOCATABLE  ::            rvec(:, :)
INTEGER(kind=StandardInteger), ALLOCATABLE  ::    nlhalo(:)
INTEGER(kind=StandardInteger) ::                  nlsize
REAL(kind=DoubleReal) ::                          uv_inv(1:3,1:3)
REAL(kind=DoubleReal) ::                          c(1:SIZE(cin,1),1:SIZE(cin,2))
REAL(kind=DoubleReal) ::                          v(1:SIZE(cin,1),1:SIZE(cin,2))
REAL(kind=DoubleReal) ::                          vhalf(1:SIZE(cin,1),1:SIZE(cin,2))
REAL(kind=DoubleReal) ::                          forces(1:SIZE(cin,1),1:SIZE(cin,2))
REAL(kind=DoubleReal) ::                          dr(1:SIZE(cin,1),1:SIZE(cin,2))
REAL(kind=DoubleReal) ::                          energy(1:3)
!############################################################

! Atom count
atom_count = SIZE(cin,1)
nl_array_size = 100000

! Allocate arrays
ALLOCATE(nlabels(1:nl_array_size, 1:4))
ALLOCATE(r(1:nl_array_size))
ALLOCATE(rvec(1:nl_array_size, 1:3))
ALLOCATE(nlhalo(1:nl_array_size))

IF(keep_history)THEN
  IF(ALLOCATED(chistory)) DEALLOCATE(chistory)
	IF(ALLOCATED(ehistory)) DEALLOCATE(ehistory)

  ALLOCATE(chistory(1:steps+1, 1:atom_count, 1:3))
	ALLOCATE(ehistory(1:steps+1))
  chistory = 0.0D0
  chistory(1, :, :) = cin(:,:)
  ehistory = 0.0D0
END IF



! Set initial coords and velocity
c(:,:) = cin(:,:)
v(:,:) = vinit(:,:)

! Pre-compute coord transform
CALL minverse3(uv, uv_inv) 
uv_inv = (1.0D0 / a0) * uv_inv

DO n = 1, steps
  ! Compute neighbour list
  IF(MODULO((n-1), rebuild) .EQ. 0)THEN
    CALL nl(clabels, c(:,:), a0, uv, rcut, nlabels, r, rvec, nlhalo, nlsize)
  END IF

  ! Calculate force/acceleration
  CALL ef(nlabels(1:nlsize, 1:4), r(1:nlsize), rvec(1:nlsize, 1:3), &
                  clabels, energy, forces)
	IF(keep_history .AND. n .EQ. 1)THEN
    ehistory(n) = energy(3)
  END IF

  ! Compute half velocity
  vhalf(:, :) = v(:, :) + 0.5D0 * forces(:, :) * dt
  
  ! New position based on half velocity
  DO k = 1, atom_count
    dr(k, :) = MATMUL(uv_inv, vhalf(k, :) * dt)
  END DO
  c(:, :) = MODULO(c(:, :) + dr(:,:), 1.0D0)

  ! Calculate time step force/acceleration
  CALL nl_update_coords(c, a0, uv, nlabels(1:nlsize, 1:4), nlhalo(1:nlsize), &
                        r(1:nlsize), rvec(1:nlsize, 1:3))
	CALL ef(nlabels(1:nlsize, 1:4), r(1:nlsize), rvec(1:nlsize, 1:3), &
											  clabels, energy, forces)
  IF(keep_history)THEN
	  ehistory(n+1) = energy(3)
	END IF

  ! Update velocity
  v(:, :) = vhalf(:, :) + 0.5D0 * forces(:, :) * dt


  IF(vdamp .NE. 1.0D0)THEN
    v(:, :) = vdamp * v(:, :)
  END IF

  IF(keep_history)THEN
    chistory(n+1, :, :) = c(:, :)
  END IF
END DO

cout(:,:) = c(:,:)
vout(:,:) = v(:,:)


DEALLOCATE(nlabels)
DEALLOCATE(r)
DEALLOCATE(rvec)
DEALLOCATE(nlhalo)


END SUBROUTINE verlet_step

!REAL(kind=DoubleReal), INTENT(INOUT) ::           energy(1:3)
!REAL(kind=DoubleReal), INTENT(INOUT) ::           forces(:, :)


!SUBROUTINE ef(nlabels, r, rvec, clabels, energy, forces)

!nl_update_coords(cin, a0, uv, nl_labels, nlhalo, r, rvec)