SUBROUTINE relax(clabels, cin, a0, uv, rcut, steps, rebuild, dt, cout)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(KIND=StandardInteger), INTENT(IN) ::      clabels(:)
REAL(kind=DoubleReal), INTENT(IN) ::              cin(:,:)
REAL(kind=DoubleReal), INTENT(IN) ::              a0
REAL(kind=DoubleReal), INTENT(IN) ::              uv(1:3,1:3)
REAL(kind=DoubleReal), INTENT(IN) ::              rcut
INTEGER(kind=StandardInteger), INTENT(IN) ::      steps
INTEGER(kind=StandardInteger), INTENT(IN) ::      rebuild
REAL(kind=DoubleReal), INTENT(IN) ::              dt
REAL(kind=DoubleReal), INTENT(INOUT) ::           cout(:,:)
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
REAL(kind=DoubleReal) ::                          energy_a = 0.0D0
!############################################################

! Atom count
atom_count = SIZE(cin,1)
nl_array_size = 100000

! Allocate arrays
ALLOCATE(nlabels(1:nl_array_size, 1:4))
ALLOCATE(r(1:nl_array_size))
ALLOCATE(rvec(1:nl_array_size, 1:3))
ALLOCATE(nlhalo(1:nl_array_size))


! Set initial coords and velocity
c(:,:) = cin(:,:)
v(:,:) = 0.0D0

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
  IF(n .GT. 1 .AND. energy(3) .GT. energy_a) THEN
    v(:,:) = 0.0D0
  END IF
  energy_a = energy(3)

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
  
  ! Update velocity
  v(:, :) = vhalf(:, :) + 0.5D0 * forces(:, :) * dt




END DO

cout(:,:) = c(:,:)


DEALLOCATE(nlabels)
DEALLOCATE(r)
DEALLOCATE(rvec)
DEALLOCATE(nlhalo)

END SUBROUTINE relax

!REAL(kind=DoubleReal), INTENT(INOUT) ::           energy(1:3)
!REAL(kind=DoubleReal), INTENT(INOUT) ::           forces(:, :)


!SUBROUTINE ef(nlabels, r, rvec, clabels, energy, forces)

!nl_update_coords(cin, a0, uv, nl_labels, nlhalo, r, rvec)