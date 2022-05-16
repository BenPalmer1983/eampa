

SUBROUTINE add_tabulated(ftype, a, b, fgroup, tab, fortran_pkey)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::     ftype
INTEGER(kind=StandardInteger), INTENT(IN) ::     a
INTEGER(kind=StandardInteger), INTENT(IN) ::     b
INTEGER(kind=StandardInteger), INTENT(IN) ::     fgroup
REAL(kind=DoubleReal), INTENT(IN) ::             tab(1:1001, 1:3)
INTEGER(kind=StandardInteger), INTENT(OUT) ::    fortran_pkey
!############################################################
INTEGER(kind=StandardInteger) :: key, n
!############################################################
fcount = fcount + 1
pot_fname(fcount) = blank
pot_fname(fcount) = "tab"
pot_is_tab(fcount) = 1
pot_ftype(fcount) = ftype
pot_a(fcount) = a
pot_b(fcount) = b
pot_fgroup(fcount) = fgroup
pot_tab(fcount, 1:1001, 1:3) = tab(1:1001, 1:3)
group_max = MAX(group_max, fgroup + 1)
IF(ftype .EQ. 1)THEN
  DO n = 0, lgmax 
    IF(pair_index_a(MIN(a, b), n) .EQ. -1)THEN
      pair_index_a(MIN(a, b), n) = MAX(a, b)
      EXIT
    END IF 
  END DO
  key = pair_key(a, b)
  pair_index_b(key) = fcount
ELSE IF(ftype .EQ. 2)THEN
  ! Independent of embedded atom
  IF(b .EQ. -1)THEN
    DO n = 0, lgmax 
      IF(dens_index_a(a, n) .EQ. -1)THEN
        dens_index_a(a, n) = fgroup
        EXIT
      END IF
    END DO
    dens_index_b(a, fgroup) = fcount
  ! Dependent on embedded atom
  ELSE
    
  END IF
ELSE IF(ftype .EQ. 3)THEN
  DO n = 0, lgmax 
    IF(embe_index_a(a, n) .EQ. -1)THEN
        embe_index_a(a, n) = fgroup
      EXIT
    END IF
  END DO 
  embe_index_b(a, fgroup) = fcount
END IF
fortran_pkey = fcount
!############################################################
END SUBROUTINE add_tabulated




SUBROUTINE add_analytic(ftype, a, b, fgroup, fname, p, pf, fortran_pkey)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) ::     ftype
INTEGER(kind=StandardInteger), INTENT(IN) ::     a
INTEGER(kind=StandardInteger), INTENT(IN) ::     b
INTEGER(kind=StandardInteger), INTENT(IN) ::     fgroup
CHARACTER*32 ::                                  fname
REAL(kind=DoubleReal), INTENT(IN) ::             p(:)
REAL(kind=DoubleReal), INTENT(IN) ::             pf(:)
INTEGER(kind=StandardInteger), INTENT(OUT) ::    fortran_pkey
!############################################################
INTEGER(kind=StandardInteger) :: key, n
!############################################################
fcount = fcount + 1
pot_fname(fcount) = fname
pot_is_tab(fcount) = 0
pot_ftype(fcount) = ftype
pot_a(fcount) = a
pot_b(fcount) = b
pot_fgroup(fcount) = fgroup
pot_p(fcount, 1:SIZE(p)) = p(:)
pot_pf(fcount, 1:SIZE(pf)) = pf(:)
pot_psize(fcount) = SIZE(p)
pot_pfsize(fcount) = SIZE(pf)
group_max = MAX(group_max, fgroup)
IF(ftype .EQ. 1)THEN
  DO n = 0, lgmax 
    IF(pair_index_a(MIN(a, b), n) .EQ. -1)THEN
      pair_index_a(MIN(a, b), n) = MAX(a, b)
      EXIT
    END IF 
  END DO
  key = pair_key(a, b)
  pair_index_b(key) = fcount
ELSE IF(ftype .EQ. 2)THEN
  ! Independent of embedded atom
  IF(b .EQ. -1)THEN
    DO n = 0, lgmax 
      IF(dens_index_a(a, n) .EQ. -1)THEN
        dens_index_a(a, n) = fgroup
        EXIT
      END IF
    END DO
    dens_index_b(a, fgroup) = fcount
  ! Dependent on embedded atom
  ELSE
    
  END IF
ELSE IF(ftype .EQ. 3)THEN
  DO n = 0, lgmax 
    IF(embe_index_a(a, n) .EQ. -1)THEN
        embe_index_a(a, n) = fgroup
      EXIT
    END IF
  END DO 
  embe_index_b(a, fgroup) = fcount
END IF
fortran_pkey = fcount
!############################################################
END SUBROUTINE add_analytic




SUBROUTINE update_tabulated(tab, fortran_pkey)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) ::            tab(1:1001, 1:3)
INTEGER(kind=StandardInteger), INTENT(IN) ::    fortran_pkey
!############################################################
!############################################################
pot_tab(fortran_pkey, 1:1001, 1:3) = tab(1:1001, 1:3)
!############################################################
END SUBROUTINE update_tabulated



SUBROUTINE update_analytic(fname, p, pf, fortran_pkey)
!############################################################
IMPLICIT NONE
!############################################################
CHARACTER*32 ::                                  fname
REAL(kind=DoubleReal), INTENT(IN) ::             p(:)
REAL(kind=DoubleReal), INTENT(IN) ::             pf(:)
INTEGER(kind=StandardInteger), INTENT(IN) ::    fortran_pkey
!############################################################
!############################################################
pot_fname(fortran_pkey) = fname
pot_p(fortran_pkey, 1:SIZE(p)) = p(:)
pot_pf(fortran_pkey, 1:SIZE(pf)) = pf(:)
pot_psize(fortran_pkey) = SIZE(p)
pot_pfsize(fortran_pkey) = SIZE(pf)
!############################################################
END SUBROUTINE update_analytic




SUBROUTINE getp(ftype, a, b, fgroup, x, y)
!############################################################
IMPLICIT NONE
INTEGER(kind=StandardInteger), INTENT(IN) ::     ftype
INTEGER(kind=StandardInteger), INTENT(IN) ::     a
INTEGER(kind=StandardInteger), INTENT(IN) ::     b
INTEGER(kind=StandardInteger), INTENT(IN) ::     fgroup
REAL(kind=DoubleReal), INTENT(IN) ::             x
REAL(kind=DoubleReal), INTENT(OUT) ::            y
!############################################################
INTEGER(kind=StandardInteger) ::                 pkey, key
!############################################################
y = 0.0D0

IF(ftype .EQ. 1)THEN
  key = pair_key(a, b)
  pkey = pair_index_b(key)
ELSE IF(ftype .EQ. 2)THEN
  
  pkey = dens_index_b(a, fgroup)

ELSE IF(ftype .EQ. 3)THEN
  pkey = embe_index_b(a, fgroup) 
END IF

IF(pot_is_tab(pkey) .EQ. 1)THEN
  IF(.NOT. (x .LT. pot_tab(pkey, 1, 1) .OR. x .GT. pot_tab(pkey, tabsize, 1)))THEN
    CALL interpolate(x, pot_tab(pkey, :, 1), pot_tab(pkey, :, 2), 4, 0, y)
  END IF
ELSE
  CALL f(pot_fname(pkey), x, pot_p(pkey, 1:pot_psize(pkey)), pot_pf(pkey, 1:pot_pfsize(pkey)), y)
END IF
!############################################################
END SUBROUTINE getp


SUBROUTINE getpgrad(ftype, a, b, fgroup, x, dydx)
!############################################################
IMPLICIT NONE
INTEGER(kind=StandardInteger), INTENT(IN) ::     ftype
INTEGER(kind=StandardInteger), INTENT(IN) ::     a
INTEGER(kind=StandardInteger), INTENT(IN) ::     b
INTEGER(kind=StandardInteger), INTENT(IN) ::     fgroup
REAL(kind=DoubleReal), INTENT(IN) ::             x
REAL(kind=DoubleReal), INTENT(OUT) ::            dydx
!############################################################
INTEGER(kind=StandardInteger) ::                 pkey, key
!############################################################
dydx = 0.0D0

IF(ftype .EQ. 1)THEN
  key = pair_key(a, b)
  pkey = pair_index_b(key)
ELSE IF(ftype .EQ. 2)THEN
  
  pkey = dens_index_b(a, fgroup)

ELSE IF(ftype .EQ. 3)THEN
  pkey = embe_index_b(a, fgroup) 
END IF

IF(pot_is_tab(pkey) .EQ. 1)THEN
  IF(.NOT. (x .LT. pot_tab(pkey, 1, 1) .OR. x .GT. pot_tab(pkey, tabsize, 1)))THEN
    CALL interpolate(x, pot_tab(pkey, :, 1), pot_tab(pkey, :, 3), 4, 0, dydx)
  END IF
ELSE
  CALL fgrad(pot_fname(pkey), x, pot_p(pkey, 1:pot_psize(pkey)), pot_pf(pkey, 1:pot_pfsize(pkey)), dydx)
END IF
!############################################################
END SUBROUTINE getpgrad


SUBROUTINE getparr(ftype, a, b, fgroup, x, y)
!############################################################
IMPLICIT NONE
INTEGER(kind=StandardInteger), INTENT(IN) ::     ftype
INTEGER(kind=StandardInteger), INTENT(IN) ::     a
INTEGER(kind=StandardInteger), INTENT(IN) ::     b
INTEGER(kind=StandardInteger), INTENT(IN) ::     fgroup
REAL(kind=DoubleReal), INTENT(IN) ::             x
REAL(kind=DoubleReal), INTENT(OUT) ::            y(1:2)
!############################################################
INTEGER(kind=StandardInteger) ::                 pkey, key
!############################################################
y = 0.0D0
IF(ftype .EQ. 1)THEN
  key = pair_key(a, b)
  pkey = pair_index_b(key)
ELSE IF(ftype .EQ. 2)THEN
  
  pkey = dens_index_b(a, fgroup)

ELSE IF(ftype .EQ. 3)THEN
  pkey = embe_index_b(a, fgroup) 
END IF

IF(pot_is_tab(pkey) .EQ. 1)THEN
  IF(.NOT. (x .LT. pot_tab(pkey, 1, 1) .OR. x .GT. pot_tab(pkey, tabsize, 1)))THEN
    CALL interpolate(x, pot_tab(pkey, :, 1), pot_tab(pkey, :, 2), 4, 0, y(1))
    CALL interpolate(x, pot_tab(pkey, :, 1), pot_tab(pkey, :, 3), 4, 0, y(2))
  END IF
ELSE
  CALL f(pot_fname(pkey), x, pot_p(pkey, 1:pot_psize(pkey)), pot_pf(pkey, 1:pot_pfsize(pkey)), y(1))
  CALL fgrad(pot_fname(pkey), x, pot_p(pkey, 1:pot_psize(pkey)), pot_pf(pkey, 1:pot_pfsize(pkey)), y(2))
END IF
!############################################################
END SUBROUTINE getparr



! KEYS
FUNCTION pair_key(a, b) RESULT (key)
! Calculates unique key for two keys (order of keys NOT important)
! (A,B) = (B,A)
!###########################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a, b
INTEGER(kind=StandardInteger) :: key
!###########################################################
key = (MAX(a+1, b+1) * (MAX(a+1, b+1) - 1)) / 2 + MIN(a+1, b+1)-1
END FUNCTION pair_key
