FUNCTION str_compare(a, b) RESULT (scresult)
!############################################################
! String Compare
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: a
CHARACTER(LEN=*), INTENT(IN) :: b
LOGICAL :: scresult
!############################################################
INTEGER(kind=StandardInteger) :: n, la, lb
!############################################################

scresult = .FALSE.
la = 0
DO n = 1, len(a)
  !print *, ICHAR(a(n:n)), a(n:n)
  IF(.NOT. (a(n:n) .EQ. " " .OR. ICHAR(a(n:n)) .EQ. 0))THEN
    la = la + 1
  END IF
END DO
lb = 0

DO n = 1, len(b)
  !print *, ICHAR(b(n:n)), b(n:n)
  IF(.NOT. (b(n:n) .EQ. " " .OR. ICHAR(b(n:n)) .EQ. 0))THEN
    lb = lb + 1
  END IF
END DO
IF(la .NE. lb)THEN
  RETURN
END IF


scresult = .TRUE.
DO n = 1, min(len(a), len(b))
  IF(a(n:n) .NE. b(n:n))THEN
    scresult = .FALSE.
    RETURN
  END IF
END DO
RETURN
END FUNCTION str_compare





