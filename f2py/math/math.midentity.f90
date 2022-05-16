SUBROUTINE midentity(n, m)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: n
REAL(kind=DoubleReal), INTENT(OUT) :: m(1:n,1:n)
!############################################################
INTEGER(kind=StandardInteger) :: i
!############################################################
m(:,:) = 0.0D0
DO i=1,n
  m(i,i) = 1.0D0
END DO
END SUBROUTINE midentity