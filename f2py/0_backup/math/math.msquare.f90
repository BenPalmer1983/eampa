SUBROUTINE msquare(a, m)
!############################################################
IMPLICIT NONE
!############################################################
INTEGER(kind=StandardInteger), INTENT(IN) :: a(:)
REAL(kind=DoubleReal), INTENT(OUT) :: m(1:SIZE(a),1:SIZE(a))
!############################################################
INTEGER(kind=StandardInteger) :: i
!############################################################
m(:,:) = 0.0D0
DO i=1,SIZE(a,1)
  m(i,i) = a(i)
END DO
END SUBROUTINE msquare