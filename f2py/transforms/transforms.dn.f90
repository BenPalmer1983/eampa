FUNCTION dn(uv, s, n) RESULT (uvout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal), INTENT(IN) :: s
REAL(kind=DoubleReal) :: uvout(1:3, 1:3)
!############################################################
REAL(kind=DoubleReal) :: d(1:3, 1:3)
INTEGER(kind=StandardInteger) :: n
!############################################################
IF(n .EQ. 1)THEN
  uvout(1:3,1:3) = d1(uv(1:3,1:3), s)
ELSE IF(n .EQ. 2)THEN
  uvout(1:3,1:3) = d2(uv(1:3,1:3), s)
ELSE IF(n .EQ. 3)THEN
  uvout(1:3,1:3) = d3(uv(1:3,1:3), s)
ELSE IF(n .EQ. 4)THEN
  uvout(1:3,1:3) = d4(uv(1:3,1:3), s)
ELSE IF(n .EQ. 5)THEN
  uvout(1:3,1:3) = d5(uv(1:3,1:3), s)
ELSE IF(n .EQ. 6)THEN
  uvout(1:3,1:3) = d6(uv(1:3,1:3), s)
ELSE IF(n .EQ. 7)THEN
  uvout(1:3,1:3) = d7(uv(1:3,1:3), s)
ELSE IF(n .EQ. 8)THEN
  uvout(1:3,1:3) = d8(uv(1:3,1:3), s)
ELSE IF(n .EQ. 9)THEN
  uvout(1:3,1:3) = d9(uv(1:3,1:3), s)
END IF
END FUNCTION dn