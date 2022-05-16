FUNCTION ctd1(uv, s) RESULT (uvout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal), INTENT(IN) :: s
REAL(kind=DoubleReal) :: uvout(1:3, 1:3)
!############################################################
REAL(kind=DoubleReal) :: d(1:3, 1:3)
!############################################################
d(1:3, 1:3) = 0.0D0
d(1,1) = 1.0D0 + s
d(2,2) = 1.0D0 - (1.0D0 - sqrt(1.0D0 / (1.0D0 + s)))
d(3,3) = 1.0D0 - (1.0D0 - sqrt(1.0D0 / (1.0D0 + s)))
uvout(1:3, 1:3) = MATMUL(d(1:3, 1:3), uv(1:3, 1:3))
END FUNCTION ctd1

FUNCTION ctd1_inv(uv, s) RESULT (uvout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal), INTENT(IN) :: s
REAL(kind=DoubleReal) :: uvout(1:3, 1:3)
!############################################################
REAL(kind=DoubleReal) :: d(1:3, 1:3)
!############################################################
d(1:3, 1:3) = 0.0D0
d(1,1) = 1.0D0 + s
d(2,2) = 1.0D0
d(3,3) = 1.0D0
d = minverse3(d)
uvout(1:3, 1:3) = MATMUL(d(1:3, 1:3), uv(1:3, 1:3))
END FUNCTION ctd1_inv