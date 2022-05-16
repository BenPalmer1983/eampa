FUNCTION orthorhombic(uv, s) RESULT (uvout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal), INTENT(IN) :: s
REAL(kind=DoubleReal) :: uvout(1:3, 1:3)
!############################################################
uvout(1:3, 1:3) = (1.0D0 + s) * uv(1:3, 1:3)
END FUNCTION orthorhombic