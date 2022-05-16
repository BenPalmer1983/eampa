



FUNCTION vecvol(uv) RESULT (vol)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal) :: vol
!############################################################
REAL(kind=DoubleReal) :: bxc (1:3)
!############################################################
bxc(1) = uv(2, 2)*uv(3, 3)-uv(2, 3)*uv(3, 2)
bxc(2) = uv(2, 3)*uv(3, 1)-uv(2, 1)*uv(3, 3)
bxc(3) = uv(2, 1)*uv(3, 2)-uv(2, 2)*uv(3, 1)

vol = uv(1,1)*bxc(1)+uv(1,2)*bxc(2)+uv(1,3)*bxc(3)

END FUNCTION vecvol


