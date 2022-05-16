SUBROUTINE matmulcell(a, b)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a(1:3,1:3)
REAL(kind=DoubleReal), INTENT(INOUT) :: b(1:3)
!############################################################
REAL(kind=DoubleReal) :: t(1:3)
!############################################################
t(1:3) = b(1:3)
b(1) = a(1,1) * t(1) + a(1,2) * t(2) + a(1,3) * t(3) 
b(2) = a(2,1) * t(1) + a(2,2) * t(2) + a(2,3) * t(3) 
b(3) = a(3,1) * t(1) + a(3,2) * t(2) + a(3,3) * t(3) 
!############################################################
END SUBROUTINE matmulcell