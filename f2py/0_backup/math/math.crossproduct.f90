FUNCTION crossproduct(a, b) RESULT (c)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: b(1:3)
!############################################################
REAL(kind=DoubleReal) :: c(1:3)
!############################################################
! Calculate cross product
c(1) = a(2)*b(3)-a(3)*b(2)
c(2) = a(3)*b(1)-a(1)*b(3)
c(3) = a(1)*b(2)-a(2)*b(1)
!############################################################
END FUNCTION crossproduct