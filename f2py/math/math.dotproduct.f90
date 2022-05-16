FUNCTION dotproduct(a, b) RESULT (c)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: b(1:3)
REAL(kind=DoubleReal) :: c
!############################################################
c = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
END FUNCTION dotproduct