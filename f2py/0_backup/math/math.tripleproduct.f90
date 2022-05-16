FUNCTION tripleproduct(a, b, c) RESULT (tp)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: b(1:3)
REAL(kind=DoubleReal), INTENT(IN) :: c(1:3)
!############################################################
REAL(kind=DoubleReal) :: tp
!############################################################
tp = dotproduct(a,crossproduct(b, c))
!############################################################
END FUNCTION tripleproduct