

SUBROUTINE coords(cin, a0, uv, cout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: cin(:,:)
REAL(kind=DoubleReal), INTENT(IN) :: a0
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3,1:3)
REAL(kind=DoubleReal), INTENT(OUT) :: cout(1:SIZE(cin,1),1:3)
!############################################################
INTEGER(kind=StandardInteger) :: n
REAL(kind=DoubleReal) :: tr(1:3,1:3)
!############################################################

tr = a0 * uv
DO n = 1, SIZE(cin,1)
  cout(n,:) = MATMUL(tr, cin(n,:))
  print *, cout(n,:)
END DO


END SUBROUTINE coords

