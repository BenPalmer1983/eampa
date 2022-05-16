SUBROUTINE gcoords(cin, labelsin, a0, uv, cout, labelsout, halo)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: cin(:,:)
INTEGER(kind=StandardInteger), INTENT(IN) :: labelsin(:)
REAL(kind=DoubleReal), INTENT(IN) :: a0
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3,1:3)
REAL(kind=DoubleReal), INTENT(OUT) :: cout(1:27*SIZE(cin,1),1:3)
INTEGER(kind=StandardInteger), INTENT(OUT) :: labelsout(1:27*SIZE(cin,1), 1:2)
INTEGER(kind=StandardInteger), INTENT(OUT) :: halo(1:27*SIZE(cin,1))
!############################################################
INTEGER(kind=StandardInteger) :: m, n, i, j, k
REAL(kind=DoubleReal) :: tr(1:3,1:3)
REAL(kind=DoubleReal) :: s(1:3)
!############################################################
tr = a0 * uv
m = 0
DO i=-1,1
  DO j=-1,1 
    DO k=-1,1
      s(1) = 1.0D0 * i
      s(2) = 1.0D0 * j
      s(3) = 1.0D0 * k
      DO n = 1, SIZE(cin,1)
        m = m + 1
        cout(m,:) = MATMUL(tr, cin(n,1:3) + s(1:3))
        labelsout(m, 1) = labelsin(n)
        labelsout(m, 2) = n
        IF(i .EQ. 0 .AND. j .EQ. 0 .AND. k .EQ. 0)THEN
          halo(m) = 0
        ELSE
          halo(m) = 1
        END IF
      END DO
    END DO
  END DO
END DO
END SUBROUTINE gcoords
