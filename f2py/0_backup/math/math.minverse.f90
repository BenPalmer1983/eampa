!============================================================
! Original subroutine found online.  
! Author Alex G. December 2009
!===========================================================


SUBROUTINE minverse(a, c)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: a(:,:)
REAL(kind=DoubleReal), INTENT(OUT) :: c(1:SIZE(a,1),1:SIZE(a,1))
!############################################################
REAL(kind=DoubleReal) :: aw(SIZE(a,1),SIZE(a,1))
REAL(kind=DoubleReal) :: L(SIZE(a,1),SIZE(a,1))
REAL(kind=DoubleReal) :: U(SIZE(a,1),SIZE(a,1))
REAL(kind=DoubleReal) :: b(SIZE(a,1))
REAL(kind=DoubleReal) :: d(SIZE(a,1))
REAL(kind=DoubleReal) :: x(SIZE(a,1))
REAL(kind=DoubleReal) :: coeff
INTEGER(kind=StandardInteger) :: n, i, j, k
!############################################################
IF(SIZE(a,1) .NE. SIZE(a,2))THEN
  RETURN
END IF
aw(:,:) = a(:,:)
L(:,:) = 0.0D0
U(:,:) = 0.0D0
b(:) = 0.0D0
n = SIZE(a,1)
! Forward elimination
DO k=1, n-1
   DO i=k+1,n
      coeff = aw(i,k)/aw(k,k)
      L(i,k) = coeff
      DO j=k+1,n
         aw(i,j) = aw(i,j)-coeff*aw(k,j)
      END DO
   END DO
END DO

! Prepare L and U matrices 
DO i=1,n
  L(i,i) = 1.0
END DO
! U matrix is the upper triangular part of A
DO j=1,n
  U(1:j,j) = aw(1:j,j)
END DO

! Compute columns of the inverse matrix C
DO k = 1,n
  b(k) = 1.0D0
  d(1) = b(1)
! Solve Ld=b using the forward substitution
  DO i=2,n
    d(i)=b(i)
    DO j=1, i-1
      d(i) = d(i) - L(i,j)*d(j)
    END DO
  END DO
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n) / U(n,n)
  DO i = n-1,1,-1
    x(i) = d(i)
    DO j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    END DO
    x(i) = x(i)/u(i,i)
  END DO
! Step 3c: fill the solutions x(n) into column k of C
  c(1:n,k) = x(1:n)
  b(k)=0.0D0
END DO
END SUBROUTINE minverse