SUBROUTINE sls_lu_decomp(A, lower, upper)
!############################################################
REAL(kind=DoubleReal), INTENT(IN), DIMENSION(:,:) :: A
REAL(kind=DoubleReal), INTENT(OUT), DIMENSION(:,:) :: lower
REAL(kind=DoubleReal), INTENT(OUT), DIMENSION(:,:) :: upper
!############################################################
INTEGER(kind=StandardInteger) :: n, i, j, k
REAL(kind=DoubleReal) :: sum_a, sum_b
!############################################################
lower = 0.0D0
upper = 0.0D0
n = size(A,1)
Do j = 1, n
  lower(j,j) = 1.0D0
  Do i = 1, j
    sum_a = 0.0D0
    Do k = 1, i
      sum_a = sum_a + upper(k,j) * lower(i,k)
    End Do
    upper(i,j) = A(i,j) - sum_a
  End Do
  Do i = j, n
    sum_b = 0.0D0
    Do k = 1, j-1
      sum_b = sum_b + upper(k,j) * lower(i,k)
    End Do
    lower(i,j) = (A(i,j) - sum_b) / upper(j,j)
  End Do
End Do
!############################################################
END SUBROUTINE sls_lu_decomp


