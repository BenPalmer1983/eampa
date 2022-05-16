SUBROUTINE sls_mpivot(A, y)
!############################################################
REAL(kind=DoubleReal), INTENT(INOUT), DIMENSION(:,:) :: A
REAL(kind=DoubleReal), INTENT(INOUT), DIMENSION(:) :: y
!############################################################
REAL(kind=DoubleReal), DIMENSION(1:SIZE(A,1),1:SIZE(A,2)) :: PivotMatrix
!############################################################
INTEGER(kind=StandardInteger) :: n, j
INTEGER(kind=StandardInteger) :: max_row
REAL(kind=DoubleReal) :: max_row_val
REAL(kind=DoubleReal) :: temp
!############################################################
Call sls_makepivot(A, PivotMatrix)
A = matmul(PivotMatrix, A)
y = matmul(PivotMatrix, y)
!############################################################
END SUBROUTINE sls_mpivot


SUBROUTINE sls_makepivot(A, PivotMatrix)
!############################################################
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(INOUT) :: A
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(INOUT) :: PivotMatrix
!############################################################
INTEGER(kind=StandardInteger) :: n, j
INTEGER(kind=StandardInteger) :: max_row
REAL(kind=DoubleReal) :: max_row_val
REAL(kind=DoubleReal) :: temp
!############################################################
!# Get identity
Call sls_midentity(PivotMatrix)
Do j = 1, size(A,1)
  max_row = j
  max_row_val = 0.0D0
  Do n = j, size(A,2)
    If (abs(A(n,j)) .GT. max_row_val) Then
      max_row_val = abs(A(n,j))
      max_row = n
    End If
  End Do
  If (j .NE. max_row) Then
    !# Swap row
    Do n = 1, size(A,2)
      temp = PivotMatrix(j,n)
      PivotMatrix(j,n) = PivotMatrix(max_row,n)
      PivotMatrix(max_row,n) = temp
    End Do
  End If
End Do
!############################################################
END SUBROUTINE sls_makepivot


SUBROUTINE sls_midentity(A)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(INOUT) :: A(:,:)
!############################################################
INTEGER(kind=StandardInteger) :: c
!############################################################
A = 0.0D0
DO c=1,size(A,1)
  A(c,c) = 1.0D0
END DO
!############################################################
END SUBROUTINE sls_midentity