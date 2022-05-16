SUBROUTINE sls_solve_ls(lower, upper, y, x)
!############################################################
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: lower
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: upper
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: y
REAL(kind=DoubleReal), DIMENSION(:), INTENT(OUT) :: x
!############################################################
REAL(kind=DoubleReal), DIMENSION(1:size(lower,1)) :: b
INTEGER(kind=StandardInteger) :: m_size
INTEGER(kind=StandardInteger) :: row, col
REAL(kind=DoubleReal) :: sum_term, partial
INTEGER :: i, partial_Sum, total_Sum
!############################################################

m_size = SIZE(lower,1)

b = 0.0D0

!# Solve Lb = y
  
  DO row = 1, m_size
    sum_term = 0.0D0
    DO col = 1, row - 1
      sum_term = sum_term + lower(row, col) * b(col)
    END DO
    b(row) = (y(row) - sum_term) / lower(row, row)
  END DO

  !# Solve Ux=b
  DO row = m_size, 1, -1 
    sum_term = 0.0D0
    DO col = m_size, row + 1, -1      
      sum_term = sum_term + upper(row, col) * x(col)
    END DO
    x(row) = (b(row) - sum_term) / upper(row, row)
  END DO


END SUBROUTINE sls_solve_ls