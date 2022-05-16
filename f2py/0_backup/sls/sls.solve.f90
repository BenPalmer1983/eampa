!#
!#   sls_solve and solve are the same
!#
!#    A x = y
!#    Solve to find x

SUBROUTINE sls_solve(A, y, x)
!###########################################
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: A
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: y
REAL(kind=DoubleReal), DIMENSION(1:size(y,1)), INTENT(OUT) :: x
!###########################################
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: A_working
REAL(kind=DoubleReal), DIMENSION(1:size(y,1)) :: y_working
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: lower
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: upper
!###########################################
A_working = A
y_working = y
x = 0.0D0
upper = 0.0D0
lower = 0.0D0
CALL sls_mpivot(A_working, y_working)
CALL sls_lu_decomp(A, lower, upper)
CALL sls_solve_ls(lower, upper, y, x)
END SUBROUTINE sls_solve

SUBROUTINE solve(A, y, x)
!###########################################
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: A
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: y
REAL(kind=DoubleReal), DIMENSION(1:size(y,1)), INTENT(OUT) :: x
!###########################################
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: A_working
REAL(kind=DoubleReal), DIMENSION(1:size(y,1)) :: y_working
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: lower
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: upper
!###########################################
A_working = A
y_working = y
x = 0.0D0
upper = 0.0D0
lower = 0.0D0
CALL sls_mpivot(A_working, y_working)
CALL sls_lu_decomp(A, lower, upper)
CALL sls_solve_ls(lower, upper, y, x)
END SUBROUTINE solve


SUBROUTINE solve_quick(A, y, x)
!###########################################
REAL(kind=DoubleReal), DIMENSION(:,:), INTENT(IN) :: A
REAL(kind=DoubleReal), DIMENSION(:), INTENT(IN) :: y
REAL(kind=DoubleReal), DIMENSION(1:size(y,1)), INTENT(OUT) :: x
!###########################################
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: A_working
REAL(kind=DoubleReal), DIMENSION(1:size(y,1)) :: y_working
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: lower
REAL(kind=DoubleReal), DIMENSION(1:size(A,1),1:size(A,2)) :: upper
!###########################################
A_working = A
y_working = y
x = 0.0D0
upper = 0.0D0
lower = 0.0D0
CALL sls_lu_decomp(A, lower, upper)
CALL sls_solve_ls(lower, upper, y, x)
END SUBROUTINE solve_quick