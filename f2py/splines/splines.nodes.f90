SUBROUTINE spline_nodes(rspline, poly, na, nb, xyarr)
!############################################################
! p coefficients
! pf r cutoffs
! they must be the same size
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: rspline
INTEGER(kind=StandardInteger), INTENT(IN) :: poly
REAL(kind=DoubleReal), INTENT(IN) :: na(:)
REAL(kind=DoubleReal), INTENT(IN) :: nb(:)
REAL(kind=DoubleReal), INTENT(INOUT) :: xyarr(1:4)
!############################################################
REAL(kind=DoubleReal) :: xmat(1:10,1:10)
REAL(kind=DoubleReal) :: ymat(1:10)
REAL(kind=DoubleReal) :: c(1:10)
!############################################################
xyarr(:) = 0.0D0
IF(poly .EQ. 3)THEN
  ! X VALUES IN p_fixed
  ! Y VALUES IN p
  ! DERIV VALUES IN dydr  
  xmat(1,1) = 1.0
  xmat(1,2) = na(1)
  xmat(1,3) = na(1)**2
  xmat(1,4) = na(1)**3
  xmat(2,1) = 1.0
  xmat(2,2) = nb(1)
  xmat(2,3) = nb(1)**2
  xmat(2,4) = nb(1)**3
  xmat(3,1) = 0.0
  xmat(3,2) = 1.0D0
  xmat(3,3) = 2.0D0 * na(1)
  xmat(3,4) = 3.0D0 * na(1)**2
  xmat(4,1) = 0.0
  xmat(4,2) = 1.0D0
  xmat(4,3) = 2.0D0 * nb(1)
  xmat(4,4) = 3.0D0 * nb(1)**2
  
  ymat(1) = na(2)
  ymat(2) = nb(2)
  ymat(3) = na(3)
  ymat(4) = nb(3)
    
  CALL sls_solve(xmat(1:4,1:4), ymat(1:4), c(1:4))

  xyarr(1) = rspline
  xyarr(2) = c(1) + c(2) * rspline + c(3) * rspline**2 + c(4) * rspline**3
  xyarr(3) = c(2) + 2 * c(3) * rspline + 3 * c(4) * rspline**2
  xyarr(4) = 2 * c(3) + 6 * c(4) * rspline

ELSE IF(poly .EQ. 5)THEN
  ! X VALUES IN p_fixed
  ! Y VALUES IN p
  ! DERIV VALUES IN dydr  
  xmat(1,1) = 1.0D0
  xmat(1,2) = na(1)
  xmat(1,3) = na(1)**2
  xmat(1,4) = na(1)**3
  xmat(1,5) = na(1)**4
  xmat(1,6) = na(1)**5
  xmat(2,1) = 1.0D0
  xmat(2,2) = nb(1)
  xmat(2,3) = nb(1)**2
  xmat(2,4) = nb(1)**3
  xmat(2,5) = nb(1)**4
  xmat(2,6) = nb(1)**5
  xmat(3,1) = 0.0D0
  xmat(3,2) = 1.0D0
  xmat(3,3) = 2.0D0 * na(1)
  xmat(3,4) = 3.0D0 * na(1)**2
  xmat(3,5) = 4.0D0 * na(1)**3
  xmat(3,6) = 5.0D0 * na(1)**4
  xmat(4,1) = 0.0D0
  xmat(4,2) = 1.0D0
  xmat(4,3) = 2.0D0 * nb(1)
  xmat(4,4) = 3.0D0 * nb(1)**2  
  xmat(4,5) = 4.0D0 * nb(1)**3  
  xmat(4,6) = 5.0D0 * nb(1)**4  
  xmat(5,1) = 0.0D0
  xmat(5,2) = 0.0D0
  xmat(5,3) = 2.0D0 
  xmat(5,4) = 6.0D0 * na(1)
  xmat(5,5) = 12.0D0 * na(1)**2
  xmat(5,6) = 20.0D0 * na(1)**3
  xmat(6,1) = 0.0D0
  xmat(6,2) = 0.0D0
  xmat(6,3) = 2.0D0 
  xmat(6,4) = 6.0D0 * nb(1)
  xmat(6,5) = 12.0D0 * nb(1)**2
  xmat(6,6) = 20.0D0 * nb(1)**3
  
  ymat(1) = na(2)
  ymat(2) = nb(2)
  ymat(3) = na(3)
  ymat(4) = nb(3)
  ymat(5) = na(4)
  ymat(6) = nb(4)
    
  CALL sls_solve(xmat(1:6,1:6), ymat(1:6), c(1:6))

  xyarr(1) = rspline
  xyarr(2) = c(1) + c(2) * rspline + c(3) * rspline**2 + c(4) * rspline**3 + &
                 c(5) * rspline**4 + c(6) * rspline**5
  xyarr(3) = c(2) + 2 * c(3) * rspline + 3 * c(4) * rspline**2 + &
                 4 * c(5) * rspline**3 + 5 * c(6) * rspline**4
  xyarr(4) = 2 * c(3) + 6 * c(4) * rspline + &
                 12 * c(5) * rspline**2 + 20 * c(6) * rspline**3
END IF


END SUBROUTINE spline_nodes