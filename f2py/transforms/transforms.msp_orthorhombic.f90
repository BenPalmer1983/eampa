! Mehl Singh Papaconstantopoulos strains 1993
! Properties of ordered intermetallic alloys: first-principles
! and approximate methods
!
! Orthorhombic strain to calculate C11 C12
!    
! Includes identity matrix:
!    I + e
! 

FUNCTION msp_orthorhombic(uv, s) RESULT (uvout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal), INTENT(IN) :: s
REAL(kind=DoubleReal) :: uvout(1:3, 1:3)
!############################################################
uvout(1:3, 1:3) = 0.0D0
uvout(1,1) = (1.0D0 + s)
uvout(2,2) = (1.0D0 - s)
uvout(3,3) = (1.0D0 + s**2 / (1.0D0 - s**2))
END FUNCTION msp_orthorhombic