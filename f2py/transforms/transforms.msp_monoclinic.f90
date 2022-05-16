! Mehl Singh Papaconstantopoulos strains 1993
! Properties of ordered intermetallic alloys: first-principles
! and approximate methods
!
! Mehl Klein Papaconstantopoulos 1993
! First principles calculations of elastic properties of metals
!
! Orthorhombic strain to calculate C44
!    
! Includes identity matrix:
!    I + e
! 

FUNCTION msp_monoclinic(uv, s) RESULT (uvout)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: uv(1:3, 1:3)
REAL(kind=DoubleReal), INTENT(IN) :: s
REAL(kind=DoubleReal) :: uvout(1:3, 1:3)
!############################################################
uvout(1:3, 1:3) = 0.0D0
uvout(1,1) = 1.0D0
uvout(2,2) = 1.0D0
uvout(3,3) = 1.0D0 + (s**2 / (4.0D0 - s**2))
uvout(1,2) = s / 2.0D0
uvout(2,1) = s / 2.0D0
END FUNCTION msp_monoclinic