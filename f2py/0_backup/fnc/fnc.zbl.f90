

SUBROUTINE zbl(r, p, p_fixed, y)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
!############################################################
REAL(kind=DoubleReal) :: qa
REAL(kind=DoubleReal) :: qb
REAL(kind=DoubleReal) :: rs
REAL(kind=DoubleReal) :: ep
REAL(kind=DoubleReal) :: x
!############################################################
IF(r .LE. 0.0D0)THEN
  y = 1.0D10
ELSE
  qa = p_fixed(1)
  qb = p_fixed(2)

  ! MENDELEV ACKLAND 2004 Development of an interatomic potential for phosphorus impurities in α-iron
  rs = 0.4683766 / (qa**(2.0D0/3.0D0)+qb**(2.0D0/3.0D0))
  x = r / rs
  ep = 0.1818D0 * exp(-3.2D0 * x) &
       + 0.5099D0 * exp(-0.9423D0 * x) &
       + 0.2802 * exp(-0.4029D0 * x) &
       + 0.02817 * exp(-0.2016D0 * x)   
  y = (qa * qb * ep) / r

  ! ACKLAND 1997 Computer simulation of point defect properties in dilute Fe—Cu alloy using a many-body interatomic potential
  !qa = p(1)
  !qb = p(2)
  !rs = 0.46834486 / sqrt(qa**(2.0D0/3.0D0)+qb**(2.0D0/3.0D0))
  !x = r / rs
  !ep = 0.1818D0 * exp(-3.2D0 * x) 
  !ep = ep + 0.5099D0 * exp(-0.9423D0 * x) 
  !ep = ep + 0.2802 * exp(-0.4029D0 * x) 
  !ep = ep + 0.02817 * exp(-0.2016D0 * x) 
  !y = ((qa * qb) / r) * ep
END IF
END SUBROUTINE zbl




SUBROUTINE zbl_v(r, p, p_fixed, y)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL zbl(r(n), p, p_fixed,  y(n))
END DO
END SUBROUTINE zbl_v







SUBROUTINE zbl_dydr(r, p, p_fixed, dydr)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr
!############################################################
REAL(kind=DoubleReal) :: qa
REAL(kind=DoubleReal) :: qb
REAL(kind=DoubleReal) :: rs
REAL(kind=DoubleReal) :: ep
REAL(kind=DoubleReal) :: epp
!############################################################
IF(r .LE. 0.0D0)THEN
  dydr = -1.0D10
ELSE
  qa = p(1)
  qb = p(2)
  rs = 0.4683766 / (qa**(2.0D0/3.0D0)+qb**(2.0D0/3.0D0))
  
  
  ep =      0.1818D0 * exp(r * (-3.2D0/rs)) 
  ep = ep + 0.5099D0 * exp(r * (-0.9423D0/rs)) 
  ep = ep + 0.2802 * exp(r * (-0.4029D0/rs)) 
  ep = ep + 0.02817 * exp(r * (-0.2016D0/rs)) 
  
  
  epp =      (-3.2D0/rs) * 0.1818D0 * exp(r * (-3.2D0/rs)) 
  epp = ep + (-0.9423D0/rs) * 0.5099D0 * exp(r * (-0.9423D0/rs)) 
  epp = ep + (-0.4029D0/rs) * 0.2802 * exp(r * (-0.4029D0/rs)) 
  epp = ep + (-0.2016D0/rs) * 0.02817 * exp(r * (-0.2016D0/rs)) 
  
  dydr = ((qa * qb) / r) * epp - ((qa * qb) / (r**2)) * ep
END IF
END SUBROUTINE zbl_dydr



SUBROUTINE zbl_dydr_v(r, p, p_fixed, y)
!############################################################
! ZBL POTENTIAL
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(1:2)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(1:1)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
!############################################################
INTEGER(kind=StandardInteger) :: n
!############################################################
DO n = 1, SIZE(r,1)
  CALL zbl_dydr(r(n), p, p_fixed,  y(n))
END DO
END SUBROUTINE zbl_dydr_v

