! VECTOR SUBROUTINE
SUBROUTINE f(f_name, r, p, p_fixed, y)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: f_name
REAL(kind=DoubleReal), INTENT(IN) :: r
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y
INTEGER(kind=StandardInteger) :: n
LOGICAL :: result
!############################################################

IF(str_compare(f_name, "zero"))THEN
  CALL zero(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_embedding"))THEN
  CALL ackland_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_spline"))THEN
  CALL ackland_spline(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "buckingham"))THEN
  CALL buckingham(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "fs_embedding"))THEN
  CALL fs_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "knot_spline"))THEN
  CALL knot_spline(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "lennard_jones"))THEN
  CALL lennard_jones(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "mendelev_embedding"))THEN
  CALL mendelev_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "mishin_density"))THEN
  CALL mishin_density(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "morse"))THEN
  CALL morse(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "slater_4s"))THEN
  CALL slater_4s(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "spline"))THEN
  CALL spline(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "spline_embedding"))THEN
  CALL spline_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "triple_embedding"))THEN
  CALL triple_embedding(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "zbl"))THEN
  CALL zbl(r, p, p_fixed, y)




END IF

END SUBROUTINE f
