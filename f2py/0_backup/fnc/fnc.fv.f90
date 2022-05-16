! VECTOR SUBROUTINE
SUBROUTINE fv(f_name, r, p, p_fixed, y)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: f_name
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: y(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
LOGICAL :: result
!############################################################

IF(str_compare(f_name, "zero"))THEN
  CALL zero_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_embedding"))THEN
  CALL ackland_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "ackland_spline"))THEN
  CALL ackland_spline_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "buckingham"))THEN
  CALL buckingham_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "fs_embedding"))THEN
  CALL fs_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "knot_spline"))THEN
  CALL knot_spline_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "lennard_jones"))THEN
  CALL lennard_jones_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "mendelev_embedding"))THEN
  CALL mendelev_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "mishin_density"))THEN
  CALL mishin_density_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "morse"))THEN
  CALL morse_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "slater_4s"))THEN
  CALL slater_4s_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "spline"))THEN
  CALL spline_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "spline_embedding"))THEN
  CALL spline_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "triple_embedding"))THEN
  CALL triple_embedding_v(r, p, p_fixed, y)
ELSE IF(str_compare(f_name, "zbl"))THEN
  CALL zbl_v(r, p, p_fixed, y)






END IF

END SUBROUTINE fv
