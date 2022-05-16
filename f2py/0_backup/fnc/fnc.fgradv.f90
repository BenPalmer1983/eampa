! VECTOR SUBROUTINE
SUBROUTINE fgradv(f_name, r, p, p_fixed, dydr)
!############################################################
! PAIR SPLINE
IMPLICIT NONE
!############################################################
CHARACTER(LEN=*), INTENT(IN) :: f_name
REAL(kind=DoubleReal), INTENT(IN) :: r(:)
REAL(kind=DoubleReal), INTENT(IN) :: p(:)
REAL(kind=DoubleReal), INTENT(IN) :: p_fixed(:)
REAL(kind=DoubleReal), INTENT(OUT) :: dydr(1:SIZE(r,1))
INTEGER(kind=StandardInteger) :: n
LOGICAL :: result
!############################################################

IF(str_compare(f_name, "zero"))THEN
  CALL zero_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "ackland_embedding"))THEN
  CALL ackland_embedding_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "ackland_spline"))THEN
  CALL ackland_spline_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "buckingham"))THEN
  CALL buckingham_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "fs_embedding"))THEN
  CALL fs_embedding_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "knot_spline"))THEN
  CALL knot_spline_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "lennard_jones"))THEN
  CALL lennard_jones_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "mendelev_embedding"))THEN
  CALL mendelev_embedding_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "mishin_density"))THEN
  CALL mishin_density_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "morse"))THEN
  CALL morse_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "slater_4s"))THEN
  CALL slater_4s_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "spline"))THEN
  CALL spline_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "spline_embedding"))THEN
  CALL spline_embedding_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "triple_embedding"))THEN
  CALL triple_embedding_grad_v(r, p, p_fixed, dydr)
ELSE IF(str_compare(f_name, "zbl"))THEN
  CALL zbl_grad_v(r, p, p_fixed, dydr)








END IF

END SUBROUTINE fgradv
