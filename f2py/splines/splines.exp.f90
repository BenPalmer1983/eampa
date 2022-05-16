! 
SUBROUTINE splines_exp(na, nb, p)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: na(1:4)
REAL(kind=DoubleReal), INTENT(IN) :: nb(1:4)
REAL(kind=DoubleReal), INTENT(INOUT) :: p(1:4)
!############################################################
REAL(kind=DoubleReal) :: rss
!############################################################

p = 0.0D0
g_na(:) = na(:)
g_nb(:) = nb(:)


CALL splines_exp_rss(p, rss)
print *, rss





!############################################################
END SUBROUTINE splines_exp


SUBROUTINE splines_exp_rss(p, rss)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: p(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: rss
!############################################################
REAL(kind=DoubleReal) :: na_test(1:4)
REAL(kind=DoubleReal) :: nb_test(1:4)
!############################################################
rss = 0.0D0
na_test = 0.0D0
nb_test = 0.0D0

na_test(1) = g_na(1)
na_test(2) = exp(p(1) + p(1)*g_na(2) + p(3)*g_na(1)**2 + p(3)*g_na(1)**3)
na_test(3) = (p(1) + 2.0D0 * p(3)*g_na(1) + 3.0D0 * p(3)*g_na(1)**3) * exp(p(1) + p(1)*g_na(2) + p(3)*g_na(1)**2 + p(3)*g_na(1)**3)
na_test(4) = 0.0D0

nb_test(1) = g_nb(1)
nb_test(2) = exp(p(1) + p(1)*g_na(2) + p(3)*g_na(1)**2 + p(3)*g_na(1)**3)
nb_test(3) = (p(1) + 2.0D0 * p(3)*g_nb(1) + 3.0D0 * p(3)*g_nb(1)**3) * exp(p(1) + p(1)*g_nb(2) + p(3)*g_nb(1)**2 + p(3)*g_nb(1)**3)
nb_test(4) = 0.0D0

rss = SUM((na_test(2:4)-g_na(2:4))**2) + SUM((nb_test(2:4)-g_nb(2:4))**2)
print *, nb_test(:)

!############################################################
END SUBROUTINE splines_exp_rss





SUBROUTINE splines_exp_fit(p_in, p_out)
!############################################################
IMPLICIT NONE
!############################################################
REAL(kind=DoubleReal), INTENT(IN) :: p_in(1:4)
REAL(kind=DoubleReal), INTENT(OUT) :: p_out(1:4)
!############################################################
REAL(kind=DoubleReal) :: J(1:4, 1:3)




!############################################################
END SUBROUTINE splines_exp_fit






