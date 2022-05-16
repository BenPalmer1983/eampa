!#############################################################
!#  HOW TO USE:
!#
!#
!#                                                        
!#
!#############################################################

MODULE fnc

USE kinds
USE sls, ONLY: sls_solve

IMPLICIT NONE

INCLUDE "fnc.globals.f90"
INCLUDE "fnc.interfaces.f90"

CONTAINS


INCLUDE "fnc.f.f90"
INCLUDE "fnc.fv.f90"
INCLUDE "fnc.fgrad.f90"
INCLUDE "fnc.fgradv.f90"
INCLUDE "fnc.str_compare.f90"
INCLUDE "fnc.fzbl.f90"
INCLUDE "fnc.interp.f90"
INCLUDE "fnc.step.f90"
INCLUDE "fnc.heaviside.f90"
INCLUDE "fnc.generic.f90"
INCLUDE "fnc.minverse.f90"
INCLUDE "fnc.cutoff.f90"


! AVAILABLE FUNCTIONS
INCLUDE "fnc.zero.f90"
INCLUDE "fnc.zero_grad.f90"
INCLUDE "fnc.ackland_embedding.f90"
INCLUDE "fnc.ackland_embedding_grad.f90"
INCLUDE "fnc.ackland_spline.f90"
INCLUDE "fnc.ackland_spline_grad.f90"
INCLUDE "fnc.buckingham.f90"
INCLUDE "fnc.buckingham_grad.f90"
INCLUDE "fnc.fs_embedding.f90"
INCLUDE "fnc.fs_embedding_grad.f90"
INCLUDE "fnc.knot_spline.f90"
INCLUDE "fnc.knot_spline_grad.f90"
INCLUDE "fnc.lennard_jones.f90"
INCLUDE "fnc.lennard_jones_grad.f90"
INCLUDE "fnc.mendelev_embedding.f90"
INCLUDE "fnc.mendelev_embedding_grad.f90"
INCLUDE "fnc.mishin_density.f90"
INCLUDE "fnc.mishin_density_grad.f90"
INCLUDE "fnc.morse.f90"
INCLUDE "fnc.morse_grad.f90"
INCLUDE "fnc.slater_4s.f90"
INCLUDE "fnc.slater_4s_grad.f90"
INCLUDE "fnc.spline.f90"
INCLUDE "fnc.spline_grad.f90"
INCLUDE "fnc.spline_embedding.f90"
INCLUDE "fnc.spline_embedding_grad.f90"
INCLUDE "fnc.triple_embedding.f90"
INCLUDE "fnc.triple_embedding_grad.f90"
INCLUDE "fnc.zbl.f90"
INCLUDE "fnc.zbl_grad.f90"





END MODULE fnc
