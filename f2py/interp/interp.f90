!#############################################################
!#  HOW TO USE:
!#
!#
!# interp.interp4(x, x_points, y_points, y)
!# interp.interp4dydx(x, x_points, y_points, y)
!#
!#
!# 
!# 
!#
!#
!#
!#
!#############################################################

MODULE interp

USE kinds

IMPLICIT NONE

INCLUDE "interp.globals.f90"
INCLUDE "interp.interfaces.f90"

CONTAINS

INCLUDE "interp.interpn.f90"
INCLUDE "interp.interpndydx.f90"
INCLUDE "interp.interpndydxn.f90"

INCLUDE "interp.interpolate.f90"
INCLUDE "interp.interpolate_row.f90"
INCLUDE "interp.fill.f90"


!INCLUDE "interp.interpolate_dydxn.f90"
!INCLUDE "interp.trap.f90"

!INCLUDE "interp.fill.f90"

!INCLUDE "interp.search_x.f90"



!INCLUDE "interp.speed_test.f90"

END MODULE interp
