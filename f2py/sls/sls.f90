MODULE sls

USE kinds

IMPLICIT NONE

INCLUDE "sls.interfaces.f90"
INCLUDE "sls.globals.f90"

CONTAINS

INCLUDE "sls.solve.f90"
INCLUDE "sls.pivot.f90"
INCLUDE "sls.lu_decomp.f90"
INCLUDE "sls.solve_ls.f90"



END MODULE sls