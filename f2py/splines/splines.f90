MODULE splines

USE kinds
USE sls, ONLY: sls_solve

IMPLICIT NONE

INCLUDE "splines.interfaces.f90"
INCLUDE "splines.globals.f90"

CONTAINS

INCLUDE "splines.nodes.f90"
INCLUDE "splines.exp.f90"




END MODULE splines