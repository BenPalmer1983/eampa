MODULE math

USE kinds

IMPLICIT NONE

INCLUDE "math.interfaces.f90"
INCLUDE "math.globals.f90"

CONTAINS

INCLUDE "math.vecvol.f90"
INCLUDE "math.minverse.f90"
INCLUDE "math.minverse3.f90"
INCLUDE "math.matmulcell.f90"
INCLUDE "math.msquare.f90"
INCLUDE "math.midentity.f90"
INCLUDE "math.crossproduct.f90"
INCLUDE "math.dotproduct.f90"
INCLUDE "math.tripleproduct.f90"
INCLUDE "math.cellvol.f90"


END MODULE math