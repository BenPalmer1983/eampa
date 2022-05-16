MODULE optimizer

USE kinds

IMPLICIT NONE

!INCLUDE "minimizer.interfaces.f90"
!INCLUDE "minimizer.globals.f90"

CONTAINS

INCLUDE "optimizer.newtongauss.f90"
INCLUDE "callback.f90"




END MODULE optimizer