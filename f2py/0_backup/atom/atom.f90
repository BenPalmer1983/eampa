!#############################################################
!
!#############################################################

MODULE atom

USE kinds
USE interp, ONLY: interpolate
USE fnc, ONLY: f, fgrad
USE math, ONLY: vecvol, minverse3

IMPLICIT NONE

INCLUDE "atom.globals.f90"
INCLUDE "atom.interfaces.f90"

CONTAINS

INCLUDE "atom.coords.f90"
INCLUDE "atom.gcoords.f90"
INCLUDE "atom.nl.f90"
INCLUDE "atom.nl_change_cell.f90"
INCLUDE "atom.nl_update_coords.f90"
INCLUDE "atom.potential.f90"
INCLUDE "atom.e.f90"
INCLUDE "atom.ef.f90"
INCLUDE "atom.efs.f90"
INCLUDE "atom.verlet_step.f90"
INCLUDE "atom.relax.f90"







END MODULE atom
