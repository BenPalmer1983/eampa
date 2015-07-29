Module finalise

! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: finalise
! Updated: 18th May 2015
! --------------------------------------------------------------!
! Description:
! Makes directories used when program runs
! Creates the output files with headers
! Defines the ProgramTime() function
! --------------------------------------------------------------!

! Setup Modules
  Use kinds
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use globals        ! declare all globals
  Use initialise     ! initialise program

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Subroutines
  Public :: runFinaliseEval        !Subroutine

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE SUBROUTINES                                                      !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

  contains

! Run all the input subroutines
  Subroutine runFinaliseEval()
! Internal subroutine variables
    If(mpiProcessID.eq.0)Then
      print *,""
      print *,""
      print *,"                           Run Time"
      print *,"----------------------------------------------------------------------"
      print *,"Globals Init Time:         ",globInitTime
      print *,"EAM Potential Load Time:   ",eamLoadTime
      print *,"Configuration input:       ",configLoadTime
      print *,""
      print *,"Neighbour list:            ",nlTime
      print *,"Neighbour list BP:         ",nlTimeBP
      print *,""
      print *,"EFS Calculations:          ",efsCalcTime
      print *,"BP Calculation:            ",evalTimeBP
      print *,""
      print *,"Total time:                ",ProgramTime()
    End If
  End Subroutine runFinaliseEval

End Module finalise
