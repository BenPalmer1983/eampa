Module opti
! --------------------------------------------------------------!
! Optimise EAM functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calls the eval and calcEAM subroutines to optimise potential functions
! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------
! Setup Modules
  Use kinds
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use globals
  Use initialise
  Use loadData
  Use output
  Use eamGen
  Use calcEAM
  Use eval
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: optiEAM
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine optiEAM()
    Implicit None   ! Force declaration of all variables
    !Integer(kind=StandardInteger) :: configID, keyStart, keyEnd, selectedProcess
    !Real(kind=DoubleReal) :: configEnergy

    
  End Subroutine optiEAM

  
  
  
! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!


End Module opti
