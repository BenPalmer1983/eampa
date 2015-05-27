Module prep

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

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
  Use initialise
  Use loadData
  Use globals
  Use output
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: runPrep

  Contains
  Subroutine runPrep()
    Implicit None   ! Force declaration of all variables
! Private variables
    Call setProcessMap()
  End Subroutine runPrep

  Subroutine setProcessMap()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j
! Init variables
    i = 0
    j = 0
! Energy/force/stress calculations
    Do i=1,configCount
      processMap(i,1) = mod(i-1,mpiProcessCount)
    End Do
! Equilibrium volume calculations
    If(calcEqVol(1:3).eq."ALL")Then
      Do i=1,configCount
        processMap(i,2) = mod(i-1,mpiProcessCount)
      End Do
    ElseIf(calcEqVol(1:3).eq."SEL")Then
      Do i=1,configCount
        If(configRefEV(i).gt.-2.0D20)Then
          j = j + 1
          processMap(i,2) = mod(j-1,mpiProcessCount)
        End If
      End Do
    Else !none
! do nothing
    End If
! Bulk modulus calculations
    Do i=1,configCount
      If(configRefBM(i).gt.-2.0D20)Then
        j = j + 1
        processMap(i,3) = mod(j-1,mpiProcessCount)
      End If
    End Do
! Output to file
    Call outputProcessMap()
  End Subroutine setProcessMap

End Module prep

