Module eamGen
! --------------------------------------------------------------!
! General EAM subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Subroutines and functions used through
! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Public Functions
  Public :: FunctionKey
  Public :: FunctionCount
  Contains
! ---------------------------------------------------------------------------------------------------
  Function FunctionKey (atomA, atomB, potType, eamTypeF, elementsCountF) RESULT (funcKey)
! Returns the potential key
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: atomA, atomB, potType, eamTypeF, elementsCountF, funcKey
    Integer(kind=StandardInteger) :: atomMin, atomMax
    Integer(kind=StandardInteger) :: nE, nF
    Integer(kind=StandardInteger) :: nPair, nDens, nEmbe
    Integer(kind=StandardInteger) :: nDden, nSden, nDemb, nSemb
! potType
! 1 "PAIR"
! 2 "DENS"
! 3 "EMBE"
! 4 "DDEN"
! 5 "SDEN"
! 6 "DEMB"
! 7 "SEMB"
! eamTypeF
! 1 EAM
! 2 TBEAM
! 3 CDEAM
    funcKey = 0
    nE = elementsCountF
! ------------------------
! EAM Function Key
! ------------------------
    If(eamTypeF.eq.1)Then
      nF = FunctionCount(0,1,nE)    !FunctionCount (potType, eamTypeF, elementsCountF)
      nPair = FunctionCount(1,1,nE)
      nDens = FunctionCount(2,1,nE)
      nEmbe = FunctionCount(3,1,nE)
      If(potType.eq.1)Then          ! PAIR
        atomMax = max(atomA,atomB)-1
        atomMin = min(atomA,atomB)-1
        funcKey = 1+atomMin+(atomMax*(atomMax+1))/2
      End If
      If(potType.eq.2)Then          ! DENS
        funcKey = nPair + atomA
      End If
      If(potType.eq.3)Then          ! EMBE
        funcKey = nPair + nDens + atomA
      End If
    End If
    If(eamTypeF.eq.2)Then
! ------------------------
! 2BEAM/TBEAM Function Key
! ------------------------
      If(elementsCountF.eq.1)Then      ! One element, separate D-band and S-band
        nF = FunctionCount(0,2,nE)    !FunctionCount (potType, eamTypeF, elementsCountF)
        nPair = FunctionCount(1,2,nE)
        nSden = FunctionCount(5,2,nE)
        nDden = FunctionCount(4,2,nE)
        nSemb = FunctionCount(7,2,nE)
        nDemb = FunctionCount(6,2,nE)
        If(potType.eq.1)Then          ! PAIR
          atomMax = max(atomA,atomB)-1
          atomMin = min(atomA,atomB)-1
          funcKey = 1+atomMin+(atomMax*(atomMax+1))/2
        End If
        If(potType.eq.5)Then          ! SDEN
          funcKey = nPair + 1
        End If
        If(potType.eq.4)Then          ! DDEN
          funcKey = nPair + nSden + atomA
        End If
        If(potType.eq.7)Then          ! SEMB
          funcKey = nPair + nSden + nDden + atomA
        End If
        If(potType.eq.6)Then          ! DEMB
          funcKey = nPair + nSden + nDden + nSemb + atomA
        End If
      Else
        nF = FunctionCount(0,2,nE)    !FunctionCount (potType, eamTypeF, elementsCountF)
        nPair = FunctionCount(1,2,nE)
        nSden = FunctionCount(5,2,nE)
        nDden = FunctionCount(4,2,nE)
        nSemb = FunctionCount(7,2,nE)
        nDemb = FunctionCount(6,2,nE)
        If(potType.eq.1)Then          ! PAIR
          atomMax = max(atomA,atomB)-1
          atomMin = min(atomA,atomB)-1
          funcKey = 1+atomMin+(atomMax*(atomMax+1))/2
        End If
        If(potType.eq.5)Then          ! SDEN
          atomMax = max(atomA,atomB)-1
          atomMin = min(atomA,atomB)-1
          funcKey = 1+nPair+atomMin+(atomMax*(atomMax-1))/2
        End If
        If(potType.eq.4)Then          ! DDEN
          funcKey = nPair + nSden + atomA
        End If
        If(potType.eq.7)Then          ! SEMB
          funcKey = nPair + nSden + nDden + atomA
        End If
        If(potType.eq.6)Then          ! DEMB
          funcKey = nPair + nSden + nDden + nSemb + atomA
        End If
      End If
    End If
  End Function FunctionKey
! ---------------------------------------------------------------------------------------------------
  Function FunctionCount (potType, eamTypeF, elementsCountF) RESULT (funcCount)
! Returns the potential key
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: potType, eamTypeF, elementsCountF, funcCount
    Integer(kind=StandardInteger) :: nE
! potType
! 0 "ALL"
! 1 "PAIR"
! 2 "DENS"
! 3 "EMBE"
! 4 "DDEN"
! 5 "SDEN"
! 6 "DEMB"
! 7 "SEMB"
! eamTypeF
! 1 EAM
! 2 TBEAM
! 3 CDEAM
    funcCount = 0
    nE = elementsCountF
! ------------------------
! EAM Function Key
! ------------------------
    If(eamTypeF.eq.1)Then
      If(potType.eq.0)Then          ! All functions
        funcCount = (nE*(nE+5))/2
      End If
      If(potType.eq.1)Then          ! PAIR
        funcCount = (nE*(nE+1))/2
      End If
      If(potType.eq.2)Then          ! DENS
        funcCount = nE
      End If
      If(potType.eq.3)Then          ! EMBE
        funcCount = nE
      End If
    End If
    If(eamTypeF.eq.2)Then
! ------------------------
! 2BEAM/TBEAM Function Key
! ------------------------
      If(elementsCountF.eq.1)Then      ! One element, separate D-band and S-band
        If(potType.eq.0)Then          ! All functions
          funcCount = (nE*(nE+7))/2
        End If
        If(potType.eq.1)Then          ! PAIR
          funcCount = (nE*(nE+1))/2
        End If
        If(potType.eq.5)Then          ! SDEN
          funcCount = 1
        End If
        If(potType.eq.4)Then          ! DDEN
          funcCount = nE
        End If
        If(potType.eq.7)Then          ! SEMB
          funcCount = nE
        End If
        If(potType.eq.6)Then          ! DEMB
          funcCount = nE
        End If
      Else
        If(potType.eq.0)Then          ! All functions
          funcCount = (nE*(nE+3))
        End If
        If(potType.eq.1)Then          ! PAIR
          funcCount = (nE*(nE+1))/2
        End If
        If(potType.eq.5)Then          ! SDEN
          funcCount = (nE*(nE-1))/2
        End If
        If(potType.eq.4)Then          ! DDEN
          funcCount = nE
        End If
        If(potType.eq.7)Then          ! SEMB
          funcCount = nE
        End If
        If(potType.eq.6)Then          ! DEMB
          funcCount = nE
        End If
      End If
    End If
  End Function FunctionCount
! ---------------------------------------------------------------------------------------------------
End Module eamGen
