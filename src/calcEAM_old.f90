Module calcEAM

! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: calcEAM
! Updated: 18th May 2015
! --------------------------------------------------------------!
! Description:
! Loop through all configurations
! Build neighbour list for each configuration
! --------------------------------------------------------------!

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
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: calcEnergies
  Public :: calcEnergy
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergies()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID, forceKeyStart, forceKeyEnd, selectedProcess
    Real(kind=DoubleReal) :: configEnergy
! Integer(kind=StandardInteger) :: processesPerEnergy
! Start time
    Call cpu_time(timeStart)
! Init variables
    configCalcEnergies = -2.1D20
    configCalcForces = -2.1D20
    configCalcStresses = -2.1D20
    configEnergy = 0.0D0
! MPI Processes per energy calculation
! Loop through configurations
    Do configID=1,configCount
      If(processMap(configID).eq.mpiProcessID)Then
        print *,"Calc energy ",configID,mpiProcessID
        Call calcEnergy(configID, configEnergy)
        print *,"Energy ",configEnergy
      End If
! If(processMap(configID).eq.mod(mpiProcessID,configCount))Then  ! build in more mpi options here
!  print *,"Calc energy ",configID,mpiProcessID
!  Call calcEnergy(configID, configEnergy, 1)
!  If(processMap(configID).eq.mpiProcessID)Then
! configCalcEnergies(configID) = configEnergy  !Store result from primary mpi process
!  End If
! End If
    End Do
! Distribute energy array
    Call M_collDouble1D(configCalcEnergies)
    Call M_distDouble1D(configCalcEnergies)
! Distribute force array
    Do configID=1,configCount
      selectedProcess = processMap(configID)
      forceKeyStart = configurationCoordsKeyG(configID,1)
      forceKeyEnd = configurationCoordsKeyG(configID,3)
      Call M_collDouble2D(configCalcForces, selectedProcess, forceKeyStart, forceKeyEnd)
    End Do
    Call M_distDouble2D(configCalcForces)
! Distribute Stress Array
    Call M_collDouble2D(configCalcStresses)
    Call M_distDouble2D(configCalcStresses)
! Output forces to file
    If(saveForcesToFile.eq.1)Then
      Call outputForcesFile()
    End If
! End time
    Call cpu_time(timeEnd)
    Call timeAcc(efsCalcTime,timeStart,timeEnd)
  End Subroutine calcEnergies
! ---------------------------------------------------------------------------------------------------
!  Subroutine calcEnergy(configID, configEnergy, forceCalcIn, pairEnergyOut, embeddingEnergyOut)
  Subroutine calcEnergy(configID, totalEnergy)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, n, nKey, i, j, k
    Integer(kind=StandardInteger) :: aType, bType, aID, bID
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, totalEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Init variables
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    totalEnergy = 0.0D0
! NL Details
    nlStart = neighbourListKey(configID,1)
    nlLength = neighbourListKey(configID,2)
    nlEnd = neighbourListKey(configID,3)
! Config details
    configStart = configurationCoordsKeyG(configID,1)
    configLength = configurationCoordsKeyG(configID,2)
    configEnd = configurationCoordsKeyG(configID,3)
! Init calculation density array for required length
    Do i=1,configLength
      calculationDensity(i) = 0.0D0
    End Do
! Init energies
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    totalEnergy = 0.0D0
! --------------------------------------------------
! Loop 1 - sum pair energy and density
! --------------------------------------------------
    print *,nlStart,nlEnd
    Do n=nlStart,nlEnd
      aType = neighbourListI(n,1)
      bType = neighbourListI(n,2)
      aID = neighbourListI(n,3)
! Pair potential
      yArray = SearchPotentialPoint(aType,bType,1,eamType,neighbourListR(n))
      pairEnergy = pairEnergy + yArray(1)
! Electron density of each A due to the electrons of the Bs around it
      yArray = SearchPotentialPoint(bType,0,2,eamType,neighbourListR(n))
      calculationDensity(aID) = calculationDensity(aID) + yArray(1)
    End Do
! --------------------------------------------------
! Loop 2 - embedding energy
! --------------------------------------------------
    Do aID=1,configLength
      aType = configurationCoordsIG(aID,1)
      yArray = SearchPotentialPoint(aType,0,3,eamType,calculationDensity(aID))
      embeddingEnergy = embeddingEnergy + yArray(1)
    End Do
    print *,pairEnergy,embeddingEnergy,(pairEnergy+embeddingEnergy),&
    ((pairEnergy+embeddingEnergy)/configLength)
  End Subroutine calcEnergy
! ---------------------------------------------------------------------------------------------------
! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!
  Function SearchPotentialPoint (atomA, atomB, potType, eamTypeF, x) RESULT (yArray)
! Get value of function at x
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: atomA, atomB, potType, eamTypeF, funcKey
    Integer(kind=StandardInteger) :: funcStart, funcLength, funcEnd
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Function key
    funcKey = FunctionKey (atomA, atomB, potType, eamTypeF)
! Function start/end
    funcStart = eamKey(funcKey,4)
    funcLength = eamKey(funcKey,5)
    funcEnd = eamKey(funcKey,6)
! Get values
    yArray = PointInterp(eamData,x,eamInterpPoints,1,funcStart,funcLength)
  End Function SearchPotentialPoint
! ---------------------------------------------------------------------------------------------------
  Function FunctionKey (atomA, atomB, potType, eamTypeF) RESULT (funcKey)
! Returns the potential key
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: atomA, atomB, potType, eamTypeF, funcKey
    Integer(kind=StandardInteger) :: atomMin, atomMax
! 1 "PAIR"
! 2 "DENS"
! 3 "EMBE"
! 4 "DDEN"
! 5 "SDEN"
! 6 "DEMB"
! 7 "SEMB"
    If(eamTypeF.eq.1)Then
      If(potType.eq.1)Then          ! PAIR
        atomMax = max(atomA,atomB)-1
        atomMin = min(atomA,atomB)-1
        funcKey = 1+atomMin+(atomMax*(atomMax+1))/2
      End If
      If(potType.eq.2)Then          ! DENS
        funcKey = eamPairCount + atomA
      End If
      If(potType.eq.3)Then          ! EMBE
        funcKey = eamPairCount + eamDensCount + atomA
      End If
    End If
    If(eamTypeF.eq.2)Then
      If(potType.eq.1)Then          ! PAIR
        atomMax = max(atomA,atomB)-1
        atomMin = min(atomA,atomB)-1
        funcKey = 1+atomMin+(atomMax*(atomMax+1))/2
      End If
      If(potType.eq.5)Then          ! SDEN
        If(elementsCount.eq.1)Then
          funcKey = eamPairCount + 1
        Else
          atomMax = max(atomA,atomB)-1
          atomMin = min(atomA,atomB)-1
          funcKey = eamPairCount+1+atomMin+(atomMax*(atomMax-1))/2
        End If
      End If
      If(potType.eq.4)Then          ! DDEN
        funcKey = eamPairCount + eamSdenCount + atomA
      End If
      If(potType.eq.7)Then          ! SEMB
        funcKey = eamPairCount + eamDensCount + eamSdenCount + eamDembCount + atomA
      End If
      If(potType.eq.6)Then          ! DEMB
        funcKey = eamPairCount + eamDensCount + eamSdenCount + atomA
      End If
    End If
  End Function FunctionKey
End Module calcEAM
