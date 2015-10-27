Module bpCalcEAM
! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Subroutines and functions to calculate energy, forces and stresses
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
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: calcEnergiesBP
  Public :: calcEnergyBP
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergiesBP()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: energyTimeStart, energyTimeEnd
! Integer(kind=StandardInteger) :: processesPerEnergy
! Start time
    Call cpu_time(energyTimeStart)
! Init variables
    Call resetE()  ! Reset energy
! MPI Processes per energy calculation
! Loop through configurations
    Do configID=1,configCountBP
      If(processMapBP(configID).eq.mpiProcessID)Then  ! build in more mpi options here
        Call calcEnergyBP(configID, configEnergy)
        configCalcEnergiesBP(configID)=&
        (configEnergy/(1.0D0*configurationCoordsKeyBP(configID,2)))
      End If
    End Do
    Call M_synchProcesses()
! Distribute energy array
    Call M_collDouble1D(configCalcEnergiesBP)
    Call M_distDouble1D(configCalcEnergiesBP)
! End time
    Call cpu_time(energyTimeEnd)
! Add time
    efsCalcTimeBP = efsCalcTimeBP + (energyTimeEnd-energyTimeStart)
  End Subroutine calcEnergiesBP
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergyBP(configID, configEnergy)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: configEnergy
    If(eamType.eq.1)Then ! EAM
      Call calcEnergy_EAM(configID, configEnergy)
    End If
    If(eamType.eq.2)Then
      Call calcEnergy_TBEAM(configID, configEnergy)
    End If
  End Subroutine calcEnergyBP
! ---------------------------------------------------------------------------------------------------
! Standard EAM Calculation
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergy_EAM(configID, configEnergy)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: aType, bType, aID_R, bID_R, aID_A, bID_A
    Integer(kind=StandardInteger) :: configID,n,i
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff
    Real(kind=DoubleReal) :: rIJ
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! -----------------
! Init
! -----------------
! Init variables
    rCutoff = bpCutoff
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    configEnergy = 0.0D0
    nlStart = neighbourListKeyBP(configID,1)
    nlLength = neighbourListKeyBP(configID,2)
    nlEnd = neighbourListKeyBP(configID,3)
    configStart = configurationCoordsKeyBP(configID,1)
    configLength = configurationCoordsKeyBP(configID,2)
    configEnd = configurationCoordsKeyBP(configID,3)
! Init arrays
    Do i=1,configLength
      calculationDensity(i) = 0.0D0     ! density at each atom, length = number of atoms in configuration
    End Do
! --------------------------------------------------
! Loop 1 - sum pair energy, pair force and electron density
! --------------------------------------------------
    Do n=nlStart,nlEnd
      aType = neighbourListIBP(n,1)
      bType = neighbourListIBP(n,2)
      aID_R = neighbourListIBP(n,3)   ! ID Relatice (pos in this config)
      bID_R = neighbourListIBP(n,4)
      aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
      bID_A = bID_R+configStart-1
      rIJ = neighbourListRBP(n)
! Check if in cutoff radius
      If(neighbourListRBP(n).le.rCutoff)Then
! Pair potential v(r) and v'(r)
        yArray = SearchPotentialPoint(aType,bType,1,eamType,elementsCount,rIJ)
        pairEnergy = pairEnergy + yArray(1)
! Electron density each A is embedded in due to the electrons of the Bs around it
        yArray = SearchPotentialPoint(bType,0,2,eamType,elementsCount,rIJ)
        calculationDensity(aID_R) = calculationDensity(aID_R) + yArray(1)
! Electron density each B is embedded in due to the electrons of the As around it
        yArray = SearchPotentialPoint(aType,0,2,eamType,elementsCount,rIJ)
        calculationDensity(bID_R) = calculationDensity(bID_R) + yArray(1)
      End If
    End Do
! --------------------------------------------------
! Loop 2 - embedding energy
! --------------------------------------------------
    Do i=1,configLength
      aID_R = i                     ! ID Relatice (pos in this config)
      aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
      aType = configurationCoordsIBP(aID_A,1)
      yArray = SearchPotentialPoint(aType,0,3,eamType,elementsCount,calculationDensity(aID_R))
      embeddingEnergy = embeddingEnergy + yArray(1)
    End Do
! Sum energies
    configEnergy = pairEnergy + embeddingEnergy
  End Subroutine calcEnergy_EAM
! ---------------------------------------------------------------------------------------------------
! Standard Two Band EAM Calculation
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergy_TBEAM(configID, configEnergy)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: aType, bType, aID_R, bID_R, aID_A, bID_A
    Integer(kind=StandardInteger) :: configID,n,i
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff, rIJ
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! -----------------
! Init
! -----------------
! Init variables
    rCutoff = bpCutoff
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    configEnergy = 0.0D0
    nlStart = neighbourListKeyBP(configID,1)
    nlLength = neighbourListKeyBP(configID,2)
    nlEnd = neighbourListKeyBP(configID,3)
    configStart = configurationCoordsKeyBP(configID,1)
    configLength = configurationCoordsKeyBP(configID,2)
    configEnd = configurationCoordsKeyBP(configID,3)
! Init arrays
    Do i=1,configLength
      calculationDensity(i) = 0.0D0     ! density at each atom, length = number of atoms in configuration
      calculationDensityS(i) = 0.0D0
    End Do 
! --------------------------------------------------
! Loop 1 - sum pair energy, pair force and electron density
! --------------------------------------------------
    Do n=nlStart,nlEnd
      aType = neighbourListIBP(n,1)
      bType = neighbourListIBP(n,2)
      aID_R = neighbourListIBP(n,3)   ! ID Relatice (pos in this config)
      bID_R = neighbourListIBP(n,4)
      aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
      bID_A = bID_R+configStart-1
      rIJ = neighbourListRBP(n)
! Check if in cutoff radius
      If(neighbourListRBP(n).le.rCutoff)Then
! Pair potential v(r) and v'(r)
        yArray = SearchPotentialPoint(aType,bType,1,eamType,elementsCount,rIJ)   ! 1: Pair
        pairEnergy = pairEnergy + yArray(1)
        configAtomEnergy(aID_A,1) = configAtomEnergy(aID_A,1) + 0.5D0 * yArray(1)
        configAtomEnergy(bID_A,1) = configAtomEnergy(bID_A,1) + 0.5D0 * yArray(1)
      End If
    End Do
! --------------------------------------------------
! Loop 2 - embedding energy
! --------------------------------------------------
    Do i=1,configLength
      aID_R = i                     ! ID Relatice (pos in this config)
      aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
      aType = configurationCoordsIG(aID_A,1)
! D-band
      yArray = SearchPotentialPoint(aType,0,6,eamType,elementsCount,calculationDensity(aID_R))  ! 6: DEMB
      embeddingEnergy = embeddingEnergy + yArray(1)
      configAtomEnergy(aID_A,2) = configAtomEnergy(aID_A,2) + yArray(1)
! S-band
      yArray = SearchPotentialPoint(aType,0,7,eamType,elementsCount,calculationDensityS(aID_R))  ! 7: SEMB
      embeddingEnergy = embeddingEnergy + yArray(1)
      configAtomEnergy(aID_A,2) = configAtomEnergy(aID_A,2) + yArray(1)
    End Do
! Sum energies
    configEnergy = pairEnergy + embeddingEnergy
  End Subroutine calcEnergy_TBEAM

  Subroutine resetE()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: i
! Init arrays
! Energy
    Do i=1,configCountBP
      configCalcEnergiesBP(i) = -2.1D20
    End Do
  End Subroutine resetE

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!
  Function SearchPotentialPoint (atomA, atomB, potType, eamTypeF, elementsCountF, x) RESULT (yArray)
! Get value of function at x
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: atomA, atomB, potType, eamTypeF, elementsCountF, funcKey
    Integer(kind=StandardInteger) :: funcStart, funcLength, funcEnd
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Function key
    funcKey = FunctionKey (atomA, atomB, potType, eamTypeF, elementsCountF)
! Function start/end
    funcStart = eamKey(funcKey,4)
    funcLength = eamKey(funcKey,5)
    funcEnd = eamKey(funcKey,6)
! Get values
    yArray = PointInterp(eamData,x,eamInterpPoints,1,funcStart,funcLength)
  End Function SearchPotentialPoint
! ---------------------------------------------------------------------------------------------------

End Module bpCalcEAM
