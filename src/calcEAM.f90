Module calcEAM
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
  Public :: calcEnergies
  Public :: calcEnergy
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergies()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID, keyStart, keyEnd, selectedProcess
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: energyTimeStart, energyTimeEnd
! Integer(kind=StandardInteger) :: processesPerEnergy
! Start time
    Call cpu_time(energyTimeStart)
! Init variables
    Call resetESF()  ! Reset energy, stress and force arrays
! MPI Processes per energy calculation
! Loop through configurations
    forceStressSwitch = 1
    Do configID=1,configCount
      If(processMap(configID).eq.mpiProcessID)Then  ! build in more mpi options here
        Call calcEnergy(configID, configEnergy)
        configCalcEnergies(configID)=&
        (configEnergy/(1.0D0*configurationCoordsKeyG(configID,2)))
      End If
    End Do
    Call M_synchProcesses()
! Distribute energy array
    Call M_collDouble1D(configCalcEnergies)
    Call M_distDouble1D(configCalcEnergies)
! Collect and distribute forces array
    Call M_synchProcesses()
    Do configID=1,configCount
      selectedProcess = processMap(configID)
      keyStart = configurationCoordsKeyG(configID,1)
      keyEnd = configurationCoordsKeyG(configID,3)
      Call M_collDouble2D(configCalcForces, selectedProcess, keyStart, keyEnd)
    End Do
    Call M_distDouble2D(configCalcForces)
! Collect and distribute atom energy array
    Call M_synchProcesses()
    Do configID=1,configCount
      selectedProcess = processMap(configID)
      keyStart = configurationCoordsKeyG(configID,1)
      keyEnd = configurationCoordsKeyG(configID,3)
      Call M_collDouble2D(configAtomEnergy, selectedProcess, keyStart, keyEnd)
    End Do
    Call M_distDouble2D(configAtomEnergy)
! Collect and distribute stress array
    Call M_synchProcesses()
    Do configID=1,configCount
      selectedProcess = processMap(configID)
      keyStart = configurationCoordsKeyG(configID,1)
      keyEnd = configurationCoordsKeyG(configID,3)
      Call M_collDouble2D(configCalcStresses)
    End Do
    Call M_distDouble2D(configCalcStresses)
! Output forces to file
    If(saveForcesToFile)Then
      Call outputForcesFile()
    End If
    Call outputAtomEnergiesFile()
    Call outputEnergyT()
! End time
    Call cpu_time(energyTimeEnd)
! Add time
    efsCalcTime = efsCalcTime + (energyTimeEnd-energyTimeStart)
  End Subroutine calcEnergies
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergy(configID, configEnergy)
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: configEnergy
    If(eamType.eq.1)Then ! EAM
      Call calcEnergy_EAM(configID, configEnergy)
    End If
    If(eamType.eq.2)Then
      Call calcEnergy_TBEAM(configID, configEnergy)
    End If
  End Subroutine calcEnergy
! ---------------------------------------------------------------------------------------------------
! Standard EAM Calculation
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergy_EAM(configID, configEnergy)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: aType, bType, aID_R, bID_R, aID_A, bID_A
    Integer(kind=StandardInteger) :: configID,n,i,j,k
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff
    Real(kind=DoubleReal) :: rIJ, fX, fY, fZ
    Real(kind=DoubleReal) :: stressNX, stressNY, stressNZ
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: densDerivAB, densDerivBA, embeDerivA, embeDerivB, forceM
! -----------------
! Init
! -----------------
! Init variables
    rCutoff = configurationsR(configID,11)
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    configEnergy = 0.0D0
    nlStart = neighbourListKey(configID,1)
    nlLength = neighbourListKey(configID,2)
    nlEnd = neighbourListKey(configID,3)
    configStart = configurationCoordsKeyG(configID,1)
    configLength = configurationCoordsKeyG(configID,2)
    configEnd = configurationCoordsKeyG(configID,3)
! Init arrays
    Do i=1,configLength
      calculationDensity(i) = 0.0D0     ! density at each atom, length = number of atoms in configuration
    End Do
    If(forceStressSwitch.eq.1)Then
! Do i=1,nlLength
!  pairForce(i) = 0.0D0            ! Pair forces, length of neighbour list
! End Do
! Do i=configStart,configEnd        ! store directly into global array
! configCalcForces(i,1) = 0.0D0   ! Force x component atom i
! configCalcForces(i,2) = 0.0D0   ! Force y component atom i
! configCalcForces(i,3) = 0.0D0   ! Force z component atom i
! configAtomEnergy(i,1) = 0.0D0   ! Atom energy from pair
! configAtomEnergy(i,2) = 0.0D0   ! Atom energy from embedding
! End Do
! Init stress array for this config
! Do i=1,9                                  ! store directly into global array
! configCalcStresses(configID,i) = 0.0D0  ! xx, xy...zy, zz for config i
! End Do
    End If
!     neighbourListI(n,1)  Atom A type
!     neighbourListI(n,2)  Atom B type
!     neighbourListI(n,3)  Atom A id
!     neighbourListI(n,4)  Atom B id
! --------------------------------------------------
! Loop 1 - sum pair energy, pair force and electron density
! --------------------------------------------------
    Do n=nlStart,nlEnd
      aType = neighbourListI(n,1)
      bType = neighbourListI(n,2)
      aID_R = neighbourListI(n,3)   ! ID Relatice (pos in this config)
      bID_R = neighbourListI(n,4)
      aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
      bID_A = bID_R+configStart-1
      rIJ = neighbourListR(n)
! Check if in cutoff radius
      If(neighbourListR(n).le.rCutoff)Then
! Pair potential v(r) and v'(r)
        yArray = SearchPotentialPoint(aType,bType,1,eamType,elementsCount,rIJ)
        pairEnergy = pairEnergy + yArray(1)
        configAtomEnergy(aID_A,1) = configAtomEnergy(aID_A,1) + 0.5D0 * yArray(1)
        configAtomEnergy(bID_A,1) = configAtomEnergy(bID_A,1) + 0.5D0 * yArray(1)
! If force-stress: store force from Pair Potential
        If(forceStressSwitch.eq.1)Then
          fX = -1.0D0*neighbourListCoords(n,10)*yArray(2)
          fY = -1.0D0*neighbourListCoords(n,11)*yArray(2)
          fZ = -1.0D0*neighbourListCoords(n,12)*yArray(2)
! Pair force on atom A
          configCalcForces(aID_A,1) = configCalcForces(aID_A,1) + fX
          configCalcForces(aID_A,2) = configCalcForces(aID_A,2) + fY
          configCalcForces(aID_A,3) = configCalcForces(aID_A,3) + fZ
! Pair force on atom B
          configCalcForces(bID_A,1) = configCalcForces(bID_A,1) - fX
          configCalcForces(bID_A,2) = configCalcForces(bID_A,2) - fY
          configCalcForces(bID_A,3) = configCalcForces(bID_A,3) - fZ
! Virial Stress pair force
          Do i=1,3  ! Loop through coordinate axis x,y,z
            j=1     ! force x
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNX = neighbourListCoords(n,6+i)*fX
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX  ! Stress doubles as -F and -Rij negatives cancel
            End If
            j=2     ! force y
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNY = neighbourListCoords(n,6+i)*fY
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            End If
            j=3     ! force z
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNZ = neighbourListCoords(n,6+i)*fZ
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            End If
! Atom i (Atom A) are always in the volume
          End Do
        End If
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
      aType = configurationCoordsIG(aID_A,1)
      yArray = SearchPotentialPoint(aType,0,3,eamType,elementsCount,calculationDensity(aID_R))
      embeddingEnergy = embeddingEnergy + yArray(1)
      configAtomEnergy(aID_A,2) = configAtomEnergy(aID_A,2) + yArray(1)
    End Do
! Sum energies
    configEnergy = pairEnergy + embeddingEnergy
! --------------------------------------------------
! Loop 3 - remaining derivative values for force
! --------------------------------------------------
    If(forceStressSwitch.eq.1)Then
      Do n=nlStart,nlEnd
        aType = neighbourListI(n,1)
        bType = neighbourListI(n,2)
        aID_R = neighbourListI(n,3)   ! ID Relatice (pos in this config)
        bID_R = neighbourListI(n,4)
        aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
        bID_A = bID_R+configStart-1
        rIJ = neighbourListR(n)
! Check if in cutoff radius
        If(neighbourListR(n).le.rCutoff)Then
! @Fi(p)/@p
          yArray = SearchPotentialPoint(aType,0,3,eamType,elementsCount,calculationDensity(aID_R))
          embeDerivA = yArray(2)
! @Pji(r)/@r
          yArray = SearchPotentialPoint(bType,0,2,eamType,elementsCount,rIJ)
          densDerivBA = yArray(2)
! @Fj(p)/@p
          yArray = SearchPotentialPoint(bType,0,3,eamType,elementsCount,calculationDensity(bID_R))
          embeDerivB = yArray(2)
! @Pij(r)/@r
          yArray = SearchPotentialPoint(aType,0,2,eamType,elementsCount,rIJ)
          densDerivAB = yArray(2)
! Forces
          forceM = -1.0D0*(embeDerivA*densDerivBA+embeDerivB*densDerivAB)
          fX = -1.0D0*neighbourListCoords(n,10)*forceM
          fY = -1.0D0*neighbourListCoords(n,11)*forceM
          fZ = -1.0D0*neighbourListCoords(n,12)*forceM
! Pair force on atom A
          configCalcForces(aID_A,1) = configCalcForces(aID_A,1) + fX
          configCalcForces(aID_A,2) = configCalcForces(aID_A,2) + fY
          configCalcForces(aID_A,3) = configCalcForces(aID_A,3) + fZ
! Pair force on atom B
          configCalcForces(bID_A,1) = configCalcForces(bID_A,1) - fX
          configCalcForces(bID_A,2) = configCalcForces(bID_A,2) - fY
          configCalcForces(bID_A,3) = configCalcForces(bID_A,3) - fZ
! Virial Stress pair force
          Do i=1,3  ! Loop through coordinate axis x,y,z
            j=1     ! force x
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNX = neighbourListCoords(n,6+i)*fX
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX  ! Stress doubles as -F and -Rij negatives cancel
            End If
            j=2     ! force y
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNY = neighbourListCoords(n,6+i)*fY
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            End If
            j=3     ! force z
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNZ = neighbourListCoords(n,6+i)*fZ
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            End If
! Atom i (Atom A) are always in the volume
          End Do
        End If
      End Do
! Multiply stress values by factor of volume of domain
      Do i=1,3
        Do j=1,3
          k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ...
          configCalcStresses(configID,k) = &
          (0.5D0/configVolume(configID))*configCalcStresses(configID,k)
          configCalcStresses(configID,k) = &
          UnitConvert(configCalcStresses(configID,k),"EVAN3","GPA")
        End Do
      End Do
    End If
  End Subroutine calcEnergy_EAM
! ---------------------------------------------------------------------------------------------------
! Standard Two Band EAM Calculation
! ---------------------------------------------------------------------------------------------------
  Subroutine calcEnergy_TBEAM(configID, configEnergy)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: aType, bType, aID_R, bID_R, aID_A, bID_A
    Integer(kind=StandardInteger) :: configID,n,i,j,k
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff
    Real(kind=DoubleReal) :: rIJ, fX, fY, fZ
    Real(kind=DoubleReal) :: stressNX, stressNY, stressNZ
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal) :: densDerivAB, densDerivBA, embeDerivA, embeDerivB, forceM
! -----------------
! Init
! -----------------
! FunctionKey (atomA, atomB, potType, eamTypeF, elementsCountF)
! print *,FunctionKey (1, 1, potType, 2)
!    print *,"1",FunctionKey (1, 1, 1, 2, 1)
!    print *,"4",FunctionKey (1, 1, 4, 2, 1)
!    print *,"5",FunctionKey (1, 1, 5, 2, 1)
!    print *,"6",FunctionKey (1, 1, 6, 2, 1)
!    print *,"7",FunctionKey (1, 1, 7, 2, 1)
! Init variables
    rCutoff = configurationsR(configID,11)
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    configEnergy = 0.0D0
    nlStart = neighbourListKey(configID,1)
    nlLength = neighbourListKey(configID,2)
    nlEnd = neighbourListKey(configID,3)
    configStart = configurationCoordsKeyG(configID,1)
    configLength = configurationCoordsKeyG(configID,2)
    configEnd = configurationCoordsKeyG(configID,3)
! Init arrays
    Do i=1,configLength
      calculationDensity(i) = 0.0D0     ! density at each atom, length = number of atoms in configuration
      calculationDensityS(i) = 0.0D0
    End Do
    If(forceStressSwitch.eq.1)Then
! Do i=1,nlLength
!  pairForce(i) = 0.0D0            ! Pair forces, length of neighbour list
! End Do
! Do i=configStart,configEnd        ! store directly into global array
!  configCalcForces(i,1) = 0.0D0   ! Force x component atom i
!  configCalcForces(i,2) = 0.0D0   ! Force y component atom i
!  configCalcForces(i,3) = 0.0D0   ! Force z component atom i
!  configAtomEnergy(i,1) = 0.0D0   ! Atom energy from pair
!  configAtomEnergy(i,2) = 0.0D0   ! Atom energy from embedding
! End Do
! Init stress array for this config
! Do i=1,9                                  ! store directly into global array
!  configCalcStresses(configID,i) = 0.0D0  ! xx, xy...zy, zz for config i
! End Do
    End If
!     neighbourListI(n,1)  Atom A type
!     neighbourListI(n,2)  Atom B type
!     neighbourListI(n,3)  Atom A id
!     neighbourListI(n,4)  Atom B id
! --------------------------------------------------
! Loop 1 - sum pair energy, pair force and electron density
! --------------------------------------------------
    Do n=nlStart,nlEnd
      aType = neighbourListI(n,1)
      bType = neighbourListI(n,2)
      aID_R = neighbourListI(n,3)   ! ID Relatice (pos in this config)
      bID_R = neighbourListI(n,4)
      aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
      bID_A = bID_R+configStart-1
      rIJ = neighbourListR(n)
! Check if in cutoff radius
      If(neighbourListR(n).le.rCutoff)Then
! Pair potential v(r) and v'(r)
        yArray = SearchPotentialPoint(aType,bType,1,eamType,elementsCount,rIJ)   ! 1: Pair
! delete If(n.le.20)Then
! delete print *,aType,bType,1,eamType,elementsCount,rIJ,yArray(1)
! delete End If
        pairEnergy = pairEnergy + yArray(1)
        configAtomEnergy(aID_A,1) = configAtomEnergy(aID_A,1) + 0.5D0 * yArray(1)
        configAtomEnergy(bID_A,1) = configAtomEnergy(bID_A,1) + 0.5D0 * yArray(1)
! If force-stress: store force from Pair Potential
        If(forceStressSwitch.eq.1)Then
          fX = -1.0D0*neighbourListCoords(n,10)*yArray(2)
          fY = -1.0D0*neighbourListCoords(n,11)*yArray(2)
          fZ = -1.0D0*neighbourListCoords(n,12)*yArray(2)
! Pair force on atom A
          configCalcForces(aID_A,1) = configCalcForces(aID_A,1) + fX
          configCalcForces(aID_A,2) = configCalcForces(aID_A,2) + fY
          configCalcForces(aID_A,3) = configCalcForces(aID_A,3) + fZ
! Pair force on atom B
          configCalcForces(bID_A,1) = configCalcForces(bID_A,1) - fX
          configCalcForces(bID_A,2) = configCalcForces(bID_A,2) - fY
          configCalcForces(bID_A,3) = configCalcForces(bID_A,3) - fZ
! Virial Stress pair force
          Do i=1,3  ! Loop through coordinate axis x,y,z
            j=1     ! force x
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNX = neighbourListCoords(n,6+i)*fX
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX  ! Stress doubles as -F and -Rij negatives cancel
            End If
            j=2     ! force y
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNY = neighbourListCoords(n,6+i)*fY
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            End If
            j=3     ! force z
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNZ = neighbourListCoords(n,6+i)*fZ
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            End If
! Atom i (Atom A) are always in the volume
          End Do
        End If
! D-Band 4 "DDEN"
! Electron density each A is embedded in due to the electrons of the Bs around it
        yArray = SearchPotentialPoint(bType,0,4,eamType,elementsCount,rIJ)    ! 4: DDEN
        calculationDensity(aID_R) = calculationDensity(aID_R) + yArray(1)
! Electron density each B is embedded in due to the electrons of the As around it
        yArray = SearchPotentialPoint(aType,0,4,eamType,elementsCount,rIJ)    ! 4: DDEN
        calculationDensity(bID_R) = calculationDensity(bID_R) + yArray(1)
! S-Band 5 "SDEN"
! Electron density each A is embedded in due to the electrons of the Bs around it
        yArray = SearchPotentialPoint(bType,aType,5,eamType,elementsCount,rIJ)  ! 5: SDEN
        calculationDensityS(aID_R) = calculationDensityS(aID_R) + yArray(1)
! Electron density each B is embedded in due to the electrons of the As around it
        yArray = SearchPotentialPoint(aType,bType,5,eamType,elementsCount,rIJ)  ! 5: SDEN
        calculationDensityS(bID_R) = calculationDensityS(bID_R) + yArray(1)
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
! --------------------------------------------------
! Loop 3 - remaining derivative values for force
! --------------------------------------------------
    If(forceStressSwitch.eq.1)Then
      Do n=nlStart,nlEnd
        aType = neighbourListI(n,1)
        bType = neighbourListI(n,2)
        aID_R = neighbourListI(n,3)   ! ID Relatice (pos in this config)
        bID_R = neighbourListI(n,4)
        aID_A = aID_R+configStart-1   ! ID Absolute (pos in in all configs)
        bID_A = bID_R+configStart-1
        rIJ = neighbourListR(n)
! Check if in cutoff radius
        If(neighbourListR(n).le.rCutoff)Then
! D-Band
! @Fi(p)/@p
          yArray = SearchPotentialPoint(aType,0,6,eamType,elementsCount,calculationDensity(aID_R))
          embeDerivA = yArray(2)
! @Pji(r)/@r
          yArray = SearchPotentialPoint(bType,0,4,eamType,elementsCount,rIJ)
          densDerivBA = yArray(2)
! @Fj(p)/@p
          yArray = SearchPotentialPoint(bType,0,6,eamType,elementsCount,calculationDensity(bID_R))
          embeDerivB = yArray(2)
! @Pij(r)/@r
          yArray = SearchPotentialPoint(aType,0,4,eamType,elementsCount,rIJ)
          densDerivAB = yArray(2)
! Forces
          forceM = -1.0D0*(embeDerivA*densDerivBA+embeDerivB*densDerivAB)
          fX = -1.0D0*neighbourListCoords(n,10)*forceM
          fY = -1.0D0*neighbourListCoords(n,11)*forceM
          fZ = -1.0D0*neighbourListCoords(n,12)*forceM
! Pair force on atom A
          configCalcForces(aID_A,1) = configCalcForces(aID_A,1) + fX
          configCalcForces(aID_A,2) = configCalcForces(aID_A,2) + fY
          configCalcForces(aID_A,3) = configCalcForces(aID_A,3) + fZ
! Pair force on atom B
          configCalcForces(bID_A,1) = configCalcForces(bID_A,1) - fX
          configCalcForces(bID_A,2) = configCalcForces(bID_A,2) - fY
          configCalcForces(bID_A,3) = configCalcForces(bID_A,3) - fZ
! Virial Stress pair force
          Do i=1,3  ! Loop through coordinate axis x,y,z
            j=1     ! force x
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNX = neighbourListCoords(n,6+i)*fX
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX  ! Stress doubles as -F and -Rij negatives cancel
            End If
            j=2     ! force y
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNY = neighbourListCoords(n,6+i)*fY
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            End If
            j=3     ! force z
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNZ = neighbourListCoords(n,6+i)*fZ
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            End If
! Atom i (Atom A) are always in the volume
          End Do
! S-Band
! @Fi(p)/@p
          yArray = SearchPotentialPoint(aType,0,7,eamType,elementsCount,calculationDensityS(aID_R))
          embeDerivA = yArray(2)
! @Pji(r)/@r
          yArray = SearchPotentialPoint(bType,0,5,eamType,elementsCount,rIJ)
          densDerivBA = yArray(2)
! @Fj(p)/@p
          yArray = SearchPotentialPoint(bType,0,7,eamType,elementsCount,calculationDensityS(bID_R))
          embeDerivB = yArray(2)
! @Pij(r)/@r
          yArray = SearchPotentialPoint(aType,0,5,eamType,elementsCount,rIJ)
          densDerivAB = yArray(2)
! Forces
          forceM = -1.0D0*(embeDerivA*densDerivBA+embeDerivB*densDerivAB)
          fX = -1.0D0*neighbourListCoords(n,10)*forceM
          fY = -1.0D0*neighbourListCoords(n,11)*forceM
          fZ = -1.0D0*neighbourListCoords(n,12)*forceM
! Pair force on atom A
          configCalcForces(aID_A,1) = configCalcForces(aID_A,1) + fX
          configCalcForces(aID_A,2) = configCalcForces(aID_A,2) + fY
          configCalcForces(aID_A,3) = configCalcForces(aID_A,3) + fZ
! Pair force on atom B
          configCalcForces(bID_A,1) = configCalcForces(bID_A,1) - fX
          configCalcForces(bID_A,2) = configCalcForces(bID_A,2) - fY
          configCalcForces(bID_A,3) = configCalcForces(bID_A,3) - fZ
! Virial Stress pair force
          Do i=1,3  ! Loop through coordinate axis x,y,z
            j=1     ! force x
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNX = neighbourListCoords(n,6+i)*fX
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNX  ! Stress doubles as -F and -Rij negatives cancel
            End If
            j=2     ! force y
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNY = neighbourListCoords(n,6+i)*fY
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNY
            End If
            j=3     ! force z
            k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ... 3,3=9
            stressNZ = neighbourListCoords(n,6+i)*fZ
            configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            If(neighbourListI(n,6).eq.1)Then  ! If atom j (atom B) is in the volume
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + stressNZ
            End If
! Atom i (Atom A) are always in the volume
          End Do
        End If
      End Do
! Multiply stress values by factor of volume of domain
      Do i=1,3
        Do j=1,3
          k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ...
          configCalcStresses(configID,k) = &
          (0.5D0/configVolume(configID))*configCalcStresses(configID,k)
          configCalcStresses(configID,k) = &
          UnitConvert(configCalcStresses(configID,k),"EVAN3","GPA")
        End Do
      End Do
    End If
  End Subroutine calcEnergy_TBEAM

  Subroutine resetESF()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: i, j
! Clears entire useable energies, stress, force and atom energy arrays
! Energy
    Do i=1,configCount
      configCalcEnergies(i) = 0.0D20
    End Do
! Stresses
    Do i=1,configCount
      Do j=1,9
        configCalcStresses(i,j) = 0.0D20
      End Do
    End Do
! Forces
    Do i=1,configsAtomTotal
      Do j=1,3
        configCalcForces(i,j) = 0.0D20
      End Do
    End Do
! Atom energies
    Do i=1,configsAtomTotal
      Do j=1,2
        configAtomEnergy(i,j) = 0.0D0
      End Do
    End Do
  End Subroutine resetESF

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

End Module calcEAM
