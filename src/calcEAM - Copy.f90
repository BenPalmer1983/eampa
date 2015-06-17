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
    configCalcEnergies = -2.1D20
    configCalcForces = -2.1D20
    configCalcStresses = -2.1D20
! MPI Processes per energy calculation
! Loop through configurations
    forceStressSwitch = 1
    Do configID=1,configCount
      If(processMap(configID).eq.mpiProcessID)Then  ! build in more mpi options here
        Call calcEnergy(configID, configEnergy)
        configCalcEnergies(configID) = configEnergy  !Store result from primary mpi process
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
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: aType, bType, aID_R, bID_R, aID_A, bID_A
    Integer(kind=StandardInteger) :: aID, bID, aRhoID, bRhoID
    Integer(kind=StandardInteger) :: configID, n, nKey, i, j, k
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff
    Real(kind=DoubleReal) :: rIJ, fX, fY, fZ
    Real(kind=DoubleReal) :: stressNX, stressNY, stressNZ
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: forceArr
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Integer(kind=StandardInteger) :: forceKeyA, forceKeyB
    Real(kind=DoubleReal) :: densDerivAB, densDerivBA, embeDerivA, embeDerivB, forceM
    Real(kind=DoubleReal) :: timeStartEFS, timeEndEFS
    Integer(kind=StandardInteger) :: processesPerEnergy, primaryProcess, embKeyInc, processID
    Integer(kind=StandardInteger), Dimension(1:10,1:2) :: procArrayEmbe
    Integer(kind=StandardInteger), Dimension(1:1000) :: neighbourCount 
!-----------------
! Init 
!-----------------
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
    neighbourCount = 0
! Init arrays
    Do i=1,configLength
      calculationDensity(i) = 0.0D0     ! density at each atom, length = number of atoms in configuration
    End Do
    If(forceStressSwitch.eq.1)Then
      Do i=1,nlLength
        pairForce(i) = 0.0D0            ! Pair forces, length of neighbour list
      End Do
      Do i=configStart,configEnd        ! store directly into global array
        configCalcForces(i,1) = 0.0D0   ! Force x component atom i
        configCalcForces(i,2) = 0.0D0   ! Force y component atom i
        configCalcForces(i,3) = 0.0D0   ! Force z component atom i
        configAtomEnergy(i,1) = 0.0D0   ! Atom energy from pair
        configAtomEnergy(i,2) = 0.0D0   ! Atom energy from embedding
      End Do
! Init stress array for this config
      Do i=1,9                                  ! store directly into global array
        configCalcStresses(configID,i) = 0.0D0  ! xx, xy...zy, zz for config i
      End Do
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
        yArray = SearchPotentialPoint(aType,bType,1,eamType,rIJ)
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
        yArray = SearchPotentialPoint(bType,0,2,eamType,rIJ)
        calculationDensity(aID_R) = calculationDensity(aID_R) + yArray(1)
! Electron density each B is embedded in due to the electrons of the As around it
        yArray = SearchPotentialPoint(aType,0,2,eamType,rIJ)
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
      yArray = SearchPotentialPoint(aType,0,3,eamType,calculationDensity(aID_R))
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
          yArray = SearchPotentialPoint(aType,0,3,eamType,calculationDensity(aID_R))
          embeDerivA = yArray(2)
! @Pji(r)/@r
          yArray = SearchPotentialPoint(bType,0,2,eamType,rIJ)
          densDerivBA = yArray(2)
! @Fj(p)/@p
          yArray = SearchPotentialPoint(bType,0,3,eamType,calculationDensity(bID_R))
          embeDerivB = yArray(2)
! @Pij(r)/@r
          yArray = SearchPotentialPoint(aType,0,2,eamType,rIJ)
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
