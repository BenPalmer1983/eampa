Module calcEAM
!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!
! Subroutines and functions to calculate energy, forces and stresses
!----------------------------------------
! Updated: 12th Aug 2014
!----------------------------------------
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
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: calcEnergies
  Public :: calcEnergy  
Contains
!---------------------------------------------------------------------------------------------------
  Subroutine calcEnergies()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID, forceKeyStart, forceKeyEnd, selectedProcess
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: energyTimeStart, energyTimeEnd
    !Integer(kind=StandardInteger) :: processesPerEnergy
! Start time
    Call cpu_time(energyTimeStart)
! Init variables    
    configCalcEnergies = -2.1D20 
    configCalcForces = -2.1D20
    configCalcStresses = -2.1D20 
! MPI Processes per energy calculation    
! Loop through configurations
    Do configID=1,configCount
      If(processMap(configID,1).eq.mod(mpiProcessID,configCount))Then  ! build in more mpi options here
        Call calcEnergy(configID, configEnergy, 1)
        If(processMap(configID,1).eq.mpiProcessID)Then
          configCalcEnergies(configID) = configEnergy  !Store result from primary mpi process
        End If
      End If  
    End Do
! Distribute energy array    
    Call M_collDouble1D(configCalcEnergies)
    Call M_distDouble1D(configCalcEnergies)
! Distribute force array    
    Do configID=1,configCount
      selectedProcess = processMap(configID,1)
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
    Call cpu_time(energyTimeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("E-F-S Calc Configs",energyTimeEnd-energyTimeStart)      
  End Subroutine calcEnergies 
!---------------------------------------------------------------------------------------------------
  Subroutine calcEnergy(configID, configEnergy, forceCalcIn, pairEnergyOut, embeddingEnergyOut)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID, n, nKey, i, j, k 
    Integer(kind=StandardInteger) :: nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Optional :: pairEnergyOut, embeddingEnergyOut
    Real(kind=DoubleReal), Dimension(1:3) :: yArray  
    Real(kind=DoubleReal), Dimension(1:3) :: forceArr
    Integer(kind=StandardInteger), Optional :: forceCalcIn
    Integer(kind=StandardInteger) :: forceCalc, forceKeyA, forceKeyB
    Real(kind=DoubleReal) :: densDerivAB, densDerivBA, embeDerivA, embeDerivB, forceM
    Real(kind=DoubleReal) :: timeStartEFS, timeEndEFS
    Integer(kind=StandardInteger) :: processesPerEnergy, primaryProcess, embKeyInc, processID
    Integer(kind=StandardInteger), Dimension(1:10,1:2) :: procArrayEmbe
! Start Time
    Call cpu_time(timeStartEFS)
! Init optional variables    
    forceCalc = 0
    If(Present(forceCalcIn))Then
      forceCalc = forceCalcIn
    End If
    If(Present(pairEnergyOut))Then
      pairEnergyOut = 0.0D0
    End If  
    If(Present(embeddingEnergyOut))Then
      embeddingEnergyOut = 0.0D0
    End If  
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
      calculationDensity(i) = 0.0D0    
    End Do  
! Set force arrays    
    If(forceCalc.eq.1)Then
! Init force arrays for this config
      Do i=1,nlLength
        pairForce(i) = 0.0D0
      End Do
      Do i=configStart,configEnd
        configCalcForces(i,1) = 0.0D0       
        configCalcForces(i,2) = 0.0D0      
        configCalcForces(i,3) = 0.0D0      
      End Do
! Init stress array for this config      
      Do i=1,9
        configCalcStresses(configID,i) = 0.0D0
      End Do  
    End If
! Multiple processes per energy calculation
    If(configCount.gt.0)Then
      processesPerEnergy = Floor(1.0D0*(mpiProcessCount/configCount))
    Else
      processesPerEnergy = 1
    End If
    primaryProcess = processMap(configID,1)
    procArrayEmbe = 0
! Set embe energy calc start-end array
    If(processesPerEnergy.eq.1)Then
      procArrayEmbe(1,1) = 1
      procArrayEmbe(1,2) = configLength
    Else
      embKeyInc = Floor(1.0D0*configLength/processesPerEnergy)
      procArrayEmbe(1,1) = 1
      procArrayEmbe(1,2) = embKeyInc
      Do i=2,processesPerEnergy
        If(i.eq.processesPerEnergy)Then
          procArrayEmbe(i,1) = procArrayEmbe(i-1,2)+1
          procArrayEmbe(i,2) = configLength
        Else  
          procArrayEmbe(i,1) = procArrayEmbe(i-1,2)+1
          procArrayEmbe(i,2) = procArrayEmbe(i-1,2)+embKeyInc          
        End If
      End Do  
    End If      
!--------------------------------------------------
! Loop 1 - sum pair energy and density	
!--------------------------------------------------
    nKey = 0
    Do n=nlStart,nlEnd    
! Check if in cutoff radius
      If(neighbourListR(n).le.rCutoff)Then
! Pair potential
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,1),neighbourListI(n,2),1,eamType,neighbourListR(n))
        pairEnergy = pairEnergy + yArray(1)
        If(forceCalc.eq.1)Then
! Force from Pair Potential
          nKey = nKey + 1
          pairForce(nKey) = yArray(2)
        End If          
! Electron density of each A due to the electrons of the Bs around it        
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,2),0,2,eamType,neighbourListR(n))
        calculationDensity(neighbourListI(n,3)) = &
                 calculationDensity(neighbourListI(n,3)) + yArray(1)     
      End If
    End Do 
!--------------------------------------------------
! Loop 2 - embedding energy	
!--------------------------------------------------   
    If(mpiEnergy.eq.1)Then
      Do j=1,processesPerEnergy
        processID = mod(mpiProcessID,configCount)+(j-1)*configCount
        If(mpiProcessID.eq.processID)Then
          Do i=procArrayEmbe(j,1),procArrayEmbe(j,2)
            yArray = SearchPotentialPoint&
               (neighbourListI(i,1),0,3,eamType,calculationDensity(neighbourListI(i,1)))
            embeddingEnergy = embeddingEnergy + yArray(1)      
          End Do        
        End If
      End Do 
      If(processesPerEnergy.gt.1)Then
        Call M_sumDouble(mod(mpiProcessID,configCount), mpiProcessID, embeddingEnergy)        
      End If 
    Else
      Do i=1,configLength
        yArray = SearchPotentialPoint&
          (neighbourListI(i,1),0,3,eamType,calculationDensity(neighbourListI(i,1)))
          embeddingEnergy = embeddingEnergy + yArray(1)      
      End Do
    End If
! Sum energies   
    configEnergy = pairEnergy + embeddingEnergy
! output components if required
    If(Present(pairEnergyOut))Then
      pairEnergyOut = pairEnergy
    End If  
    If(Present(embeddingEnergyOut))Then
      embeddingEnergyOut = embeddingEnergy
    End If     
!--------------------------------------------------
! Loop 3 - remaining derivative values for force	
!--------------------------------------------------    
    If(forceCalc.eq.1)Then
    nKey = 0    
    Do n=nlStart,nlEnd          
! Check if in cutoff radius
      If(neighbourListR(n).le.rCutoff)Then
        nKey = nKey + 1
! Density derivatives
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,2),0,2,eamType,neighbourListR(n))
        densDerivAB = yArray(2)     !derivative of the density  
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,1),0,2,eamType,neighbourListR(n))
        densDerivBA = yArray(2)  
! Embedding derivatives
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,1),0,3,eamType,neighbourListR(n))
        embeDerivA = yArray(2)   
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,2),0,3,eamType,neighbourListR(n))
        embeDerivB = yArray(2)         
! Force, embedding dens A
        forceM = -1.0D0*(pairForce(nKey)+embeDerivA*densDerivAB+embeDerivB*densDerivBA)
        forceArr(1) = forceM*neighbourListCoords(n,10)
        forceArr(2) = forceM*neighbourListCoords(n,11)
        forceArr(3) = forceM*neighbourListCoords(n,12)     
! Force key
        forceKeyA = (configStart-1)+neighbourListI(n,3)
        forceKeyB = (configStart-1)+neighbourListI(n,4)
! Add force to atom A, subtract from atom B
        configCalcForces(forceKeyA,1)=configCalcForces(forceKeyA,1)+forceArr(1)
        configCalcForces(forceKeyA,2)=configCalcForces(forceKeyA,2)+forceArr(2)
        configCalcForces(forceKeyA,3)=configCalcForces(forceKeyA,3)+forceArr(3)
        configCalcForces(forceKeyB,1)=configCalcForces(forceKeyB,1)-forceArr(1)
        configCalcForces(forceKeyB,2)=configCalcForces(forceKeyB,2)-forceArr(2)
        configCalcForces(forceKeyB,3)=configCalcForces(forceKeyB,3)-forceArr(3)
! Store stress          
        If(neighbourListI(n,6).eq.1)Then  !Only contributions of forces ON atoms within the volume
          Do i=1,3
            Do j=1,3
              k = 3*(i-1)+j  ! 1,1=1 1,2=2 1,3=3 2,1=4 ...
              configCalcStresses(configID,k) = configCalcStresses(configID,k) + &
              neighbourListCoords(n,6+i)*forceArr(j)  
            End Do
          End Do  
        End If 
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
! End Time
    Call cpu_time(timeEndEFS)        
! Store Time    
    Call storeTime(6,timeEndEFS-timeStartEFS)   
  End Subroutine calcEnergy 
!---------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
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
  End function SearchPotentialPoint    
!---------------------------------------------------------------------------------------------------
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
  End function FunctionKey
End Module calcEAM  