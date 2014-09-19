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


  Subroutine calcEnergies()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID, forceKeyStart, forceKeyEnd, selectedProcess
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: energyTimeStart, energyTimeEnd
! Start time
    Call cpu_time(energyTimeStart)
! Init variables    
    configCalcEnergies = -2.1D20 
    configCalcForces = -2.1D20
    configCalcStresses = -2.1D20  
! Loop through configurations
    Do configID=1,configCount
      If(processMap(configID,1).eq.mpiProcessID)Then
        Call calcEnergy(configID, configEnergy, 1)
        configCalcEnergies(configID) = configEnergy
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
! End time
    Call cpu_time(energyTimeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("E-F-S Calc Configs",energyTimeEnd-energyTimeStart)      
  End Subroutine calcEnergies 
  
  

  Subroutine calcEnergy(configID, configEnergy, forceCalcIn)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID, n, i, j, k  
    Integer(kind=StandardInteger) :: nlStart, nlEnd
    Real(kind=DoubleReal) :: configEnergy
    Real(kind=DoubleReal) :: rCutoff
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:3) :: yArray  
    Real(kind=DoubleReal), Dimension(1:3) :: forceArr
    Integer(kind=StandardInteger), Optional :: forceCalcIn
    Integer(kind=StandardInteger) :: forceCalc, forceKeyA, forceKeyB
    Real(kind=DoubleReal) :: densDerivAB, densDerivBA, embeDerivA, embeDerivB, forceM
    Real(kind=DoubleReal) :: timeStartEFS, timeEndEFS
! Start Time
    Call cpu_time(timeStartEFS)
! Init variables
    rCutoff = configurationsR(configID,11)
    pairEnergy = 0.0D0
    embeddingEnergy = 0.0D0
    configEnergy = 0.0D0
    nlStart = neighbourListKey(configID,1)
    nlEnd = neighbourListKey(configID,3)
    calculationDensity = 0.0D0    
    forceCalc = 0
    If(Present(forceCalcIn))Then
      forceCalc = forceCalcIn
    End If
    If(forceCalc.eq.1)Then
! Init force array for this config    
      Do i=configurationCoordsKeyG(configID,1),configurationCoordsKeyG(configID,3)
        configCalcForces(i,1) = 0.0D0       
        configCalcForces(i,2) = 0.0D0      
        configCalcForces(i,3) = 0.0D0      
      End Do
! Init stress array for this config      
      Do i=1,9
        configCalcStresses(configID,i) = 0.0D0
      End Do  
    End If
!--------------------------------------------------
! Loop 1 - sum pair energy and density	
!--------------------------------------------------
    Do n=nlStart,nlEnd    
! Check if in cutoff radius
      If(neighbourListR(n).le.rCutoff)Then
! Pair potential
        yArray = SearchPotentialPoint&
                 (neighbourListI(n,1),neighbourListI(n,2),1,eamType,neighbourListR(n))
        pairEnergy = pairEnergy + yArray(1)
        If(forceCalc.eq.1)Then
! Force from Pair Potential
          forceArr(1) = yArray(2)*neighbourListCoords(n,10)  ! x
          forceArr(2) = yArray(2)*neighbourListCoords(n,11)  ! y
          forceArr(3) = yArray(2)*neighbourListCoords(n,12)  ! z
          forceKeyA = configurationCoordsKeyG(configID,1)+neighbourListI(n,3)-1
          forceKeyB = configurationCoordsKeyG(configID,1)+neighbourListI(n,4)-1
          configCalcForces(forceKeyA,1)=configCalcForces(forceKeyA,1)-forceArr(1)
          configCalcForces(forceKeyA,2)=configCalcForces(forceKeyA,2)-forceArr(2)
          configCalcForces(forceKeyA,3)=configCalcForces(forceKeyA,3)-forceArr(3)
          configCalcForces(forceKeyB,1)=configCalcForces(forceKeyB,1)+forceArr(1)
          configCalcForces(forceKeyB,2)=configCalcForces(forceKeyB,2)+forceArr(2)
          configCalcForces(forceKeyB,3)=configCalcForces(forceKeyB,3)+forceArr(3)
! Store stress          
          If(neighbourListI(n,6).eq.1)Then
            Do i=1,3
              Do j=1,3
                k = 3*(i-1)+j
              End Do
            End Do  
          End If 
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
    Do i=1,configurationCoordsKeyG(configID,2)
      yArray = SearchPotentialPoint&
               (neighbourListI(i,1),0,3,eamType,calculationDensity(neighbourListI(i,1)))
      embeddingEnergy = embeddingEnergy + yArray(1)      
    End Do      
! Sum energies   
    configEnergy = pairEnergy + embeddingEnergy
!--------------------------------------------------
! Loop 3 - remaining derivative values for force	
!--------------------------------------------------    
    If(forceCalc.eq.1)Then
    Do n=nlStart,nlEnd          
! Check if in cutoff radius
      If(neighbourListR(n).le.rCutoff)Then
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
        forceM = (embeDerivA*densDerivAB+embeDerivB*densDerivBA)
        forceArr(1) = forceM*neighbourListCoords(n,10)
        forceArr(2) = forceM*neighbourListCoords(n,11)
        forceArr(3) = forceM*neighbourListCoords(n,12)     
! Force key
        forceKeyA = configurationCoordsKeyG(configID,1)+neighbourListI(n,3)-1
        forceKeyB = configurationCoordsKeyG(configID,1)+neighbourListI(n,4)-1
! Add forces
        configCalcForces(forceKeyA,1)=configCalcForces(forceKeyA,1)-forceArr(1)
        configCalcForces(forceKeyA,2)=configCalcForces(forceKeyA,2)-forceArr(2)
        configCalcForces(forceKeyA,3)=configCalcForces(forceKeyA,3)-forceArr(3)
        configCalcForces(forceKeyB,1)=configCalcForces(forceKeyB,1)+forceArr(1)
        configCalcForces(forceKeyB,2)=configCalcForces(forceKeyB,2)+forceArr(2)
        configCalcForces(forceKeyB,3)=configCalcForces(forceKeyB,3)+forceArr(3)
! Store stress          
          If(neighbourListI(n,6).eq.1)Then
            Do i=1,3
              Do j=1,3
                k = 3*(i-1)+j
              End Do
            End Do  
          End If 
      End If 
    End Do   
    End If
! End Time
    Call cpu_time(timeEndEFS)        
! Store Time    
    Call storeTime(1,timeEndEFS-timeStartEFS)   
  End Subroutine calcEnergy 


  
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