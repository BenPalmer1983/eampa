Module eval
! --------------------------------------------------------------!
! Evaluate configurations/EAM
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calls the calcEAM subroutines and evaluates the results
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
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: evalEAM
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine evalEAM()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID
    
    !Real(kind=DoubleReal) :: configEnergy
   
! Calculate config energies, 
    Call calcEnergies()
! Loop through configs
    totalRSS = 0.0D0
    Do configID=1,configCount
      Call evalEAM_RSS(configID)
      totalRSS = totalRSS + configRSS(configID,10)
    End Do
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Total RSS: ",totalRSS
    End If
    
  End Subroutine evalEAM

  
  Subroutine evalEAM_RSS(configID)
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID, i, j, coordStartG, coordEndG
    Real(kind=DoubleReal) :: energyRSS, stressRSS, forceRSS
    Logical :: calcRSS
! Clear config RSS arrays
    energyRSS = 0.0D0   
    stressRSS = 0.0D0   
    forceRSS = 0.0D0    
    Do i=1,10
      configRSS(configID,i) = 0.0D0                      ! 1 energy, 2 forces, 3 stresses, 10 All
    End Do
! Energy
    If(configRefEnergies(configID).gt.-2.0D20.and.configCalcEnergies(configID).gt.-2.0D20)Then
      energyRSS = rssWeighting(1)*&
      (configRefEnergies(configID)-configCalcEnergies(configID))**2
    End If
! Stress
    calcRSS = .true.
    Do i=1,9
      If(configRefStresses(configID,i).lt.-2.0D20)Then
        calcRSS = .false.
      End If
      If(configCalcStresses(configID,i).lt.-2.0D20)Then
        calcRSS = .false.
      End If
    End Do
    If(calcRSS)Then
      Do i=1,9
        stressRSS = stressRSS+&
        (configRefStresses(configID,i)-configCalcStresses(configID,i))**2
      End Do
      stressRSS = rssWeighting(2)*stressRSS
    End If
! Forces
    calcRSS = .true.
    coordStartG = configurationCoordsKeyG(configID,1)
    coordEndG = configurationCoordsKeyG(configID,3)
    Do i=coordStartG,coordEndG
      Do j=1,3
        If(configRefForces(i,j).lt.-2.0D20)Then
          calcRSS = .false.
        End If
        If(configCalcForces(i,j).lt.-2.0D20)Then
          calcRSS = .false.
        End If
      End Do
    End Do
    If(calcRSS)Then
      Do i=coordStartG,coordEndG
        Do j=1,3
          forceRSS = forceRSS+&
          (configRefForces(i,j)-configCalcForces(i,j))**2
        End Do
      End Do
      forceRSS = rssWeighting(3)*forceRSS
    End If
! Save rss values
    configRSS(configID,1) = energyRSS
    configRSS(configID,2) = stressRSS
! Sum config RSS and store in slot 10
    Do i=1,9
      configRSS(configID,10) = configRSS(configID,10) + configRSS(configID,i)                       ! 1 energy, 2 forces, 3 stresses, 10 All
    End Do
  End Subroutine evalEAM_RSS
  
! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!


End Module eval
