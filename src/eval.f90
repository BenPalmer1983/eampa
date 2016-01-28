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
! Calculate config energies,
    Call calcEnergies()  ! calcEAM.f90
! Loop through configs and calculate RSS
    totalRSS = 0.0D0        
    crCount = 0
    rssConfigsArrTotal%total = 0.0D0
    rssConfigsArrTotal%energy = 0.0D0
    rssConfigsArrTotal%force = 0.0D0
    rssConfigsArrTotal%stress = 0.0D0
    Do configID=1,configCount
! RSS
      Call evalEAM_RSS(configID)
! Totals      
      rssConfigsArrTotal%total = rssConfigsArrTotal%total+rssConfigsArr(configID)%total
      rssConfigsArrTotal%energy = rssConfigsArrTotal%energy+rssConfigsArr(configID)%energy
      rssConfigsArrTotal%force = rssConfigsArrTotal%force+rssConfigsArr(configID)%force
      rssConfigsArrTotal%stress = rssConfigsArrTotal%stress+rssConfigsArr(configID)%stress
      !print *,"Total  ",rssConfigsArrTotal%total
      !print *,"Energy ",rssConfigsArrTotal%energy
      !print *,"Force  ",rssConfigsArrTotal%force
      !print *,"Stress ",rssConfigsArrTotal%stress
! Total RSS
      totalRSS = totalRSS + rssConfigsArr(configID)%total
      !print *,configID,configCalcEnergies(configID)
    End Do
! Distribute value
    !Call M_distDouble(totalRSS)
  End Subroutine evalEAM
! ---------------------------------------------------------------------------------------------------
  Subroutine evalEAM_RSS(configID)
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID, i, j, coordStartG, coordEndG, configAtoms
    Real(kind=DoubleReal) :: energyRSS, stressRSS, forceRSS
    Real(kind=DoubleReal) :: refVal, calcVal
    Logical :: calcRSS
! Init    
    configAtoms = configurationCoordsKeyG(configID,2)
! Clear config RSS arrays
    energyRSS = 0.0D0
    stressRSS = 0.0D0
    forceRSS = 0.0D0
! Energy
    If(configRefEnergies(configID).gt.-2.0D20.and.configCalcEnergies(configID).gt.-2.0D20)Then
      energyRSS = rssWeighting(1)*&
      (configRefEnergies(configID)*configAtoms-&
      configCalcEnergies(configID)*configAtoms)**2
      crCount = crCount + 1
      calcRef(crCount,1) = configCalcEnergies(configID)
      calcRef(crCount,2) = configRefEnergies(configID)
    End If
! Forces
    calcRSS = .true.
    coordStartG = configurationCoordsKeyG(configID,1)
    coordEndG = configurationCoordsKeyG(configID,3)
    Do i=coordStartG,coordEndG  ! check that all ref-calc forces exist
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
      refVal = 0.0D0
      calcVal = 0.0D0
      Do i=coordStartG,coordEndG
        Do j=1,3
          forceRSS = forceRSS+&
          (configRefForces(i,j)-configCalcForces(i,j))**2
          refVal = refVal + abs(configRefForces(i,j))
          calcVal = calcVal + abs(configCalcForces(i,j))
! R Array
          crCount = crCount + 1
          calcRef(crCount,1) = configCalcForces(i,j)
          calcRef(crCount,2) = configRefForces(i,j)
        End Do
      End Do
      forceRSS = rssWeighting(2)*forceRSS
! 
      crCount = crCount + 1
      calcRef(crCount,1) = calcVal
      calcRef(crCount,2) = refVal
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
        (configRefStresses(configID,i)-&
        configCalcStresses(configID,i))**2 ! Calc rss in eV/ang^3
        !print *,configRefStresses(configID,i),configCalcStresses(configID,i)
      End Do
      stressRSS = rssWeighting(3)*stressRSS
    End If
! Save rss values
    rssConfigsArr(configID)%energy = energyRSS*configWeighting(configID)
    rssConfigsArr(configID)%force = forceRSS*configWeighting(configID)
    rssConfigsArr(configID)%stress = stressRSS*configWeighting(configID)
    rssConfigsArr(configID)%total = energyRSS + forceRSS + stressRSS    
    !print *,rssConfigsArr(configID)%total
  End Subroutine evalEAM_RSS

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

End Module eval
