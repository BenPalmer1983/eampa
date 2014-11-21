Module calcEval
!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!
! Evaluate configuration subroutine
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
  Use readEAM
  Use calcEAM  
  Use bulkProperties  
  Use testEAM
  Use output
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: evaluate
  Public :: calcRSS
  Public :: resetCalcVars  
Contains
!---------------------------------------------------------------------------------------------------
! Run Evaluate
!--------------------------------------------------------------------------------------------------- 
  Subroutine evaluate()
    Implicit None   ! Force declaration of all variables
! Private variables      
    Real(kind=DoubleReal) :: timeStartEval, timeEndEval
! Start time
    Call cpu_time(timeStartEval)  
! Reset RSS values  
    totalRSS = 0.0D0  
    eosFitRSS = 0.0D0 
    configRSS = 0.0D0
    configTotalRSS = 0.0D0   
! Just energy, stress, force
    If(eampaRunType(1:4).eq."OPTI".or.eampaRunType(1:4).eq."EVAL")Then      
      Call calcEnergies()      ! Calculate energies, stresses and forces of input configurations
    End If
! Full evaluation options
    If(eampaRunType(1:4).eq."OPTF".or.eampaRunType(1:4).eq."OPTE".or.eampaRunType(1:4).eq."EVAF")Then
      Call calcEnergies()      ! Calculate energies, stresses and forces of input configurations
      Call runTestEAM()        ! Calculate bulk properties for FCC and BCC
    End If
! EAM Testing evaluation options
    If(eampaRunType(1:4).eq."OPTT".or.eampaRunType(1:4).eq."EVAT")Then
      Call runTestEAM()        ! Calculate bulk properties for FCC and BCC
    End If
    Call calcRSS()           ! Calculate RSS of stresses, forces and energies
! Output to file
    If(mpiProcessID.eq.0)Then 
      Open(UNIT=25,FILE=Trim(outputDirectory)//"/rssLog.dat",&
      status="old",position="append",action="write") 
      If(eampaRunType(1:4).eq."OPTI".or.eampaRunType(1:4).eq."EVAL")Then 
        write(25,"(E16.8,E16.8)") &
        totalRSS, configTotalRSS
      End If
      If(eampaRunType(1:4).eq."OPTF".or.eampaRunType(1:4).eq."OPTE".or.eampaRunType(1:4).eq."EVAF")Then
        write(25,"(E16.8,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8)") &
        totalRSS,configTotalRSS,testingALatRSS,testingEMinRSS,&
        testingBMRSS,testingECRSS,eosFitRSS
      End If
      If(eampaRunType(1:4).eq."OPTT".or.eampaRunType(1:4).eq."EVAT")Then
        write(25,"(E16.8,E16.8,E16.8,E16.8,E16.8,E16.8)") &
        totalRSS,testingALatRSS,testingEMinRSS,&
        testingBMRSS,testingECRSS,eosFitRSS
      End If
      If(eampaRunType(1:4).eq."OPES")Then
        write(25,"(E16.8,E16.8)") &
        totalRSS,eosFitRSS
      End If
      Close(25)
    End If   
! Output evaluate results to output file
    Call outputEvaluate()  
! End time
    Call cpu_time(timeEndEval)
! Store Time    
    Call storeTime(2,timeEndEval-timeStartEval)
  End Subroutine evaluate  
!---------------------------------------------------------------------------------------------------
! Reset RSS subroutines
!---------------------------------------------------------------------------------------------------
  Subroutine resetCalcVars()
    configCalc = -2.1D20 
    configCalcForces = -2.1D20
    configCalcStresses = -2.1D20  
  End Subroutine resetCalcVars  
  Subroutine calcRSS()  
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: i, j, configID, fsKey, feKey, xSize, ySize
! Reset variables
    configRSS = 0.0D0
    totalRSS = 0.0D0
    configTotalRSS = 0.0D0    
    If(eampaRunType(1:4).eq."OPTI".or.eampaRunType(1:4).eq."EVAI".or.&
    eampaRunType(1:4).eq."OPTF".or.eampaRunType(1:4).eq."EVAF".or.&
    eampaRunType(1:4).eq."OPTE")Then
! Loop through configs    
      Do configID=1,configCount
! Energy RSS - all atoms
        If(configCalcEnergies(configID).gt.-2.0D20.and.configRefEnergies(configID).gt.-2.0D20)Then
          configRSS(configID,1) = rssWeighting(1)*&
          (configCalcEnergies(configID)-&
          configRefEnergies(configID)*configurationCoordsKeyG(configID,2))**2
          configRSS(configID,1) = configWeighting(configID)*configRSS(configID,1)
        End If
! Force RSS - all atoms
        fsKey = configurationCoordsKeyG(configID,1)
        feKey = configurationCoordsKeyG(configID,3)
        If(configRefForces(fsKey,1).gt.-2.0D20.and.configCalcForces(fsKey,1).gt.-2.0D20)Then    
        Do i=fsKey,feKey
            configRSS(configID,2) = configRSS(configID,2) + &
            (configRefForces(i,1)-configCalcForces(i,1))**2 + &
            (configRefForces(i,2)-configCalcForces(i,2))**2 + &
            (configRefForces(i,3)-configCalcForces(i,3))**2
          End Do
          configRSS(configID,2) = rssWeighting(2)*configRSS(configID,2)
          configRSS(configID,2) = configWeighting(configID)*configRSS(configID,2)
        End If
! Stress - volume of atoms     
        If(configRefStresses(configID,1).gt.-2.0D20.and.configCalcStresses(configID,1).gt.-2.0D20)Then    
          Do i=1,9
            configRSS(configID,3) = configRSS(configID,3) + &
            (configCalcStresses(configID,i)-configRefStresses(configID,i))**2
          End Do
          configRSS(configID,3) = rssWeighting(3)*configRSS(configID,3)
          configRSS(configID,3) = configWeighting(configID)*configRSS(configID,3)
        End If    
      End Do
! Total RSS
      xSize = size(configRSS,1)
      ySize = size(configRSS,2)
      Do i=1,xSize
        configRSS(i,ySize) = 0    ! Last column is total RSS for config
        Do j=1,(ySize-1)
          totalRSS = totalRSS + configRSS(i,j)
          configRSS(i,ySize) = configRSS(i,ySize) + configRSS(i,j)
        End Do
      End Do 
      configTotalRSS = totalRSS
    End If
! Apply weightings to testing RSS
    If(eampaRunType(1:4).eq."OPTF".or.eampaRunType(1:4).eq."EVAF".or.&
    eampaRunType(1:4).eq."OPTT".or.eampaRunType(1:4).eq."EVAT".or.&
    eampaRunType(1:4).eq."OPTE")Then
      testingALatRSS = testingALatRSS*rssWeighting(4)
      testingEMinRSS = testingEMinRSS*rssWeighting(5)   
      testingBMRSS = testingBMRSS*rssWeighting(6)   
      testingECRSS = testingECRSS*rssWeighting(7) 
      eosFitRSS = eosFitRSS*rssWeighting(8) 
! Add to total    
      totalRSS = totalRSS+testingALatRSS+testingEMinRSS+testingBMRSS+&
                 testingECRSS+eosFitRSS
    End If  
    If(eampaRunType(1:4).eq."OPES")Then
      totalRSS = totalRSS+eosFitRSS
    End If
  End Subroutine calcRSS 
End Module calcEval  