Module calcEval

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Read user input file 

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
  
  
  
 
  
!------------------------------------
! Run Evaluate
!------------------------------------  
  
  Subroutine evaluate()
    Implicit None   ! Force declaration of all variables
! Private variables      
    Real(kind=DoubleReal) :: timeStartEval, timeEndEval
! Start time
    Call cpu_time(timeStartEval)    
! Run energy/force/stress calculations
    Call calcEnergies()  
! Full evaluation options
    If(eampaRunType(1:4).eq."OPTF".or.eampaRunType(1:4).eq."EVAF")Then
! Run bulk property calculations
      Call calcEquilibrium()
    
    End If
   
! Calculate RSS between reference and calculated values
    Call calcRSS()
! Output to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"RSS: ",totalRSS
    End If
    
! Output evaluate results to output file
    Call outputEvaluate()  
! End time
    Call cpu_time(timeEndEval)
! Store Time    
    Call storeTime(2,timeEndEval-timeStartEval)
  End Subroutine evaluate  


  
  
!------------------------------------
! Reset RSS subroutines
!------------------------------------
  
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
! Loop through configs    
    Do configID=1,configCount
! Energy RSS
      If(configCalcEnergies(configID).gt.-2.0D20.and.configRefEnergies(configID).gt.-2.0D20)Then
        configRSS(configID,1) = rssWeighting(1)*&
        (configCalcEnergies(configID)-&
        configRefEnergies(configID)*configurationCoordsKeyG(configID,2))**2
      End If
! Force RSS  
      fsKey = configurationCoordsKeyG(configID,1)
      feKey = configurationCoordsKeyG(configID,3)
      !If(configRefForces(fsKey,1).gt.-2.0D20.and.configCalcForces(fsKey,1).gt.-2.0D20)Then    
      !  Do i=fsKey,feKey
      !    configRSS(configID,2) = configRSS(configID,2) + &
      !    (configRefForces(i,1)-configCalcForces(i,1))**2 + &
      !    (configRefForces(i,2)-configCalcForces(i,2))**2 + &
      !    (configRefForces(i,3)-configCalcForces(i,3))**2
      !  End Do
      !  configRSS(configID,2) = rssWeighting(2)*configRSS(configID,2)
      !  print *,configRSS(configID,2),rssWeighting(2)
      !End If
! Eq Vol RSS      
      If(configCalcEV(configID).gt.-2.0D20.and.configRefEV(configID).gt.-2.0D20)Then
        configRSS(configID,4) = rssWeighting(4)*&
        (configCalcEV(configID)-configRefEV(configID))**2
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
! Print out results    
    If(mpiProcessID.eq.0)Then
      Do i=1,configCount
        !print *,"Config: ",i,configRSS(i,ySize)
      End Do
    End If
    
    !configRef = -2.1D20                    
    !configCalc = -2.1D20                      
    !configRefForces = -2.1D20               
    !configCalcForces = -2.1D20
    !configRefStresses = -2.1D20
    !configCalcStresses = -2.1D20      
  End Subroutine calcRSS 
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  

    
  
End Module calcEval  