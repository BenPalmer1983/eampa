Module opti
! --------------------------------------------------------------!
! Optimise EAM functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calls the eval and calcEAM subroutines to optimise potential functions
! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------
! Setup Modules
  Use kinds
  Use types
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use plot
  Use globals
  Use initialise
  Use loadData
  Use output
  Use readEAM
  Use analytic
  Use makePotential
  Use eamGen
  Use eamOpt
  Use calcEAM
  Use bpCalcEAM
  Use eval
  Use evalBP
  Use output
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: optiEAM
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine optiEAM()
! Assumptions/Requirements
! Potential will always have a ZBL hard core
! Last node will be y(x)=0, y'(x)=0, y''(x)=0 for pair and density
! First and last node of pair functions are FIXED
! Last node of density functions are FIXED
!   
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(saConfig) :: saConfigLive
    Integer(kind=StandardInteger) :: i
    Integer(kind=StandardInteger) :: LM_nodes
! output
    If(TerminalPrint())Then
      print *,""
      Call printBR(90,"=")
      print *,"                           Optimise"
      Call printBR(90,"=")
      print *,""
      print *,""
    End If    
    
! Record input potential
      
    Call saveEamFile("P_01_input.pot")    ! output.f90
    Call saveEamNodes("P_01_input.nodes") ! output.f90
    Call eamCharts("P_01_input")          ! readEAM.f90
    Call runCalcs(0,optRunType,.false.)             ! use input potential and output EoS charts runCalcs(splineNodesIn,runTypeIn,eosChartIn)
    
    Print *,""
    Print *,""
    print *,"==================================================================="
    print *,"Input Potential RSS"
    print *,"==================================================================="
    If(TerminalPrint())Then
      print *,"RSS: ",totalRSS
    End If  
    Call runCalcs(0,1,.false.)             ! use input potential and output EoS charts runCalcs(splineNodesIn,runTypeIn,eosChartIn)      
! Output Bulk Properties       
    Call outputBpT()  

    
        
    
! ----------------------------------
! Make an EAM to start from
! ----------------------------------    
    If(optFrom.eq.1)Then
      ! do nothing, use input potential
      If(TerminalPrint())Then
        print *,""
        print *,"Using input potential as starting point"
        print *,""
      End If  
    End If
    If(optFrom.eq.2)Then
      ! make a potential
      If(TerminalPrint())Then
        print *,""
        print *,"Making a starting potential"
        print *,""
      End If  
      Call MakeStartingPotential() 
    End If
    !Call saveEamFile("opt_001.pot")  
    !Call eamCharts("P_02_start")
    
! ----------------------------------
! Set spline nodes, 
! ----------------------------------    
    
    Call setEamNodes()                         ! set spline nodes from eamKey/eamData   
    optEmbeddingFit = 0                        ! Force embedding functional to fit form - on    
    optDensityFit = 0                          ! Spline poly for density - off   
    Call setEamSpline(.true.)                  ! make EAM potential by splining nodes
    !Call RemoveIrrelevantPoints() 
    
    
    Call saveEamFile("P_02_start.pot")
    Call saveEamNodes("P_02_input.nodes")      ! output.f90
    Call eamCharts("P_02_start")
    
    
    Print *,""
    Print *,""
    print *,"==================================================================="
    print *,"Splined Input Potential RSS"
    print *,"==================================================================="
    Call runCalcs(1,optRunType,.false.)
    If(TerminalPrint())Then
      print *,"RSS: ",totalRSS
    End If    
! Output Bulk Properties    
    Call runCalcs(1,1,.false.)   
    Call outputBpT() 
    
    
     
! Print
    print *,"----------------------------------"
    print *,"Starting Optimisation"
    print *,"----------------------------------"
    ! ZBL + exp(poly) spline + polyspline for pair
    optEmbeddingFit = 0    ! Force embedding functional to fit form     
    optDensityFit = 0      ! Spline poly for density
    Call runCalcs(1,optRunType)     ! Call bulk properties
    If(TerminalPrint())Then
      print *,"Spline converted input RSS ",totalRSS
    End If  
    
! -----------------------------------------    
! Simulated annealing    
! -----------------------------------------    
    
    ! ZBL + exp(poly) spline + polyspline for pair
    
    
    
    saConfigLive = saConfigIn(1)
    Call saOpt_A(saConfigLive)   ! Analytic

    !saConfigLive = saConfigIn
    !Call saOpt(saConfigLive)
    
    
    
    
    
    
  End Subroutine optiEAM
  
   
  
! --------------------------------------------------------------------------------------------------- 
  Subroutine runCalcs(splineNodesIn,runTypeIn,eosChartIn)
! Calculate stress/energy/force of configuration/s   
    Implicit None   ! Force declaration of all variables
! Private variables 
    Integer(kind=StandardInteger), Optional :: splineNodesIn
    Integer(kind=StandardInteger), Optional :: runTypeIn
    Logical, Optional :: eosChartIn
    Integer(kind=StandardInteger) :: splineNodes
    Integer(kind=StandardInteger) :: runType
    Logical :: eosChart
! Optional Arguments
    splineNodes = 0
    If(Present(splineNodesIn))Then
      splineNodes = splineNodesIn
    End If   
    runType = 0
    If(Present(runTypeIn))Then
      runType = runTypeIn
    End If  
    eosChart = .false.
    If(Present(eosChartIn))Then
      eosChart = eosChartIn
    End If 
! Make spline    
    If(splineNodes.eq.1)Then
      Call setEamSpline(.true.)  ! readEAM.f90
    End If        
! Turn off verbose    
    quietOverride = .true.     
    Call evalEAM()               ! eval.f90
! Bulk Properties 
    If(runType.ge.1)Then
      Call evalBulkProperties(.false.,eosChart)  ! evalBP.f90
    End If
    quietOverride = .false. 
! Save data to file    
    Call outputOptLine()
  End Subroutine runCalcs 
  
  Subroutine calcRefRSS(crRSS)
! Must run after runCalcs
    Implicit None   ! Force declaration of all variables
! Private variables 
    Real(kind=DoubleReal) :: crRSS
    Integer(kind=StandardInteger) :: i
    crRSS = 0.0D0
    Do i=1,crCount
      crRSS = crRSS + (calcRef(i,1)-calcRef(i,2))**2
    End Do
  End Subroutine calcRefRSS 
  
  Subroutine completeNodeDerivs()
! Must run after runCalcs
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Loop through functions/nodes
    Do i=1,eamFunctionCount
      Do j=splineNodesKey(i,4),splineNodesKey(i,6)      
        If(splineNodesKey(i,3).eq.1)Then ! pair
! Do not alter satr node, start node + 1 or end node y'(x) and y''(x)       
          If(j.ne.splineNodesKey(i,4).and.j.ne.(splineNodesKey(i,4)+1).and.j.ne.splineNodesKey(i,6))Then
            yArray = PointInterp(splineNodesData,splineNodesData(j,1),3,2,splineNodesKey(i,4),splineNodesKey(i,5))
            splineNodesData(j,3) = yArray(2)  
            splineNodesData(j,4) = yArray(3)
          End If
        End If 
        If(splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)Then ! dens/dden/sden
          If(j.ne.splineNodesKey(i,6))Then
            yArray = PointInterp(splineNodesData,splineNodesData(j,1),4,2,splineNodesKey(i,4),splineNodesKey(i,5))  
            splineNodesData(j,3) = yArray(2)  
            splineNodesData(j,4) = yArray(3)
          End If
        End If   
        If(splineNodesKey(i,3).eq.3.or.splineNodesKey(i,3).eq.6.or.splineNodesKey(i,3).eq.7)Then ! embe/demb/semb
          yArray = PointInterp(splineNodesData,splineNodesData(j,1),4,2,splineNodesKey(i,4),splineNodesKey(i,5))  
          splineNodesData(j,3) = yArray(2)  
          splineNodesData(j,4) = yArray(3)   
        End If      
      End Do
    End Do
  End Subroutine completeNodeDerivs 
! ---------------------------------------------------------------------------------------------------  
!
! Simulated Annealing Functions 
!  
! --------------------------------------------------------------------------------------------------- 
  Subroutine saOpt(saConfigLive)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Type(saConfig) :: saConfigLive
! Start SA
    If(TerminalPrint())Then
      print *,""
      print *,"Simulated Annealing Starting"
      print *,"SA Loop: ",saConfigLive%varLoops
      print *,""
    End If
! Initial calculation
    Call runCalcs(1,optRunType)       ! opt.f90
! Store optimum (starting) nodes + rss
    splineNodesKeyOpt = splineNodesKey
    splineNodesDataOpt = splineNodesData
    optimumRSS = totalRSS   ! optimum rss - of opt nodes currently being varied
    bestRSS = totalRSS      ! rss of best found so far
! ---------------------------------------
! Run 1 - Force/Stress/Energy only
! ---------------------------------------
! Loop through and decrease temperature
    Call saOpt_VarLoop(saConfigLive)  
! Output SA results 
    Call saOpt_Output()    
  End Subroutine saOpt
! --------------------------------------------------------------------------------------------------- 
  Subroutine saOpt_VarLoop(saConfigLive)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! In vars    
    Type(saConfig) :: saConfigLive
! Private variables  
    Real(kind=DoubleReal) :: temperature, temperatureChange, temperatureBase
    Integer(kind=StandardInteger) :: i, totalLoops
    Real(kind=DoubleReal) :: aProb, randDouble
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: probTest
    Logical :: accept, bad, titlePrint
    If(TerminalPrint())Then
      print *,""
      print *,"Start Simulated Annealing - Spline Based Tabulated Function"
      print *,"Temp: ",saConfigLive%tempStart," to ",saConfigLive%tempEnd
      print *,"Variation: ",saConfigLive%maxVar," to ",saConfigLive%minVar
      print *,""
    End If    
! Base
    temperatureBase = 1.5D0
! Temperature change per variation loop
    temperatureChange = 1.0D0 / saConfigLive%varLoops
! Total loops    
    totalLoops = saConfigLive%varLoops
! Loop and vary nodes
    Do i=1,saConfigLive%varLoops
! Decrease temperature
      temperature = &
      (saConfigLive%tempStart)/((saConfigLive%tempStart/saConfigLive%tempEnd)**&
      (1.0D0*((i-1)/(1.0D0*(saConfigLive%varLoops-1)))))
! Decrease maximum vary amount
      varyAmount = &
      (saConfigLive%maxVar)/((saConfigLive%maxVar/saConfigLive%minVar)**&
      (1.0D0*((i-1)/(1.0D0*(saConfigLive%varLoops-1)))))
! Load optimum nodes    
      splineNodesKey = splineNodesKeyOpt
      splineNodesData = splineNodesDataOpt
! Vary the spline nodes on the root process
      If(mpiProcessID.eq.0)Then
        Call saOpt_VarNodes(saConfigLive%maxVar)
      End If  
! Distribute spline nodes from root to workers
      Call M_distDouble2D(splineNodesData)     
! Run calculations      
      Call runCalcs(1,optRunType)
! Acceptance - SA - Root process only
      bad = .false.
      accept = .false.
      If(mpiProcessID.eq.0)Then
        If(totalRSS.lt.optimumRSS)Then  ! improvement
          accept = .true.
        Else  
          probTest = totalRSS-bestRSS
          !totalRSS/optimumRSS
          aProb = exp((-1.0D0*(probTest))/temperature)
          randDouble = RandomLCG()
          If(randDouble.le.aProb)Then
            accept = .true.
            bad = .true.
          End If
        End If
      End If  
      Call M_distLogical(accept)
      Call M_distLogical(bad)
! Store if new value is accepted
      If(accept)Then
        optimumRSS = totalRSS
        splineNodesKeyOpt = splineNodesKey
        splineNodesDataOpt = splineNodesData
! attempt to optimise further f'(x) and f''(x)
        !Call saOpt_VarDerivs(5)        
        If(optimumRSS.lt.bestRSS)Then
          bestRSS = totalRSS  ! store best only if the best rss, not just a bac accepted result
          splineNodesKeyBest = splineNodesKeyOpt
          splineNodesDataBest = splineNodesDataOpt
        End If  
      End If
! Print out      
      If(TerminalPrint())Then
        titlePrint = .false.
        If(Mod(i,10).eq.1)Then
          titlePrint = .true.
        End If
        If(accept)Then
          If(bad)Then
            !Call printSummary(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Bad)",titlePrint)
            Call printSummaryShort(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Bad)",titlePrint)
          Else
            !Call printSummary(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Good)")
            Call printSummaryShort(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Good)",titlePrint)
          End If          
        Else 
          !Call printSummary(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"")
          Call printSummaryShort(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"",titlePrint)
        End If        
      End If       
    End Do  
  End Subroutine saOpt_VarLoop
  
    
  Subroutine saOpt_VarNodes(varyAmountMax,funcIn)
! Vary nodes  
    Implicit None   ! Force declaration of all variables
! In vars    
    Real(kind=DoubleReal) :: varyAmountMax
    Integer(kind=StandardInteger), optional :: funcIn
    Integer(kind=StandardInteger) :: func
! Private variables    
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: varyAmount
! Optional    
    func = 1
    If(Present(funcIn))Then
      func = funcIn
    End If
! Vary nodes 
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
          !varyAmount = (0.5D0-RandomLCG())*varyAmountMax
          varyAmount = (0.5D0-RandomLCG())*varyAmountMax*splineNodesData(j,func+1)
          splineNodesData(j,func+1) = splineNodesData(j,func+1)+varyAmount 
        End Do
      Else
        Exit
      End If
    End Do    
  End Subroutine saOpt_VarNodes
  
  
  Subroutine saOpt_VarDerivs(attempts)
! Vary f'(r) and f''(r) to find a better fit 
    Implicit None   ! Force declaration of all variables
! In vars    
    Integer(kind=StandardInteger) :: attempts
! Private variables    
    Integer(kind=StandardInteger) :: i    
    Logical :: betterNodes
! optimumRSS  
    Do i=1,attempts
      If(mpiProcessID.eq.0)Then
        Call saOpt_VarNodes(0.001D0,2)
        Call saOpt_VarNodes(0.001D0,3)     
      End If
! Distribute spline nodes from root to workers
      Call M_distDouble2D(splineNodesData)     
! Run calculations      
      Call runCalcs(1,optRunType)    
      betterNodes = .false.      
! Store if better
      If(totalRSS.lt.optimumRSS)Then
        optimumRSS = totalRSS
        splineNodesKeyOpt = splineNodesKey
        splineNodesDataOpt = splineNodesData  
        betterNodes = .true.     
      End If  
! Print out
      If(TerminalPrint())Then
        If(betterNodes)Then
          print *,"   var f'(x) f''(x)  ",totalRSS," *"
        Else
          print *,"   var f'(x) f''(x)  ",totalRSS
        End If
      End If
    End Do 
  End Subroutine saOpt_VarDerivs
  
  
  Subroutine saOpt_VarNode(varyAmount, nodeToVary)
! Vary nodes  
    Implicit None   ! Force declaration of all variables
! In vars    
    Real(kind=DoubleReal) :: varyAmount
    Integer(kind=StandardInteger) :: nodeToVary
! Private variables    
    Integer(kind=StandardInteger) :: i, j
! Vary nodes 
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
          If(j.eq.nodeToVary)Then
            splineNodesData(j,2) = splineNodesData(j,2)+varyAmount
          End If
        End Do
      Else
        Exit
      End If
    End Do       
  End Subroutine saOpt_VarNode
  
  
  Subroutine saOpt_Output()
! Vary f'(r) and f''(r) to find a better fit 
    Implicit None   ! Force declaration of all variables
! Optimum Spline 
    Print *,""
    Print *,""
    print *,"==================================================================="
    print *,"Optimum RSS"
    print *,"==================================================================="
! Test optimum
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt  
    Call runCalcs(1,1)
    If(TerminalPrint())Then
      print *,totalRSS,optimumRSS
    End If    
! Output Bulk Properties       
    Call outputBpT()   
    Call saveEamFile("opt_003_sa_opt.pot")
    Call eamCharts("opt_")
! Best of all tried
    !Print *,""
    !Print *,""
    !print *,"==================================================================="
    !print *,"Best RSS"
    !print *,"==================================================================="
! Test best    
    !splineNodesKey = splineNodesKeyBest
    !splineNodesData = splineNodesDataBest  
    !Call runCalcs(1,1)
    !If(TerminalPrint())Then
    !  print *,totalRSS,optimumRSS
    !End If 
! Output Bulk Properties       
    !Call outputBpT() 
    !Call saveEamFile("opt_003_sa_best.pot")
    !Call eamCharts("best_")  
  End Subroutine saOpt_Output
  
! ---------------------------------------------------------------------------------------------------  
!
! Simulated Annealing Functions - Analytic Function
!  
! --------------------------------------------------------------------------------------------------- 
  Subroutine saOpt_A(saConfigLive)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Type(saConfig) :: saConfigLive
! Start SA
    If(TerminalPrint())Then
      print *,""
      print *,"Simulated Annealing Starting"
      print *,"SA Loop: ",saConfigLive%varLoops
      print *,""
    End If
! Initial calculation
    Call runCalcs(0,optRunType,.false.)       ! runCalcs(splineNodesIn,runTypeIn,eosChartIn)
! Store optimum (starting) nodes + rss
    apfDataOpt = apfData
    !splineNodesKeyOpt = splineNodesKey
    !splineNodesDataOpt = splineNodesData
    optimumRSS = totalRSS   ! optimum rss - of opt nodes currently being varied
    bestRSS = totalRSS      ! rss of best found so far
! ---------------------------------------
! Run 1 - Force/Stress/Energy only
! ---------------------------------------
! Loop through and decrease temperature
    Call saOpt_VarLoop_A(saConfigLive)  
! Output SA results 
    Call saOpt_Output_A()    
  End Subroutine saOpt_A  
  
  ! --------------------------------------------------------------------------------------------------- 
  Subroutine saOpt_VarLoop_A(saConfigLive)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! In vars    
    Type(saConfig) :: saConfigLive
! Private variables  
    Real(kind=DoubleReal) :: temperature, temperatureChange, temperatureBase
    Integer(kind=StandardInteger) :: i, totalLoops
    Real(kind=DoubleReal) :: aProb, randDouble
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: probTest
    Logical :: accept, bad, titlePrint
    If(TerminalPrint())Then
      print *,""
      print *,"Start Simulated Annealing - Analytic Function"
      print *,"Temp: ",saConfigLive%tempStart," to ",saConfigLive%tempEnd
      print *,"Variation: ",saConfigLive%maxVar," to ",saConfigLive%minVar
      print *,""
    End If    
! Base
    temperatureBase = 1.5D0
! Temperature change per variation loop
    temperatureChange = 1.0D0 / saConfigLive%varLoops
! Total loops    
    totalLoops = saConfigLive%varLoops
! Loop and vary nodes
    Do i=1,saConfigLive%varLoops
! Decrease temperature
      temperature = &
      (saConfigLive%tempStart)/((saConfigLive%tempStart/saConfigLive%tempEnd)**&
      (1.0D0*((i-1)/(1.0D0*(saConfigLive%varLoops-1)))))
! Decrease maximum vary amount
      varyAmount = &
      (saConfigLive%maxVar)/((saConfigLive%maxVar/saConfigLive%minVar)**&
      (1.0D0*((i-1)/(1.0D0*(saConfigLive%varLoops-1)))))
! Load optimum nodes    
      splineNodesKey = splineNodesKeyOpt
      splineNodesData = splineNodesDataOpt
! Vary the spline nodes on the root process
      If(mpiProcessID.eq.0)Then
        Call varyAnalytic(apfDataOpt,varyAmount)        
      End If  
! Distribute spline nodes from root to workers
      Call M_distDouble2D(splineNodesData)     
! Run calculations      
      Call runCalcs(0,optRunType,.false.)       ! runCalcs(splineNodesIn,runTypeIn,eosChartIn)
! Acceptance - SA - Root process only
      bad = .false.
      accept = .false.
      If(mpiProcessID.eq.0)Then
        If(totalRSS.lt.optimumRSS)Then  ! improvement
          accept = .true.
        Else  
          probTest = totalRSS-bestRSS
          !totalRSS/optimumRSS
          aProb = exp((-1.0D0*(probTest))/temperature)
          randDouble = RandomLCG()
          If(randDouble.le.aProb)Then
            accept = .true.
            bad = .true.
          End If
        End If
      End If  
      Call M_distLogical(accept)
      Call M_distLogical(bad)
! Store if new value is accepted
      If(accept)Then
        optimumRSS = totalRSS
        apfDataOpt = apfData        
        If(optimumRSS.lt.bestRSS)Then
          bestRSS = totalRSS  ! store best only if the best rss, not just a bac accepted result
          splineNodesKeyBest = splineNodesKeyOpt
          splineNodesDataBest = splineNodesDataOpt
        End If  
      End If
! Print out      
      If(TerminalPrint())Then
        titlePrint = .false.
        If(Mod(i,10).eq.1)Then
          titlePrint = .true.
        End If
        If(accept)Then
          If(bad)Then
            !Call printSummary(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Bad)",titlePrint)
            Call printSummaryShort(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Bad)",titlePrint)
          Else
            !Call printSummary(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Good)")
            Call printSummaryShort(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"*(Good)",titlePrint)
          End If          
        Else 
          !Call printSummary(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"")
          Call printSummaryShort(i,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,"",titlePrint)
        End If        
      End If       
    End Do  
  End Subroutine saOpt_VarLoop_A
  
  Subroutine saOpt_Output_A()
! Vary f'(r) and f''(r) to find a better fit 
    Implicit None   ! Force declaration of all variables
! Optimum Spline 
    Print *,""
    Print *,""
    print *,"==================================================================="
    print *,"Optimum RSS"
    print *,"==================================================================="
! Test optimum
    apfData = apfDataOpt
    Call updateAnalytic()    
    Call runCalcs(0,optRunType,.false.)       ! runCalcs(splineNodesIn,runTypeIn,eosChartIn)
    If(TerminalPrint())Then
      print *,totalRSS,optimumRSS
    End If    
! Output Bulk Properties       
    Call outputBpT()   
    Call saveEamFile("opt_003_sa_opt.pot")
    Call eamCharts("opt_")
! Best of all tried
    !Print *,""
    !Print *,""
    !print *,"==================================================================="
    !print *,"Best RSS"
    !print *,"==================================================================="
! Test best    
    !splineNodesKey = splineNodesKeyBest
    !splineNodesData = splineNodesDataBest  
    !Call runCalcs(1,1)
    !If(TerminalPrint())Then
    !  print *,totalRSS,optimumRSS
    !End If 
! Output Bulk Properties       
    !Call outputBpT() 
    !Call saveEamFile("opt_003_sa_best.pot")
    !Call eamCharts("best_")  
  End Subroutine saOpt_Output_A
  
  
  
  
! ---------------------------------------------------------------------------------------------------
  Subroutine responseOptimise() 
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i
    Do i=1,4
      Call responseOptimiseProcess()
    End Do      
  End Subroutine responseOptimise
  ! ---------------------------------------------------------------------------------------------------
  Subroutine responseOptimiseProcess() 
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j, n
! make response matrix
    Call responseMatrix()  ! stored in Dimension(1:10000,1:2) :: splineNodesResponse
! Loop
    Do n=1,10
! Load spline nodes
      splineNodesKey = splineNodesKeyOpt
      splineNodesData = splineNodesDataOpt    
      If(mpiProcessID.eq.0)Then
! Vary nodes 
        Do i=1,size(splineNodesKey,1)
          If(splineNodesKey(i,1).gt.0)Then
            Do j=splineNodesKey(i,4), splineNodesKey(i,6)
              splineNodesData(j,2) = splineNodesData(j,2) * &
              (1.0D0 - 0.0001D0*splineNodesResponse(j,2)*RandomLCG())
            End Do
          Else
            Exit
          End If
        End Do        
! complete the data - y'(x) and y''(x)
        Call completeEAMSpline(splineNodesKey, splineNodesData)
      End If
! Distribute spline nodes from root to workers
      Call M_distDouble2D(splineNodesData)     
! Spline the nodes to make new eam functions        
      Call setEamSpline()     
      If(optForceZBL)Then
        Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
      End If
! Run calculations      
      Call runCalcs()
! Output       
      Call outputBpT()
! Store if better
      If(totalRSS.lt.optimumRSS)Then      
        splineNodesKeyOpt = splineNodesKey
        splineNodesDataOpt = splineNodesData
      End If
    End Do  

  End Subroutine responseOptimiseProcess
! --------------------------------------------------------------------------------------------------- 
  Subroutine responseMatrix()
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: x, y, yP, varMax
! Init    
    varMax = 0.0D0
! Reset matrix
    splineNodesResponse = 0.0D0
! build matrix
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
! Load spline nodes
          splineNodesKey = splineNodesKeyOpt
          splineNodesData = splineNodesDataOpt
! vary node by 10%      
          x = splineNodesData(j,1)    
          y = splineNodesData(j,2)
          yP = splineNodesData(j,2)*1.01D0
          splineNodesData(j,2) = yP
 ! Spline the nodes to make new eam functions        
          Call setEamSpline()     
          If(optForceZBL)Then
            Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
          End If
! Run calculations      
          Call runCalcs()    
          splineNodesResponse(j,1) = x
          splineNodesResponse(j,2) = (totalRSS-optimumRSS)/optimumRSS  ! negative is an improvement
          If(mpiProcessID.eq.0)Then
            !print *,j,x,y,yP,totalRSS,optimumRSS,(totalRSS-optimumRSS),splineNodesResponse(i,2)     
            !print *,j,x,(totalRSS-optimumRSS),splineNodesResponse(i,2) 
          End If
        End Do
      Else
        Exit
      End If
    End Do 
! Fill in gaps    
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then
        !splineNodesResponse = FillSplineResponse(splineNodesResponse,splineNodesKey(i,4),splineNodesKey(i,6))        
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
          If(varMax.lt.abs(splineNodesResponse(j,2)))Then
            varMax = abs(splineNodesResponse(j,2))
          End If
        End Do
      Else
        Exit
      End If
    End Do 
! Absolute value -1 to 1    
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
          splineNodesResponse(j,2) = splineNodesResponse(j,2) / varMax
          !If(mpiProcessID.eq.0)Then
          !  print *,j,splineNodesResponse(j,1),splineNodesResponse(j,2)
          !End If
        End Do
      Else
        Exit
      End If
    End Do     
  End Subroutine responseMatrix
! --------------------------------------------------------------------------------------------------- 
  Subroutine printSummaryShort(calcI,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,tag,titleRowIn)
    Implicit None   ! Force declaration of all variables
! In
    Integer(kind=StandardInteger) :: calcI
    Real(kind=DoubleReal) :: temperature,varyAmount,maxDensity,optimumRSS,totalRSS
    Character(*) :: tag
    Logical, Optional :: titleRowIn
! Private variables  
    Integer(kind=StandardInteger) :: configID
    Character(Len=6) :: printLine
    Logical :: titleRow
! Optional    
    titleRow = .false.
    If(Present(titleRowIn))Then
      titleRow = titleRowIn
    End If
    If(titleRow)Then
      print "(A8,A6, A14, A14, A14)", &
      " Result ",&
      " Conf ",&
      "     RSS  ",&
      "    Energy ",&
      "    Force "
      
    End If
! Int to char
    print "(A8, I6, ES14.5E3, ES14.5E3, ES14.5E3)", &
    tag,&
    calcI,&
    totalRSS,&
    rssConfigsArrTotal%energy,&
    rssConfigsArrTotal%force
    
! Output    
    !print "(A6,A7,ES14.5E3,A4,ES14.5E3,A7,A2,ES14.5E3,ES14.5E3,ES14.5E3)",&
    !adjustl(trim(printLine))," [opt: ",optimumRSS,"]  [",totalRSS,adjustl(trim(tag)),"] ",&
    !rssConfigsArrTotal%energy,rssConfigsArrTotal%force,rssBPArrTotal%aLat
  End Subroutine printSummaryShort    
  
  
! --------------------------------------------------------------------------------------------------- 
  Subroutine printSummary(calcI,temperature,varyAmount,maxDensity,optimumRSS,totalRSS,tag)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! In
    Integer(kind=StandardInteger) :: calcI
    Real(kind=DoubleReal) :: temperature,varyAmount,maxDensity,optimumRSS,totalRSS
    Character(*) :: tag
! Private variables  
    Integer(kind=StandardInteger) :: configID
    Character(Len=8) :: printLine
    
    Write(printLine,"(I8)") calcI
    print *,adjustl(trim(printLine))," [opt: ",optimumRSS,"]  [",totalRSS,adjustl(trim(tag)),"]"
    print *,"------------------------------------------------------------------------------------"
    print *,"   ",temperature,varyAmount,maxDensity
    Do configID=1,configCount
      print *,"   ","Config ",configID," E ",rssConfigsArr(configID)%energy," F ",rssConfigsArr(configID)%force,&
      " S ",rssConfigsArr(configID)%stress," T ",rssConfigsArr(configID)%total
    End Do
    Do configID=1,configCountBP 
      print *,"   ","Alat ",rssBPArr(configID)%alat,&
      " V0 ",rssBPArr(configID)%v0,&
      " E0 ",rssBPArr(configID)%e0,&
      " B0 ",rssBPArr(configID)%b0,&
      " EoS ",rssBPArr(configID)%eos,&
      " Tot ",rssBPArr(configID)%total
    End Do
    Do configID=1,configCountBP 
      print *,"   ","Alat ",calcBulkProperties(configID)%alat,&
      " V0 ",calcBulkProperties(configID)%v0,&
      " E0 ",calcBulkProperties(configID)%e0,&
      " B0 ",calcBulkProperties(configID)%b0,&
      " Bp0 ",calcBulkProperties(configID)%bp0,&      
      " C11 ",calcBulkProperties(configID)%c11,&
      " C12 ",calcBulkProperties(configID)%c12,&
      " C44 ",calcBulkProperties(configID)%c44
    End Do
    print *,""
    
  End Subroutine printSummary    
    
    
    
  
  
! ---------------------------------------------------------------------------------------------------  
!
! Starting Potential
!  
! ---------------------------------------------------------------------------------------------------     
  
  Subroutine MakeStartingPotential() 
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal) :: x, y, z, bestX, bestY, bestZ, bestRSS 
! Loop through functions/nodes
    Do i=1,4
      Do j=1,4
        Do k=1,4
          x = 0.0D0+i*0.15D0
          y = 0.0D0+j*0.15D0
          z = 0.0D0+k*0.15D0
          Call makeEAMRun(x,y,z)                          ! Make a new EAM potential
          Call setEamNodes()                         ! set spline nodes from eamKey/eamData    
          Call runCalcs(1,1) 
          print *,"rss: ",totalRSS
          If(i.eq.1.and.j.eq.1.and.k.eq.1)Then
            bestRSS = totalRSS
            bestX = x
            bestY = y
            bestZ = z
          Else
            If(totalRSS.lt.bestRSS)Then
              bestRSS = totalRSS
              bestX = x
              bestY = y
              bestZ = z              
            End If
          End If
        End Do  
      End Do  
    End Do  
    Call makeEAMRun(bestX,bestY,bestZ)                          ! Make a new EAM potential
    Call setEamNodes()                         ! set spline nodes from eamKey/eamData    
    Call runCalcs(1,1) 
    print *,"rss: ",totalRSS
  End Subroutine MakeStartingPotential 
  
! ---------------------------------------------------------------------------------------------------  
!
! Reduce Nodes
!  
! ---------------------------------------------------------------------------------------------------       
  
  Subroutine RemoveIrrelevantPoints() 
! Reduce nodes in potential
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: startRSS, trialRSS
    
    Call runCalcs(1,1) 
    startRSS = totalRSS
    Do i=1,eamFunctionCount
      Do j=splineNodesKey(i,4),splineNodesKey(i,6)  
        splineNodesData(j,2) = splineNodesData(j,2) + 1.0D-6
        Call runCalcs(1,1)
        splineNodesData(j,2) = splineNodesData(j,2) - 1.0D-6
        trialRSS = totalRSS
        print *,startRSS,trialRSS,((startRSS-trialRSS)/1.0D-6)
        !splineNodesKey
    
      End Do
    End Do
  End Subroutine RemoveIrrelevantPoints 
  
  
  
  
  
  
! ---------------------------------------------------------------------------------------------------  
!
! Noise
!  
! ---------------------------------------------------------------------------------------------------   

  Subroutine NodeNoise()
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: randomDP
! Loop through functions/nodes
    Do i=1,eamFunctionCount
      Do j=splineNodesKey(i,4),splineNodesKey(i,6)      
        If(splineNodesKey(i,3).eq.1)Then ! pair    
          If(j.ne.splineNodesKey(i,6))Then
            randomDP = RandomLCG()
            splineNodesData(j,2) = splineNodesData(j,2) * (1+(0.5D0-randomDP)/20D0)
            randomDP = RandomLCG()
            splineNodesData(j,3) = splineNodesData(j,3) * (1+(0.5D0-randomDP)/20D0)  
            randomDP = RandomLCG()
            splineNodesData(j,4) = splineNodesData(j,4) * (1+(0.5D0-randomDP)/20D0)
          End If
        End If 
        If(splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)Then ! density    
          If(j.ne.splineNodesKey(i,6))Then
            randomDP = RandomLCG()
            splineNodesData(j,2) = splineNodesData(j,2) * (1+(0.5D0-randomDP)/20D0)
            randomDP = RandomLCG()
            splineNodesData(j,3) = splineNodesData(j,3) * (1+(0.5D0-randomDP)/20D0)  
            randomDP = RandomLCG()
            splineNodesData(j,4) = splineNodesData(j,4) * (1+(0.5D0-randomDP)/20D0)
          End If
        End If 
        If(splineNodesKey(i,3).eq.3.or.splineNodesKey(i,3).eq.6.or.splineNodesKey(i,3).eq.7)Then ! embe/demb/semb
          randomDP = RandomLCG()
          splineNodesData(j,2) = splineNodesData(j,2) * (1+(0.5D0-randomDP)/20D0)
          randomDP = RandomLCG()
          splineNodesData(j,3) = splineNodesData(j,3) * (1+(0.5D0-randomDP)/20D0)  
          randomDP = RandomLCG()
          splineNodesData(j,4) = splineNodesData(j,4) * (1+(0.5D0-randomDP)/20D0)  
        End If      
      End Do
    End Do  
  End Subroutine NodeNoise 
  
  
  
  
     
    
    
    
! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

End Module opti
