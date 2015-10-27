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
      print *,""
      print *,"                           Optimise"
      print *,"----------------------------------------------------------------------"
      print *,""
      print *,""
    End If    
    
! Record input potential
    If(TerminalPrint())Then
      Print *, "Save input potential"
    End If
      
    Call saveEamFile("P_01_input.pot")    ! output.f90
    Call saveEamNodes("P_01_input.nodes") ! output.f90
    Call eamCharts("P_01_input")          ! readEAM.f90
    Call runCalcs(0,1,.true.)             ! use input potential and output EoS charts
    
! ----------------------------------
! Make an EAM to start from
! ----------------------------------    

    !Call makeEAMRun(1.0D0,1.0D0,1.0D0)                          ! Make a new EAM potential
    !Call eamCharts("Opt_1_          ")
    
    !Call setEamSpline() 
    !Call eamCharts("Opt_2_          ")
    
       
    !Call outputRssT()


    !Call MakeStartingPotential() 
    !Call saveEamFile("opt_001.pot")  
    
! ----------------------------------
! Set spline nodes, 
! ----------------------------------    
    
    Call setEamNodes()                         ! set spline nodes from eamKey/eamData   
    optEmbeddingFit = 1                        ! Force embedding functional to fit form - on    
    optDensityFit = 0                          ! Spline poly for density - off   
    Call setEamSpline(.true.)                  ! make EAM potential by splining nodes
    Call saveEamFile("P_02_start.pot")
    Call saveEamNodes("P_02_input.nodes") ! output.f90
    Call eamCharts("P_02_start")
     
! Print
    print *,"----------------------------------"
    print *,"Starting Optimisation"
    print *,"----------------------------------"
    ! ZBL + exp(poly) spline + polyspline for pair
    optEmbeddingFit = 1    ! Force embedding functional to fit form     
    optDensityFit = 0      ! Spline poly for density
    Call runCalcs(1,0)     ! Call bulk properties
    If(TerminalPrint())Then
      print *,"Spline converted input RSS ",totalRSS
    End If  
    !Call outputRssT()
    !Call outputBpT()
        
    !Call saveEamFile("opt_002.pot")    
    
! -----------------------------------------    
! Resize input functions   
! -----------------------------------------       
    
    
    
    
    
! -----------------------------------------    
! Simulated annealing    
! -----------------------------------------    
    
    ! ZBL + exp(poly) spline + polyspline for pair

    saConfigLive = saConfigIn
    Call saOpt(saConfigLive)
    
    
! Force embedding functional as a spline    
    !optEmbeddingFit = 0   
    !saConfigLive = saConfigIn
    !Call saOpt(saConfigLive)   
 
! -----------------------------------------    
! Levenberg Marquardt Algorithm
! -----------------------------------------    

    !LM_nodes = 0
    !Do i=1,eamFunctionCount
    !  If(splineNodesKey(i,3).eq.1)Then ! pair
    !    LM_nodes = LM_nodes + splineNodesKey(i,5) - 2      
    !  End If
    !  If(splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)Then ! dens/dden/sden
    !    LM_nodes = LM_nodes + splineNodesKey(i,5) - 1      
    !  End If
    !  If(splineNodesKey(i,3).eq.3.or.splineNodesKey(i,3).eq.6.or.splineNodesKey(i,3).eq.7)Then ! embe/demb/semb
    !    LM_nodes = LM_nodes + splineNodesKey(i,5)    
    !  End If
    !End Do
    !LM_nodes = 3 * LM_nodes
    
    !If(TerminalPrint())Then
    !  print *,""
    !  print *,"LMA Starting"
    !  print *,"Parameters: ",LM_nodes
    !End If  
    
    !Call LM_Opt(splineTotalNodes,LM_nodes)   
    
    
    !Call LM_Opt(crCount,splineTotalNodes)  ! must run calc once before now, to get crCount
    
    
! Set splines
    !Call setEamSpline()        ! Spline the nodes to make new eam functions
    
    !print *,crCount,splineTotalNodes
    
    !If(mpiProcessID.eq.0)Then
    !  Do i=1,crCount
    !    print *,i,calcRef(i,1),calcRef(i,2)
    !  End Do  
    !End If
    
    !If(TerminalPrint())Then
    !  print *,"Input EAM RSS value:"
    !  print *,totalRSS
    !End If
    
    
    
    
    !Call outputBpT()
    !Call saveEamFile("opt_001_input.pot")
    
! Spline functions
    !Call setEamNodes()         ! Choose the spline nodes from the full eam data for the functions
    !Call completeEAMSpline(splineNodesKey, splineNodesData) ! Fill in y'(x) and y''(x)
    !splineNodesKeyOpt = splineNodesKey
    !splineNodesDataOpt = splineNodesData

    !Call runCalcs()            ! evalEAM() + evalBulkProperties()
    !startRSS = totalRSS
    !optimumRSS = totalRSS

    
    
    
    
    
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
    !Do i=1,eamFunctionCount
    !  print *,"-----------------------------------------"
    !  print *,i,splineNodesKey(i,4),splineNodesKey(i,6)   
    !  print *,"-----------------------------------------"
   !   Do j=splineNodesKey(i,4),splineNodesKey(i,6)      
   !     print *,splineNodesData(j,1),splineNodesData(j,2),splineNodesData(j,3),splineNodesData(j,4)   
   !   End Do
   ! End Do
   ! print *,""
    
    
    Do i=1,2
      Do j=1,2
        Do k=1,2
          x = 0.0D0+i*0.3D0
          y = 0.0D0+j*0.3D0
          z = 0.0D0+k*0.3D0
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
    optRunType = 1
    Call runCalcs(1,optRunType) 
! Store optimum (starting) nodes + rss
    splineNodesKeyOpt = splineNodesKey
    splineNodesDataOpt = splineNodesData
    optimumRSS = totalRSS   ! optimum rss - of opt nodes currently being varied
    bestRSS = totalRSS      ! rss of best found so far
! ---------------------------------------
! Run 1 - Force/Stress/Energy only
! ---------------------------------------
    optRunType = 1
! Loop through and decrease temperature
    Call saOpt_VarLoop(saConfigLive)  
! Output SA results 
    Call saOpt_Output()
! -------------------------------
! Run 2 - Bulk properties + Force/Stress/Energy
! -------------------------------     
    !optRunType = 1
! Loop through and decrease temperature
    !Call saOpt_VarLoop(saConfigLive)  
! Output SA results 
    !Call saOpt_Output() 
    
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
    Logical :: accept, bad
    If(TerminalPrint())Then
      print *,""
      print *,"Start Simulated Annealing"
      print *,"Temp: ",saConfigLive%temp," to ",saConfigLive%tempEnd
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
      (saConfigLive%temp)/((saConfigLive%temp/saConfigLive%tempEnd)**&
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
        If(accept)Then
          If(bad)Then
            Call printSummary(i,temperature,varyAmount,maxDensity,&
            optimumRSS,totalRSS,"*(Bad)")
            !print *,i,temperature,varyAmount,maxDensity,": [",optimumRSS,"]  ",totalRSS,"*(Bad)"
          Else
            Call printSummary(i,temperature,varyAmount,maxDensity,&
            optimumRSS,totalRSS,"*(Good)")
            !print *,i,temperature,varyAmount,maxDensity,": [",optimumRSS,"]  ",totalRSS,"*(Good)"
          End If          
        Else 
          Call printSummary(i,temperature,varyAmount,maxDensity,&
          optimumRSS,totalRSS,"")
          !print *,i,temperature,varyAmount,maxDensity,": [",optimumRSS,"]  ",totalRSS
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
    Print *,""
    Print *,""
    print *,"==================================================================="
    print *,"Best RSS"
    print *,"==================================================================="
! Test best    
    splineNodesKey = splineNodesKeyBest
    splineNodesData = splineNodesDataBest  
    Call runCalcs(1,1)
    If(TerminalPrint())Then
      print *,totalRSS,optimumRSS
    End If 
! Output Bulk Properties       
    Call outputBpT() 
    Call saveEamFile("opt_003_sa_best.pot")
    Call eamCharts("best_")
  
  End Subroutine saOpt_Output
  
  
! ---------------------------------------------------------------------------------------------------  
!
! Levenberg-Marquardt Algorithm
!  
! --------------------------------------------------------------------------------------------------- 
  Subroutine LM_Opt(sampleCount,totalNodes)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: sampleCount,totalNodes
    Integer(kind=StandardInteger) :: f, i, n, k, m, d
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: refRSS, testRSS
    Logical :: unfixedP, calcJR
! Matrix
    Real(kind=DoubleReal), Dimension(1:sampleCount,1:totalNodes) :: J
    Real(kind=DoubleReal), Dimension(1:sampleCount) :: refVals, derivRefVals, R
!    Real(kind=DoubleReal), Dimension(1:totalNodes,1:sampleCount) :: JT     ! Transpose Jacobian
!    Real(kind=DoubleReal), Dimension(1:totalNodes,1:totalNodes) :: JTJ, JTJ_Diag  
!    Real(kind=DoubleReal), Dimension(1:totalNodes) :: JTR, P_B
    Real(kind=DoubleReal), Dimension(1:totalNodes) :: P, xP, xP_Last
    Integer(kind=StandardInteger), Dimension(1:totalNodes,1:3) :: xP_map
! LMA Dampening
    Real(kind=DoubleReal) :: lambda
! Other LMA vars
    Integer(kind=StandardInteger) :: deadRows, liveRows
    Logical :: isDead
    Character(len=500) :: fileRow
    
! Init matrices
    refVals = 0.0D0
    calcRef = 0.0D0
    derivRefVals = 0.0D0
    P = 0.0D0
    R = 0.0D0
    J = 0.0D0
    xP = 0.0D0
    xP_map = 0    
    lambda = 0.1D0
! Set x parameter matrix
    m = 0
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then  
        Do k=splineNodesKey(i,4), splineNodesKey(i,6)
          unfixedP = .true.
          If(splineNodesKey(i,3).eq.1.and.&
          (k.eq.splineNodesKey(i,4).or.k.eq.splineNodesKey(i,6)))Then
            unfixedP = .false.
          End If
          If((splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)&
          .and.k.eq.splineNodesKey(i,6))Then
            unfixedP = .false.
          End If
          If(unfixedP)Then
            Do n=2,4  ! y(x), y'(x), y''(x)
              m = m + 1
              xP(m) = splineNodesData(k,n) 
              xP_map(m,1) = k
              xP_map(m,2) = i
              xP_map(m,3) = n
            End Do  
          End If
        End Do
      End If
    End Do  
    Call runCalcs(1)
    Call calcRefRSS(testRSS)
    If(mpiProcessID.eq.0)Then
      print *,"LMA ",totalRSS, testRSS
    End If
! -----------------
! Start LMA Loop
!
    calcJR = .true.
    Do n=1,1
      If(calcJR)Then
! Build reference and residual array
        Call runCalcs(1)
        refRSS = totalRSS
        If(mpiProcessID.eq.0)Then
          print *,"ref rss: ",refRSS
        End If
        Do i=1,sampleCount
          refVals(i) = calcRef(i,2)       ! Reference values
          derivRefVals(i) = calcRef(i,1)  ! Calculated values
          R(i) = calcRef(i,1)-calcRef(i,2)
        End Do
! Calculate derivatives - vary each node 0.01%, and estimate each derivative numerically
        deadRows = 0
        Do i=1,totalNodes
          k = xP_map(i,1) ! data key
          f = xP_map(i,2) ! function
          d = xP_map(i,3) ! y(x), y'(x), y''(x)
          varyAmount = abs(0.01D0*xP(i))
          If(varyAmount.eq.0.0D0)Then
            varyAmount = 0.0001D0
          End If
          splineNodesData(k,d) = xP(i)+varyAmount ! Perturb node
          !Call CompleteNodeData(splineNodesData, splineNodesKey(f,4), splineNodesKey(f,6))
          Call runCalcs(1)
          Call calcRefRSS(testRSS)
          print *,i,totalRSS, testRSS
          splineNodesData(k,d) = xP(i)            ! reset node
          !print *,i," node key ",k,f,d,totalRSS,"          ",xP(i),varyAmount
          Do m=1,500
            fileRow(m:m) = " "
          End Do
          isDead = .true.
          Do m=1,sampleCount
            J(m,i) = (calcRef(m,1) - derivRefVals(m))/varyAmount  ! Store calculated deriv straight into J matrix
            If(J(m,i).ne.0.0D0)Then
              isDead = .false.
            End If
          End Do   
          If(isDead)Then
            deadRows = deadRows + 1
          End If          
        End Do
        print *,"Dead rows: ",deadRows
        liveRows = totalNodes - deadRows
      End If  
!-----------------------------------      
! calculate change matrix
!----------------------------------- 
      Call LM_OptCalc(sampleCount,totalNodes,liveRows,lambda,J,R,P) 
! store last set of points
      xP_Last = xP 
! update with change matrix P
      Do i=1,totalNodes
        xP(i) = xP(i) + P(i)
        If(mpiProcessID.eq.1)Then
          print *,i,P(i)
        End If  
      End Do
! Store last spline points
      splineNodesKeyTemp = splineNodesKey   
      splineNodesDataTemp = splineNodesData  
      Call runCalcs(1)
      Call calcRefRSS(testRSS)
      print *,"A",totalRSS, testRSS
! Update spline points
      Call LM_UpdateNodes(P, xP_map)
! Calculate RSS      
      Call runCalcs(1)
      Call calcRefRSS(testRSS)
      print *,"B",totalRSS, testRSS
      
      If(mpiProcessID.eq.1)Then
        Do i=1,28
          print *,i,splineNodesDataTemp(i,1),splineNodesDataTemp(i,2),splineNodesDataTemp(i,3),splineNodesDataTemp(i,4)
        End Do
        print *,""
        print *,""
        Do i=1,28
          print *,i,splineNodesData(i,1),splineNodesData(i,2),splineNodesData(i,3),splineNodesData(i,4)
        End Do
      End If
      
!
! End LMA Loop
! -----------------     
    End Do
    
    
  End Subroutine LM_Opt
  
  Subroutine LM_UpdateNodes(P, xP_map)
! Updates the spline nodes with the new parameters calculated by LMA iteration
    Implicit None   ! Force declaration of all variables
! Private variables  
    Real(kind=DoubleReal), Dimension(:) :: P
    Integer(kind=StandardInteger), Dimension(:,:) :: xP_map
    Integer(kind=StandardInteger) :: i, k, d
! Update
    Do i=1,size(xP_map,1)
        If(Isnan(P(i)))Then
          ! skip
        Else
          k = xP_map(i,1)
          d = xP_map(i,3) 
          splineNodesData(k,d) = splineNodesData(k,d)+P(i)
        End If  
      End Do    
  End Subroutine LM_UpdateNodes
  
  Subroutine LM_UpdateNodesR(P, xP_map)
! Updates the spline nodes with the new parameters calculated by LMA iteration (Reverse)
    Implicit None   ! Force declaration of all variables
! Private variables  
    Real(kind=DoubleReal), Dimension(:) :: P
    Integer(kind=StandardInteger), Dimension(:,:) :: xP_map
    Integer(kind=StandardInteger) :: i, k, d
! Update
    Do i=1,size(xP_map,1)
        If(Isnan(P(i)))Then
          ! skip
        Else
          k = xP_map(i,1)
          d = xP_map(i,3) 
          splineNodesData(k,d) = splineNodesData(k,d)-P(i)
        End If  
      End Do    
  End Subroutine LM_UpdateNodesR
  
  
  Subroutine LM_OptCalc(sampleCount,totalNodes,liveRows,lambda,J,R,P)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
! Input integers + dps
    Integer(kind=StandardInteger) :: sampleCount, totalNodes, liveRows    
    Real(kind=DoubleReal) :: lambda
! Input matrices
    Real(kind=DoubleReal), Dimension(:,:) :: J
    Real(kind=DoubleReal), Dimension(:) :: R
    Real(kind=DoubleReal), Dimension(:) :: P
! Working matrices    
    Real(kind=DoubleReal), Dimension(1:sampleCount,1:liveRows) :: Jw
    Real(kind=DoubleReal), Dimension(1:sampleCount) :: Rw
    Real(kind=DoubleReal), Dimension(1:liveRows) :: Pw
! Matrices
    Real(kind=DoubleReal), Dimension(1:liveRows,1:sampleCount) :: JT
    Real(kind=DoubleReal), Dimension(1:liveRows,1:liveRows) :: JTJ, JTJ_Diag
    Real(kind=DoubleReal), Dimension(1:liveRows) :: JTR
! Other vars
    Integer(kind=StandardInteger), Dimension(1:liveRows) :: liveToTotalMap
    Logical :: isDead
    Integer(kind=StandardInteger) :: i, n, k
    
! Make working matrices - J
    k = 0
    Do i=1,totalNodes  
      isDead = .true.
      Do n=1,sampleCount
        If(J(n,i).ne.0.0D0)Then
          isDead = .false.
        End If
      End Do
      If(isDead)Then
        ! Skip
      Else
        k = k + 1
        liveToTotalMap(k) = i ! i = key in P, k = key in Pw
        Do n=1,sampleCount
          Jw(n,k) = J(n,i)
        End Do
        !If(mpiProcessID.eq.0)Then
        !  print *,k,i,liveRows,totalNodes
        !End If  
      End If
    End Do
    
! Rw is the same as R
    Rw = R
    
!***********     
! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
!***********      
    JT = TransposeMatrix(Jw)
    JTJ = matmul(JT,Jw)   
    JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
    JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ  
    JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
    Pw = matmul(JTJ,JTR) 
! Map to P matrix
    P = 0.0D0
    Do i=1,liveRows
      P(liveToTotalMap(i)) = Pw(i)
    End Do
    
  End Subroutine LM_OptCalc
  
  
  
  
  
  
  
  
  Subroutine LM_OptC(sampleCount,totalNodes)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: sampleCount,totalNodes
    Integer(kind=StandardInteger) :: f, i, n, k, m, d
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: refRSS
    Logical :: unfixedP, calcJR
! Matrix
    Real(kind=DoubleReal), Dimension(1:sampleCount,1:totalNodes) :: J
    Real(kind=DoubleReal), Dimension(1:sampleCount) :: refVals, derivRefVals, R
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:sampleCount) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:totalNodes) :: JTJ, JTJ_Diag  
    Real(kind=DoubleReal), Dimension(1:totalNodes) :: JTR, P, P_B, xP, xP_Last
    Integer(kind=StandardInteger), Dimension(1:totalNodes,1:3) :: xP_map
! LMA Dampening
    Real(kind=DoubleReal) :: lambda
    
    Character(len=14) :: tempStr
    Character(len=500) :: fileRow
    
! Init matrices
    refVals = 0.0D0
    calcRef = 0.0D0
    derivRefVals = 0.0D0
    P = 0.0D0
    R = 0.0D0
    J = 0.0D0
    xP = 0.0D0
    xP_map = 0    
    lambda = 0.1D0
! Set x parameter matrix
    m = 0
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then  
        Do k=splineNodesKey(i,4), splineNodesKey(i,6)
          unfixedP = .true.
          If(splineNodesKey(i,3).eq.1.and.&
          (k.eq.splineNodesKey(i,4).or.k.eq.splineNodesKey(i,6)))Then
            unfixedP = .false.
          End If
          If((splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)&
          .and.k.eq.splineNodesKey(i,6))Then
            unfixedP = .false.
          End If
          If(unfixedP)Then
            Do n=2,4  ! y(x), y'(x), y''(x)
              m = m + 1
              xP(m) = splineNodesData(k,n) 
              xP_map(m,1) = k
              xP_map(m,2) = i
              xP_map(m,3) = n
            End Do  
          End If
        End Do
      End If
    End Do  
    Call runCalcs(1)
    If(mpiProcessID.eq.0)Then
      print *,"LMA ",totalRSS
    End If
! -----------------
! Start LMA Loop
!
    calcJR = .true.
    Do n=1,1
      If(calcJR)Then
! Build reference and residual array
        Call runCalcs(1)
        refRSS = totalRSS
        If(mpiProcessID.eq.0)Then
          print *,"ref rss: ",refRSS
        End If
        Do i=1,sampleCount
          refVals(i) = calcRef(i,2)       ! Reference values
          derivRefVals(i) = calcRef(i,1)  ! Calculated values
          R(i) = calcRef(i,1)-calcRef(i,2)
        End Do
! Calculate derivatives - vary each node 0.1%, and estimate each derivative numerically
        Do i=1,totalNodes
          k = xP_map(i,1) ! data key
          f = xP_map(i,2) ! function
          d = xP_map(i,3) ! y(x), y'(x), y''(x)
          varyAmount = abs(0.001D0*xP(i))
          If(varyAmount.eq.0.0D0)Then
            varyAmount = 0.00001D0
          End If
          splineNodesData(k,d) = xP(i)+varyAmount ! Perturb node
          !Call CompleteNodeData(splineNodesData, splineNodesKey(f,4), splineNodesKey(f,6))
          Call runCalcs(1)
          print *,i,totalRSS
          splineNodesData(k,d) = xP(i)            ! reset node
          !print *,i," node key ",k,f,d,totalRSS,"          ",xP(i),varyAmount
          Do m=1,500
            fileRow(m:m) = " "
          End Do
          Do m=1,sampleCount
            J(m,i) = (calcRef(m,1) - derivRefVals(m))/varyAmount  ! Store calculated deriv straight into J matrix
            write(tempStr, "(E12.5)") J(m,i)
            fileRow = trim(fileRow)//" "//trim(adjustl(tempStr))
          End Do      
          !print *,i,fileRow    
        End Do
      End If  
!-----------------------------------      
! calculate change matrix
!----------------------------------- 
      !***********     
      ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
      !***********      
! Transpose Jacobian
      !JT = TransposeMatrix(J)
      !JTJ = matmul(JT,J)
      !JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      !JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ  
      !JTR = matmul(JT,R)
      !JTR = -1.0D0*JTR ! Recycle JTR var
      !P_B = SolveLinearSet(JTJ,JTR)
      
      JT = TransposeMatrix(J)
      JTJ = matmul(JT,J)   
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ  
      JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
      P = matmul(JTJ,JTR) 
      
      If(mpiProcessID.eq.0)Then
        Do i=1,totalNodes
          print *,i,P(i),P_B(i)
        End Do
      End If
! Store last loop values
      !lastRSS = refRSS 
      !J_Last = J
      !R_Last = R      
      xP_Last = xP      
! Store last spline points
      splineNodesKeyTemp = splineNodesKey   
      splineNodesDataTemp = splineNodesData         
! Update parameters      
      Do i=1,totalNodes
        If(Isnan(P(i)))Then
          ! skip
        Else
          xP = xP + P(i)
          k = xP_map(i,1)
          d = xP_map(i,3) 
          splineNodesData(k,d) = splineNodesData(k,d)+P(i)
          print *,i,P(i)
        End If  
      End Do      
! Run to calc rss
      Call runCalcs(1) ! run calc.
      print *, n, lambda, totalRSS      
!-----------------------------------      
! adjust lambda parameter
!-----------------------------------    
      If(totalRSS.lt.refRSS)Then
        lambda = lambda * 0.2D0
        calcJR = .true.
      Else
        lambda = lambda * 1.5D0
        calcJR = .false.
        xp = xP_Last
        splineNodesKey = splineNodesKeyTemp
        splineNodesData = splineNodesDataTemp
      End If
!
! End LMA Loop
! -----------------     
    End Do
    
    
  End Subroutine LM_OptC
  
  
  
  Subroutine LM_OptOldB(sampleCount,totalNodes)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: sampleCount,totalNodes
    Integer(kind=StandardInteger) :: f, i, n, k, m
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: refRSS
    Logical :: unfixedP, calcJR
! Matrix
    Real(kind=DoubleReal), Dimension(1:sampleCount,1:totalNodes) :: J
    Real(kind=DoubleReal), Dimension(1:sampleCount) :: refVals, derivRefVals, R
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:sampleCount) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:totalNodes) :: JTJ, JTJ_Diag  
    Real(kind=DoubleReal), Dimension(1:totalNodes) :: JTR, P, P_B, xP, xP_Last
    Integer(kind=StandardInteger), Dimension(1:totalNodes,1:2) :: xP_map
! LMA Dampening
    Real(kind=DoubleReal) :: lambda
    
! Init matrices
    refVals = 0.0D0
    calcRef = 0.0D0
    derivRefVals = 0.0D0
    P = 0.0D0
    R = 0.0D0
    J = 0.0D0
    xP = 0.0D0
    xP_map = 0    
    lambda = 0.1D0
! Set x parameter matrix    
    m = 0
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then  
        Do k=splineNodesKey(i,4), splineNodesKey(i,6)
          unfixedP = .true.
          If(splineNodesKey(i,3).eq.1.and.&
          (k.eq.splineNodesKey(i,4).or.k.eq.splineNodesKey(i,6)))Then
            unfixedP = .false.
          End If
          If((splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)&
          .and.k.eq.splineNodesKey(i,6))Then
            unfixedP = .false.
          End If
          If(unfixedP)Then
            m = m + 1
            xP(m) = splineNodesData(k,2) 
            xP_map(m,1) = k
            xP_map(m,2) = k
          End If
        End Do
      End If
    End Do  
    Call runCalcs(1)
    If(mpiProcessID.eq.0)Then
      print *,"LMA ",totalRSS
    End If
! -----------------
! Start LMA Loop
!
    calcJR = .true.
    Do n=1,1
      If(calcJR)Then
! Build reference and residual array
        Call runCalcs(1)
        refRSS = totalRSS
        Do i=1,sampleCount
          refVals(i) = calcRef(i,2)       ! Reference values
          derivRefVals(i) = calcRef(i,1)  ! Calculated values
          R(i) = calcRef(i,1)-calcRef(i,2)
        End Do
! Calculate derivatives - vary each node 0.1%, and estimate each derivative numerically
        Do i=1,totalNodes
          k = xP_map(i,1) ! data key
          f = xP_map(i,2) ! function
          varyAmount = abs(0.001D0*xP(i))
          If(varyAmount.eq.0.0D0)Then
            varyAmount = 0.00001D0
          End If
          splineNodesData(k,2) = xP(i)+varyAmount ! Perturb node
          Call CompleteNodeData(splineNodesData, splineNodesKey(f,4), splineNodesKey(f,6))
          Call runCalcs(1)
          splineNodesData(k,2) = xP(i)            ! reset node
          print *,i," node key ",k,f,totalRSS,"          ",xP(i),varyAmount
          Do m=1,sampleCount
            J(m,i) = (calcRef(m,1) - derivRefVals(m))/varyAmount  ! Store calculated deriv straight into J matrix
          End Do          
        End Do
      End If  
!-----------------------------------      
! calculate change matrix
!----------------------------------- 
      !***********     
      ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
      !***********      
! Transpose Jacobian
      JT = TransposeMatrix(J)
      JTJ = matmul(JT,J)
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ  
      JTR = matmul(JT,R)
      JTR = -1.0D0*JTR ! Recycle JTR var
      P = SolveLinearSet(JTJ,JTR)
      
      JT = TransposeMatrix(J)
      JTJ = matmul(JT,J)   
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ  
      JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
      P_B = matmul(JTJ,JTR) 
      
      If(mpiProcessID.eq.0)Then
        Do i=1,totalNodes
          print *,i,P(i),P_B(i)
        End Do
      End If
! Store last loop values
      !lastRSS = refRSS 
      !J_Last = J
      !R_Last = R      
      xP_Last = xP      
! Store last spline points
      !splineNodesKeyTemp = splineNodesKey   
      !splineNodesDataTemp = splineNodesData         
! Update parameters      
      Do i=1,totalNodes
        xP = xP + P(i)
        k = xP_map(i,1)
        !splineNodesData(k,2) = splineNodesData(k,2)+P(i)
      End Do      
! Run to calc rss
      Call runCalcs(1) ! run calc.
      print *, n, lambda, totalRSS      
!-----------------------------------      
! adjust lambda parameter
!-----------------------------------    
      If(totalRSS.lt.refRSS)Then
        lambda = lambda * 0.2D0
        calcJR = .true.
      Else
        lambda = lambda * 1.5D0
        calcJR = .false.
        xp = xP_Last
        !splineNodesKey = splineNodesKeyTemp
        !splineNodesData = splineNodesDataTemp
      End If
!
! End LMA Loop
! -----------------     
    End Do
    
    
  End Subroutine LM_OptOldB
  
  Subroutine LM_OptOld(sampleCount,totalNodes)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: sampleCount,totalNodes
    Integer(kind=StandardInteger) :: i, n, k, m
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: refRSS, lastRSS
! Matrix
    Real(kind=DoubleReal), Dimension(1:sampleCount,1:totalNodes) :: J, J_Last
    Real(kind=DoubleReal), Dimension(1:sampleCount) :: refVals, derivRefVals, R, R_Last
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:sampleCount) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:totalNodes) :: JTJ, JTJ_Diag  
    Real(kind=DoubleReal), Dimension(1:totalNodes) :: JTR, P
    Character(len=2) :: numberStr
    !Character(len=255) :: row
! LMA Dampening
    Real(kind=DoubleReal) :: lambda
    
! Init matrices
    refVals = 0.0D0
    calcRef = 0.0D0
    derivRefVals = 0.0D0
    R = 0.0D0
    J = 0.0D0
    
    
    lambda = 10.0D0
    Do n=1,4
      write(numberStr,"(I2)"), n
! Build reference and residual array
      Call runCalcs(1)
      refRSS = totalRSS
      Call saveEamFile("lm_"//trim(adjustl(numberStr))//".pot")
      Do i=1,sampleCount
        refVals(i) = calcRef(i,2)       ! Reference values
        derivRefVals(i) = calcRef(i,1)  ! Calculated values
        R(i) = calcRef(i,1)-calcRef(i,2)
      End Do
      print *,n,lambda,refRSS
! Calculate derivatives - vary each node 0.1%, and estimate each derivative numerically
      Do i=1,size(splineNodesKey,1)
        If(splineNodesKey(i,1).gt.0)Then    
          Do k=splineNodesKey(i,4), splineNodesKey(i,6)
            varyAmount = abs(0.001D0*splineNodesData(k,2))
            If(varyAmount.eq.0.0D0)Then
              varyAmount = 0.0001D0
            End If
            splineNodesData(k,2) = splineNodesData(k,2)+varyAmount ! vary node
            Call CompleteNodeData(splineNodesData, splineNodesKey(i,4), splineNodesKey(i,6))  ! Fill in the y'(x) and y''(x) values
            Call runCalcs(1) ! run calc
            Do m=1,sampleCount
              J(m,k) = (calcRef(m,1) - derivRefVals(m))/varyAmount  ! Store calculated deriv straight into J matrix
              !print *,m,k,calcRef(m,1), derivRefVals(m), (calcRef(m,1) - derivRefVals(m)), varyAmount, J(m,k)
            End Do
            splineNodesData(k,2) = splineNodesData(k,2)-varyAmount !reset node
          End Do
        Else
          Exit  ! Exit, all functions cycled through
        End If          
      End Do
! Choose whether to accept update or increase/decrease lambda      
      If(n.gt.1)Then
        If(refRSS.gt.lastRSS)Then  ! If worse...reject, and increase lambda
! Discard changes, increase lambda
          J = J_Last
          R = R_Last
          splineNodesKey = splineNodesKeyTemp
          splineNodesData = splineNodesDataTemp
          refRSS = lastRSS
          lastRSS = -1.0D0
          lambda = lambda * 1.5D0
        End If
        If(refRSS.lt.lastRSS)Then  ! If better...accept, and decrease lambda
          lambda = lambda * 0.2D0
        End If
      End If
! calculate change matrix
      !***********     
      ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
      !***********      
! Transpose Jacobian
      JT = TransposeMatrix(J)
      JTJ = matmul(JT,J)
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
      JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
      JTR = matmul(JT,R)
      JTR = -1.0D0*JTR ! Recycle JTR var
      P = matmul(JTJ,JTR) 
! Store last loop values
      lastRSS = refRSS 
      J_Last = J
      R_Last = R   
! Store last spline points
      splineNodesKeyTemp = splineNodesKey   
      splineNodesDataTemp = splineNodesData           
! Update parameters      
      Do i=1,totalNodes
        splineNodesData(i,2) = splineNodesData(i,2)+P(i)
      End Do
      Do i=1,size(splineNodesKey,1)
        If(splineNodesKey(i,1).gt.0)Then    
          Call CompleteNodeData(splineNodesData, splineNodesKey(i,4), splineNodesKey(i,6))
        Else
          Exit  ! Exit, all functions cycled through
        End If
      End Do
! Run to calc rss
      Call runCalcs(1) ! run calc.
      print *, totalRSS
      !Do i=1,totalNodes 
      !  print *,i,": ",P(i)
      !End Do
    End Do
    
    
  End Subroutine LM_OptOld
  
  
  
  
  Subroutine LM_OptProcess(totalNodes)
! LMA
! 1 err value, totalNodes = parameters
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, n, k, totalNodes
    Real(kind=DoubleReal), Dimension(1:1,1:totalNodes) :: J, J_Last
    Real(kind=DoubleReal), Dimension(1:1) :: R, R_Last
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:1) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:totalNodes) :: JTJ, JTJ_Diag  
    Real(kind=DoubleReal), Dimension(1:totalNodes) :: JTR, P
    Real(kind=DoubleReal) :: varyAmount, grad
    Real(kind=DoubleReal) :: refRSS, lastRSS
    Real(kind=DoubleReal) :: lambda
! Load optimum nodes    
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt
    
    If(mpiProcessID.eq.0)Then
      print *,""
      print *,"LMA"
      print *,""
    End If

    lambda = 1.0D0
    Do n=1,3
      Call runCalcs()
      refRSS = totalRSS
      If(mpiProcessID.eq.0)Then
        print *,n,"   ",refRSS
      End If
! Build R
      R(1) = refRSS
      Do i=1,size(splineNodesKey,1)
        If(splineNodesKey(i,1).gt.0)Then    
          Do k=splineNodesKey(i,4), splineNodesKey(i,6)
! Reload optimum nodes    
            splineNodesKey = splineNodesKeyOpt
            splineNodesData = splineNodesDataOpt
! Perturb node by e=0.1% on root process
            varyAmount = abs(0.001D0*splineNodesData(k,2))
            splineNodesData(k,2) = splineNodesData(k,2)+varyAmount

            Call CompleteNodeData(splineNodesData, splineNodesKey(i,4), splineNodesKey(i,6))  ! Fill in the y'(x) and y''(x) values
            Call setEamSpline()        ! Spline the nodes to make new eam functions
            If(optForceZBL)Then
              Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
            End If 
! Run calculation        
            Call runCalcs()
! Gradient        
            grad = (totalRSS-refRSS)/varyAmount
! Build J            
            J(1,k) = grad
          End Do
        Else
          Exit  ! Exit, all functions cycled through
        End If
      End Do
! Choose whether to accept update or increase/decrease lambda      
      If(n.gt.1)Then
        print *,n,refRSS,lastRSS
! Delayed gratification scheme - 1.5*lambda or 0.2*lambda     
        If(refRSS.gt.lastRSS)Then  ! If worse...reject, and increase lambda
! Discard changes, increase lambda
          J = J_Last
          R = R_Last
          splineNodesKey = splineNodesKeyTemp
          splineNodesData = splineNodesDataTemp
          Call setEamSpline()        ! Spline the nodes to make new eam functions
          If(optForceZBL)Then
            Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
          End If
          refRSS = lastRSS
          lastRSS = -1.0D0
          lambda = lambda * 1.5D0
        End If
        If(refRSS.lt.lastRSS)Then  ! If better...accept, and decrease lambda
          lambda = lambda * 0.2D0
        End If
      End If
! calculate change matrix
      !***********     
      ! P = (JTJ+L*diag(JTJ))^(-1)(-1*JTR)   
      !***********      
! Transpose Jacobian
      JT = TransposeMatrix(J)
      JTJ = matmul(JT,J)
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
      JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
      JTR = matmul(JT,R)
      JTR = -1.0D0*JTR ! Recycle JTR var
      P = matmul(JTJ,JTR)  
! Store last loop values
      lastRSS = refRSS 
      J_Last = J
      R_Last = R   
! Store last spline points
      splineNodesKeyTemp = splineNodesKey   
      splineNodesDataTemp = splineNodesData           
! Update parameters      
      Do i=1,totalNodes
        splineNodesData(i,2) = splineNodesData(i,2)+P(i)
      End Do
      If(mpiProcessID.eq.0)Then
        print *,P(1),P(2),P(3),P(4),P(5),P(6),P(7),P(8),P(30),P(40)
      End If
      Do i=1,size(splineNodesKey,1)
        If(splineNodesKey(i,1).gt.0)Then    
          Call CompleteNodeData(splineNodesData, splineNodesKey(i,4), splineNodesKey(i,6))
        Else
          Exit  ! Exit, all functions cycled through
        End If
      End Do
      Call setEamSpline()        ! Spline the nodes to make new eam functions
      If(optForceZBL)Then
        Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
      End If
    End Do
    
  
  
    
  End Subroutine LM_OptProcess
  
  
  
  
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
      " Bp0 ",calcBulkProperties(configID)%bp0
    End Do
    print *,""
    
    !rssBPArr(configID)%total = rssBPArr(configID)%alat+rssBPArr(configID)%v0+&
    !rssBPArr(configID)%e0+rssBPArr(configID)%b0+rssBPArr(configID)%bp0+&
    !rssBPArr(configID)%c11+rssBPArr(configID)%c12+rssBPArr(configID)%c44+&
    !rssBPArr(configID)%eos    

  End Subroutine printSummary    
    
    
! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

End Module opti
