Module optimise

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
  Use readEAM
  Use calcEAM 
  Use testEAM 
  Use calcEval  
  Use bulkProperties 

! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private    
! Public Subroutines
  Public :: runOptimise
  
Contains  



  Subroutine runOptimise()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Real(kind=DoubleReal) :: saTempIn, saMaxVariationIn
    Real(kind=DoubleReal), Dimension(1:20) :: rssWeightingIn
    Integer(kind=StandardInteger) :: saTempLoopsIn, saVarLoopsIn
    Integer(kind=StandardInteger) :: optLoopsIn
    Character(len=4) :: eampaRunTypeIn
! Store input variables
    saTempIn = saTemp
    saMaxVariationIn = saMaxVariation
    saTempLoopsIn = saTempLoops
    rssWeightingIn = rssWeighting
    saVarLoopsIn = saVarLoops
    optLoopsIn = optLoops
    eampaRunTypeIn = eampaRunType
! Choose which type of optimisation process
    If(eampaRunType(1:4).eq."OPTI".or.eampaRunType(1:4).eq."OPTF".or.&
    eampaRunType(1:4).eq."OPTT")Then
!----------------------------------------------------
! Config only, Full or Test only optimisation    
!----------------------------------------------------
      Call runOptimiseProcess()
    End If
    If(eampaRunType(1:4).eq."OPTE")Then
!----------------------------------------------------
! Extensive optimisation       
!----------------------------------------------------
! Start optimising potential using test EoS only, force embedding fit to 6th order poly
      eampaRunType = "OPES"
      forceEmbeFitOpt = 1
      saTemp = 50.0
      saMaxVariation = 0.1D0
      saTempLoops = 5
      saVarLoops = 50
      optLoops = 5
      rssWeighting = 0.0D0
      rssWeighting(8) = 1.0D0      ! e0s
      Call runOptimiseProcess(1)   
! Start optimising potential using test (BM, EoS etc) only, force embedding fit to 6th order poly     
      eampaRunType = "OPTT"
      forceEmbeFitOpt = 1
      saTemp = 5.0
      saMaxVariation = 0.02D0
      saTempLoops = 3
      saVarLoops = 50
      optLoops = 5
      rssWeighting = 0.0D0
      rssWeighting(4) = 100.0D0   ! alat
      rssWeighting(5) = 100.0D0   ! emin
      rssWeighting(6) = 1.0D0      ! bm
      rssWeighting(7) = 3.0D0      ! ecs
      rssWeighting(8) = 1.0D0      ! e0s
      Call runOptimiseProcess(2)        
! Start optimising potential using test (BM, EoS etc) only, force embedding fit to 6th order poly     
      eampaRunType = "OPTT"
      forceEmbeFitOpt = 0
      saTemp = 1.0
      saMaxVariation = 0.02D0
      saTempLoops = 0
      saVarLoops = 0
      optLoops = 10
      rssWeighting = 0.0D0
      rssWeighting(4) = 100.0D0   ! alat
      rssWeighting(5) = 100.0D0   ! emin
      rssWeighting(6) = 1.0D0      ! bm
      rssWeighting(7) = 3.0D0      ! ecs
      rssWeighting(8) = 1.0D0      ! e0s
      Call runOptimiseProcess(2)          
      
      
      !rssWeighting(4) = 1000.0D0   ! alat
      !rssWeighting(5) = 1000.0D0   ! emin
      !rssWeighting(6) = 1.0D0      ! bm
      !rssWeighting(7) = 0.2D0      ! ecs
! Reload input settings
      
      
    End If
    
  End Subroutine runOptimise
  
  Subroutine runOptimiseProcess(loadInputIn)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Character(len=255) :: eamTempPotPath
    Character(len=64) :: fileName
    Character(len=32) :: fileNumber
    Real(kind=DoubleReal) :: timeStartOpt, timeEndOpt
    Real(kind=DoubleReal) :: inputRSS
    Integer(kind=StandardInteger) :: saveForcesToFileTemp
    Integer(kind=StandardInteger) :: loadInput
    Integer(kind=StandardInteger), optional :: loadInputIn
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Total Spline Nodes: ",splineTotalNodes
    End If
! Start Time
    Call cpu_time(timeStartOpt)
! Optional arguments    
    loadInput = 1
    If(Present(loadInputIn))Then
      loadInput = loadInputIn
    End If
! Init variables
    fileNumber = BlankString(fileNumber)  
! Forces to file
    saveForcesToFileTemp = saveForcesToFile
    saveForcesToFile = 0
! Print out to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Start Optimise EAM"
    End If
! Make directory for temp potentials to debug (if required)
    eamTempPotPath = Trim(outputDirectory)//"/pots"
    Call makeDir(eamTempPotPath)
! Load original input potential    
    If(loadInput.eq.1)Then
      eamKey = eamKeyInput
      eamData = eamDataInput
! Get starting nodes    
      Call setEamNodes()
! Rescale Embedding Function    
      Call scaleEmbedding()
! Fix nodes near ZBL join    
      If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then    
        Call fixZBLNodes()
      End If
! Save starting nodes as optimum    
      splineNodesKeyOpt = splineNodesKey
      splineNodesDataOpt = splineNodesData   
! Force ZBL if required
      If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
        Call eamZblHardCore()
      End If  
    End If
! Load optimum potential    
    If(loadInput.eq.2)Then    
      Call forceEmbeddingFit()                              ! Force embedding function to fit 6th order poly
      Call setEamSpline()                                   ! Make spline EAM from nodes
      If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then  ! Apply ZBL hard core if required    
        Call eamZblHardCore()
      End If 
    End If
!   
!    
!-------------------------------------  
! Evaluate Input EAM    
!-------------------------------------
    fileName = "potInput.pot"
    Call saveEamFile(fileName)
    Call evaluate()    
    inputRSS = totalRSS
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"Input RSS: ",totalRSS,"*"
    End If  
! Test input potential
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,""
      Print *,""
      Print *,""
      Print *,"**********************************************************************"
      Print *,"Input Potential Test Results"
      Print *,"**********************************************************************"
    End If 
    Call runTestAnalysis()  
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"**********************************************************************"
    End If      
!
!  
!-------------------------------------      
! Spline EAM - Evaluate Splined EAM
!-------------------------------------
    Call forceEmbeddingFit()                                ! Force embedding function to fit 6th order poly
    Call setEamSpline()                                     ! Make EAM functions spline of nodes
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then    ! Apply the zbl hard core
      Call eamZblHardCore()
    End If  
    fileName = "potSpline.pot"  
    Call saveEamFile(fileName) 
    Call evaluate()         
    splineNodesKeyOpt = splineNodesKey
    splineNodesDataOpt = splineNodesData   
    optimumRSS = totalRSS  
    startRSS = totalRSS
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,""
      Print *,""
      Print *,""
      Print *,"**********************************************************************"
      Print *,"Splined Input Potential Test Results"
      Print *,"**********************************************************************"
    End If 
    Call runTestAnalysis()  
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"**********************************************************************"
    End If  
    
!
!  
!-------------------------------------      
! Reduce nodes
!-------------------------------------    
! Response of each node to a small perturbation 
    If(reduceNodes.eq.1)Then
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
        print *,"-------------------------------------------------"
        print *,"Perturbation response and reduce nodes"
        print *,"-------------------------------------------------"
      End If
      Call nodeResponseList()
! Reduce the number of nodes
      Call reduceNodeList()
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
        print *,"Reduced Spline Nodes: ",splineTotalNodes
      End If
      Call setEamSpline()
      If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
        Call eamZblHardCore()
      End If      
      fileName = "reducedPotSpline.pot"  
      Call saveEamFile(fileName) 
    End If    
!
!  
!---------------------------------
! Vary/Jumble all nodes
!---------------------------------    
! Jumble nodes
    If(jumbleNodesOpt.eq.1)Then
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
        print *,"-------------------------------------------------"
        print *,"Jumble starting nodes ",ProgramTime ()
        print *,"-------------------------------------------------"
      End If
      Call jumbleNodes()
! Store jumbled nodes as optimum (starting point)
      Call evaluate()         
      splineNodesKeyOpt = splineNodesKey
      splineNodesDataOpt = splineNodesData   
      optimumRSS = totalRSS
      startRSS = totalRSS
    End If  
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"Starting RSS: ",startRSS,"*"
    End If  
    If(saTempLoops.gt.0.and.saVarLoops.gt.0)Then
! Save potential    
      fileName = "preSA.pot"  
      Call saveEamFile(fileName) 
!
!  
!---------------------------------
! Simulated Annealing
!---------------------------------   
! Run simulated annealing
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
        print *,"-------------------------------------------------"
        print *,"Simulated annealing ",ProgramTime ()
        print *,"Temperature: ",saTemp
        print *,"T Loops: ",saTempLoops
        print *,"V Loops: ",saVarLoops
        print *,"Total Vars: ",(saVarLoops*saTempLoops)
        print *,"-------------------------------------------------"
      End If    
      Call simulatedAnnealing(saTemp,saTempLoops,saVarLoops)    
    End If  
!
!  
!---------------------------------
! Fine Tune Points
!--------------------------------- 
    If(optLoops.gt.0)Then
! Save potential    
      fileName = "preFT.pot"  
      Call saveEamFile(fileName) 
! Fine tune points    
      Call M_synchProcesses()
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
        print *,"-------------------------------------------------"
        print *,"Fine tune nodes individually ",ProgramTime ()
        print *,"Loops: ",optLoops
        print *,"Total Vars: ",(optLoops*splineTotalNodes)
        print *,"-------------------------------------------------"
      End If
      Call fineTuneNodes(optLoops)
    End If    
!---------------------------------    
! Test optimum
!---------------------------------
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")   
      write(999,"(A70)") "----------------------------------------------------------------------"       
      write(999,"(A22)") "Optimisation finished."  
      write(999,"(A22)") "Final calculation.    "
      write(999,"(A70)") "----------------------------------------------------------------------"       
      Close(999)
    End If 
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt   
! Force embedding function to fit a 6th order polynomial
    Call forceEmbeddingFit()
! Make spline EAM from nodes
    Call setEamSpline()
! Apply ZBL hard core if required      
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      Call eamZblHardCore()
    End If  
    fileName = BlankString(fileName)
    fileName = "optimum.pot"
    Call saveEamFile(fileName)
! Reload force save option    
    saveForcesToFile = saveForcesToFileTemp
! Run evaluation      
    Call evaluate()  
! Output    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
      Print *,"Optimum EAM: ",totalRSS
    End If
! End Time
    Call cpu_time(timeEndOpt)  
! Store Time    
    Call storeTime(5,timeEndOpt-timeStartOpt)    
! Output    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
      Print *,"Optimisation time: ",(timeEndOpt-timeStartOpt)
    End If   
!---------------------------------    
! Analyse potential
!---------------------------------     
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,""
      Print *,""
      Print *,""
      Print *,"**********************************************************************"
      Print *,"Optimised Potential Test Results"
      Print *,"**********************************************************************"
    End If 
    Call runTestAnalysis()  
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"**********************************************************************"
    End If  
  End Subroutine runOptimiseProcess
!---------------------------------------------------------------------------------------------------  
  Subroutine fineTuneNodes(totalLoops)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID, i, j, loops, nodeToVary, totalLoops
    Integer(kind=StandardInteger) :: varsSinceImprov
    Integer(kind=StandardInteger), Dimension(1:splineTotalNodes) :: splineNodeList
    Integer(kind=StandardInteger), Dimension(1:splineTotalNodes,1:2) :: splineImprovement  ! successful var counts, successful + counts
    Real(kind=DoubleReal), Dimension(1:splineTotalNodes,1:2) :: splineImprovementDP        ! total rss reduction from node, overall variation
    Real(kind=DoubleReal) :: unvariedNode, variedNode
    Character(len=64) :: fileName
! Init variables    
    splineImprovement = 0  ! count of which node move made an improvement
    splineImprovementDP = 0.0D0  
    varsSinceImprov = 0
! Run first evaluation to get rss of optimum nodes
    splineNodesKey = splineNodesKeyOpt                    ! Set working nodes as optimum
    splineNodesData = splineNodesDataOpt
    Call forceEmbeddingFit()                              ! Force embedding function to fit 6th order poly
    Call setEamSpline()                                   ! Make spline EAM from nodes
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then  ! Apply ZBL hard core if required    
      Call eamZblHardCore()
    End If    
    Call evaluate()                                       ! Run evaluation     
    optimumRSS = totalRSS                                 ! Set optimum RSS 
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
      Print *,"Fine Tuning Start RSS (Opt): ",optimumRSS
    End If
! Start loop    
    j = 0
    Do loops=1,totalLoops
! Make a shuffled list of node points to vary
      Call makeSplineNodeList(splineNodeList)
! Rescale Embedding Function    
      splineNodesKey = splineNodesKeyOpt      ! Load nodes
      splineNodesData = splineNodesDataOpt   
      Call setEamSpline()                     ! Make spline EAM from nodes
      Call scaleEmbeddingOpt()                ! Rescale the embedding function for next loop
      splineNodesKeyOpt = splineNodesKey      ! Reload opt nodes
      splineNodesDataOpt = splineNodesData       
! Loop through nodes
      Do i=1,size(splineNodeList,1)
        nodeToVary = splineNodeList(i)      
        j = j + 1
! Set optimum spline points as working spline points
        splineNodesKey = splineNodesKeyOpt
        splineNodesData = splineNodesDataOpt
        If(varyFixedNodes.eq.1.or.&                       ! Update if node not fixed, or force vary fixed is on
        (varyFixedNodes.eq.0.and.&
        splineNodesData(nodeToVary,6).lt.0.5D0))Then    
! Vary point    
          unvariedNode = splineNodesData(nodeToVary,2)
          Call varyNodePoint(nodeToVary)
          variedNode = splineNodesData(nodeToVary,2)
! Save nodes to file
          If(loops.eq.1)Then
            Call outputSplineNodes("ft.nodes", 1)
          Else
            Call outputSplineNodes("ft.nodes", 0)
          End If            
! Force embedding function to fit a 6th order polynomial
          Call forceEmbeddingFit()
! Make spline EAM from nodes
          Call setEamSpline()
! Apply ZBL hard core if required      
          If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
            Call eamZblHardCore()
          End If  
! Run evaluation      
          Call evaluate()  
          If(totalRSS.lt.optimumRSS)Then 
! New spline nodes are an improvement
            varsSinceImprov = 0
! Record node specific data
            splineImprovement(nodeToVary,1) = splineImprovement(nodeToVary,1) + 1
            If(variedNode.gt.unvariedNode)Then
              splineImprovement(nodeToVary,2) = splineImprovement(nodeToVary,2) + 1
            End If
            splineImprovementDP(nodeToVary,1) = splineImprovementDP(nodeToVary,1)&
              + (optimumRSS-totalRSS)
            splineImprovementDP(nodeToVary,2) = splineImprovementDP(nodeToVary,2)&
              + (variedNode-unvariedNode)
! Set optimum nodes and rss
            splineNodesKeyOpt = splineNodesKey
            splineNodesDataOpt = splineNodesData  
            optimumRSS = totalRSS
            Do configID=1,configCount
              configVolumeOpt(configID) = configCalcEV(configID)
            End Do
! Save best potential so far            
            fileName = "potTempOptFT.pot"  
            Call saveEamFile(fileName)             
! Print to terminal
            If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
              Print *,"Variation RSS",j,i,totalRSS,(optimumRSS/startRSS),"*"
            End If 
          Else
! No improvement
            varsSinceImprov = varsSinceImprov + 1          
            If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
              Print *,"Variation RSS",j,i,totalRSS,(optimumRSS/startRSS)
            End If 
          End If   
        Else  
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
            Print *,"Variation RSS",j,i," (skipped, fixed)"
          End If 
        End If    
        If(varsSinceImprov.eq.100)Then  !Break out of loops
           Exit
        End If        
      End Do 
      If(varsSinceImprov.eq.100)Then  !Break out of loops
         Exit
      End If   
    End Do
! Print spline node variation data    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then  
      Print *,"-------------------------------------------------"
      Print *,"Node RSS improvement summary"
      Print *,"-------------------------------------------------"
      Do i=1,size(splineNodeList,1)
        print *,i,splineImprovement(i,1),splineImprovement(i,2),&
        splineImprovementDP(i,1),splineImprovementDP(i,2)
      End Do
    End If   
  End Subroutine fineTuneNodes
!---------------------------------------------------------------------------------------------------   
  Subroutine nodeResponseList()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i, runEvalSwitch
    Real(kind=DoubleReal) :: unperturbedRSS
    Do i=0,splineTotalNodes
! Set optimum spline points as working spline points
      splineNodesKey = splineNodesKeyOpt
      splineNodesData = splineNodesDataOpt     
      If(i.eq.0)Then
        ! no variation
      Else
        splineNodesData(i,2) = splineNodesData(i,2) + 0.001D0
      End If  
      runEvalSwitch = 0
      If(i.eq.0)Then
        runEvalSwitch = 1
      Else
        If(splineNodesData(i,6).lt.0.5D0)Then
          runEvalSwitch = 1
        End If
      End If
      If(runEvalSwitch.eq.1)Then
! Make spline EAM from nodes
        Call setEamSpline()
! Apply ZBL hard core if required      
        If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
          Call eamZblHardCore()
        End If 
! Run evaluation      
        Call evaluate()  
! Store
        If(i.eq.0)Then
          unperturbedRSS = totalRSS
        Else
          splineNodesResponse(i,1) = 0.001D0
          splineNodesResponse(i,2) = abs(unperturbedRSS-totalRSS)
! Print to terminal if set      
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
            Print *,i,abs(unperturbedRSS-totalRSS)
          End If
        End If  
      End If  
    End Do
  End Subroutine nodeResponseList
!---------------------------------------------------------------------------------------------------   
  Subroutine reduceNodeList()
    Implicit None   ! Force declaration of all variables
! Private variables    
    !Integer(kind=StandardInteger), Dimension(:) :: splineNodeList
    Integer(kind=StandardInteger) :: i, j, k, m, nodesStart
    Integer(kind=StandardInteger), Dimension(1:50,1:6) :: splineNodesKeyTemp
    Real(kind=DoubleReal), Dimension(1:10000,1:6) :: splineNodesDataTemp    
! Init arrays
    splineNodesKeyTemp = -1
    splineNodesDataTemp = 0.0D0
! Remove nodes with no response
    k = 0
    nodesStart = 1
    Do i=1,size(splineNodesKey)
      If(splineNodesKey(i,1).eq.-1)Then
        Exit
      End If
! Loop through nodes
      m = 0
      Do j=splineNodesKey(i,4),splineNodesKey(i,6)
        If(j.eq.splineNodesKey(i,4).or.&               ! Start node
        j.eq.splineNodesKey(i,6).or.&                  ! End node
        splineNodesResponse(j,2).gt.0.0D0.or.&         ! Node that creates a response when perturbed
        splineNodesDataTemp(j,6).gt.0.5D0)Then         ! Node that has been fixed
          k = k + 1
          m = m + 1
          splineNodesDataTemp(k,1) = splineNodesData(j,1)
          splineNodesDataTemp(k,2) = splineNodesData(j,2)
          splineNodesDataTemp(k,3) = splineNodesData(j,3)
          splineNodesDataTemp(k,4) = splineNodesData(j,4)
          splineNodesDataTemp(k,5) = 1.0D0 * m
          splineNodesDataTemp(k,6) = splineNodesData(j,6)
        End If
      End Do
! Store key data
      splineNodesKeyTemp(i,1) = splineNodesKey(i,1)
      splineNodesKeyTemp(i,2) = splineNodesKey(i,2)
      splineNodesKeyTemp(i,3) = splineNodesKey(i,3)
      splineNodesKeyTemp(i,4) = nodesStart
      splineNodesKeyTemp(i,5) = k - nodesStart + 1
      splineNodesKeyTemp(i,6) = k
      nodesStart = k + 1
    End Do
! Store reduced set of nodes    
    splineNodesKey = -1
    splineNodesData = 0.0D0
    splineNodesKey = splineNodesKeyTemp
    splineNodesData = splineNodesDataTemp
! Number of nodes    
    splineTotalNodes = k
! Output spline nodes
    Call outputSplineNodes("reduced.nodes")
  End Subroutine reduceNodeList      
!---------------------------------------------------------------------------------------------------   
  Subroutine jumbleNodes()
    Implicit None   ! Force declaration of all variables
! Private variables  
    Character(len=32) :: jumblePotName
    Integer(kind=StandardInteger) :: i, j, nodeToVary  
    Integer(kind=StandardInteger), Dimension(1:splineTotalNodes) :: splineNodeList
    Real(kind=DoubleReal) :: sigma, maxVariation
    Real(kind=DoubleReal) :: unvariedNode, variedNode
! Init variables     
    sigma = varyNodeOptions(2)    
    jumblePotName = BlankString(jumblePotName)
    jumblePotName = "jumble.pot"      
! Print to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Jumble all nodes:"
    End If
    Do i=1,10
! Make shuffled node list    
      Call makeSplineNodeList(splineNodeList)
! vary every node
      Do j=1,splineTotalNodes
! vary node value         
        nodeToVary = splineNodeList(j)
        If(varyFixedNodes.eq.1.or.&                       ! Update if node not fixed, or force vary fixed is on
        (varyFixedNodes.eq.0.and.&
        splineNodesData(nodeToVary,6).lt.0.5D0))Then     
          unvariedNode = splineNodesData(nodeToVary,2)    ! y value only, fix x value
          maxVariation = 0.1D0
          If(maxVariation.gt.(0.01D0*unvariedNode))Then
            maxVariation = (0.01D0*unvariedNode)  ! vary 0.1% of original value
          End If          
          variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)   
          splineNodesData(nodeToVary,2) = variedNode
! Update y'(x) and y''(x) using 4 point interp
          Call nodeGrad(nodeToVary)     
        End If          
      End Do
! Distribute array  
      Call M_distDouble2D(splineNodesData)  
! Force embedding function to fit a 6th order polynomial
      Call forceEmbeddingFit()
! Make spline EAM from nodes
      Call setEamSpline()
! Apply ZBL hard core if required      
      If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
        Call eamZblHardCore()
      End If        
! Run evaluation      
      Call evaluate() 
! Print to terminal if set      
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
        Print *,"Jumble",i,totalRSS
      End If      
    End Do
    Call saveEamFile(jumblePotName) 
  End Subroutine jumbleNodes  
!---------------------------------------------------------------------------------------------------   
  Subroutine fixZBLNodes()
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j, functionCounter
    Real(kind=DoubleReal) :: splineRightX
! Init variables    
    functionCounter = 0
    splineRightX = zblHardCore(2)
! Loop through EAM functions
    Do i=1,size(splineNodesKey,1)
      functionCounter = functionCounter + 1
      If(splineNodesKey(i,1).gt.0.and.splineNodesKey(i,3).eq.1)Then        ! If pair function        
        Do j=splineNodesKey(i,4),splineNodesKey(i,6)
          splineNodesData(j,6) = 1.0D0
          If(splineNodesData(j,1).gt.splineRightX)Then
            Exit
          End If
        End Do      
      End If  
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do
  End Subroutine fixZBLNodes  
!---------------------------------------------------------------------------------------------------   
  Subroutine makeSplineNodeList(splineNodeList)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger), Dimension(:) :: splineNodeList
    splineNodeList = 0
! Make shuffled array of nodes on root process only
    If(mpiProcessID.eq.0)Then
      splineNodeList = IntegerList(1,size(splineNodeList,1),size(splineNodeList,1))
    End If
! Distribute list to all processes
    Call M_distInt1D(splineNodeList)
  End Subroutine makeSplineNodeList
!---------------------------------------------------------------------------------------------------   
  Subroutine varyNodePoint(nodeToVary)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: nodeToVary, varyType 
    Real(kind=DoubleReal) :: unvariedNode, variedNode
    Real(kind=DoubleReal) :: sigma, maxVariation, maxVariationIn, maxLargeVariation
    Real(kind=DoubleReal) :: maxLargeVariationChance, randNumber
! Init variables    
    unvariedNode = splineNodesData(nodeToVary,2)    !y value only, fix x value
    variedNode = 0.0D0    
    varyType = Int(varyNodeOptions(1)) 
    sigma = varyNodeOptions(2) 
    maxVariation = varyNodeOptions(3) 
    maxLargeVariation = varyNodeOptions(4) 
    maxLargeVariationChance = varyNodeOptions(5)
! Vary on root process only, then distribute
    If(mpiProcessID.eq.0)Then  
! 1. Vary by fixed amount    
      If(varyType.eq.1)Then  
        variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)      
      End If   
! 2. vary by fixed amount, allow chance of larger fixed variation      
      If(varyType.eq.2)Then  
        Call RANDOM_NUMBER(randNumber)
        If(randNumber.le.maxLargeVariationChance)Then
          variedNode = RandomVaryPoint(unvariedNode, maxLargeVariation, sigma)   
        Else
          variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)  
        End If        
      End If
! 3. Vary by fixed amount    
      If(varyType.eq.3)Then  
        maxVariation = maxVariation * unvariedNode
        variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)      
      End If         
! 4. vary by fixed amount, allow chance of larger fixed variation      
      If(varyType.eq.4)Then  
        maxVariation = maxVariation * unvariedNode
        maxLargeVariation = maxLargeVariation * unvariedNode
        Call RANDOM_NUMBER(randNumber)
        If(randNumber.le.maxLargeVariationChance)Then
          variedNode = RandomVaryPoint(unvariedNode, maxLargeVariation, sigma)   
        Else
          variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)  
        End If        
      End If
! 5. Vary relative to nearby nodes   
      If(varyType.eq.5)Then  
        If(nodeToVary.eq.1)Then
          maxVariation = maxVariation * &
          abs(splineNodesData(2,2)-splineNodesData(1,2))
        Else If(nodeToVary.eq.splineTotalNodes)Then
          maxVariation = maxVariation * &
          abs(splineNodesData(splineTotalNodes,2)-&
          splineNodesData(splineTotalNodes-1,2))        
        Else
          maxVariation = maxVariation * &
          0.5D0*(abs(splineNodesData(splineTotalNodes,2)-&
          splineNodesData(splineTotalNodes-1,2))+&
          abs(splineNodesData(splineTotalNodes,2)-&
          splineNodesData(splineTotalNodes+1,2)))         
        End If
        variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)      
      End If  
! 6. Vary relative to nearby nodes   
      If(varyType.eq.6)Then  
        maxVariationIn = maxVariation
        If(nodeToVary.eq.1)Then
          maxVariation = maxVariation * &
          abs(splineNodesData(2,2)-splineNodesData(1,2))
        Else If(nodeToVary.eq.splineTotalNodes)Then
          maxVariation = maxVariation * &
          abs(splineNodesData(splineTotalNodes,2)-&
          splineNodesData(splineTotalNodes-1,2))        
        Else
          maxVariation = maxVariation * &
          0.5D0*(abs(splineNodesData(nodeToVary,2)-&
          splineNodesData(nodeToVary-1,2))+&
          abs(splineNodesData(nodeToVary,2)-&
          splineNodesData(nodeToVary+1,2)))        
        End If
        maxLargeVariation = maxLargeVariation*(maxVariation/maxVariationIn)
        Call RANDOM_NUMBER(randNumber)
        If(randNumber.le.maxLargeVariationChance)Then
          variedNode = RandomVaryPoint(unvariedNode, maxLargeVariation, sigma)   
        Else
          variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)  
        End If    
      End If    
      splineNodesData(nodeToVary,2) = variedNode
! Update y'(x) and y''(x) using 4 point interp
      Call nodeGrad(nodeToVary)
    End If  
! Distribute array  
    Call M_distDouble2D(splineNodesData)  
  End Subroutine varyNodePoint
!---------------------------------------------------------------------------------------------------   
  Subroutine nodeGrad(nodeToVary)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: nodeToVary, nodePos, interpOrder
    Integer(kind=StandardInteger) :: i, j, startPos, endPos, doLoop
    Real(kind=DoubleReal) :: x
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: interpPoints
    Real(kind=DoubleReal), Dimension(1:3) :: yArray    
! Update y'(x) and y''(x) using 4 point interp
    nodePos = Int(splineNodesData(nodeToVary,5))
    interpPoints = 0.0D0      
    startPos = nodeToVary-nodePos+1
! Loop through points
    doLoop = 1
    i = 0
    Do While(doLoop.eq.1)
      j = startPos + i
      i = i + 1
! Store data points
      interpPoints(i,1) = splineNodesData(j,1)
      interpPoints(i,2) = splineNodesData(j,2)
! break out of loop
      If(splineNodesData(j,5).gt.splineNodesData(j+1,5))Then
        endPos = j
        doLoop = 0
      End If
    End Do
    interpOrder = 4  ! Number of data points for interpolation
    If(i.lt.interpOrder)Then
      interpOrder = i
    End If
    x = splineNodesData(nodeToVary,1)
    yArray = PointInterp(interpPoints,x,4,2,1,(endPos-startPos+1))
! Update y'(x) and y''(x)
    splineNodesData(j,3) = yArray(2)
    splineNodesData(j,4) = yArray(3)
  End Subroutine nodeGrad
!---------------------------------------------------------------------------------------------------   
  Subroutine scaleEmbedding()
! Rescale the embedding function to match the density function  
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i, j, k, functionType, atomA
    Integer(kind=StandardInteger) :: eamStart, eamLength, eamEnd
    Real(kind=DoubleReal) :: rhoMin,rhoMax,embeMin,embeMax,xEmbe,xEmbeInc
    If(embeRescale.eq.1)Then
! Loop through EAM functions
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then        
          functionType = eamKey(i,3)
          eamStart = eamKey(i,4)
          eamLength = eamKey(i,5)
          eamEnd = eamKey(i,6)
          If(functionType.eq.2)Then   ! Dens
! Find range of density values
            atomA = eamKey(i,1)
            rhoMin = 2.1D20
            rhoMax = -2.1D20
            Do j=eamStart,eamEnd
              If(eamData(j,2).lt.rhoMin)Then
                rhoMin = eamData(j,2)
              End If
              If(eamData(j,2).gt.rhoMax)Then
                rhoMax = eamData(j,2)
              End If
            End Do
! Loop through functions to find matching embedding
            Do j=1,size(eamKey,1)
              If(eamKey(j,3).eq.3.and.eamKey(j,1).eq.atomA)Then  ! If embe function and the right atom species
! Assume density starts at zero
                rhoMin = 0.0D0
                xEmbe = 0.0D0
                xEmbeInc = rhoMax/(1.0D0*(splineNodesKey(j,5)-1.0D0))
! Loop through nodes
                Do k = splineNodesKey(j,4),splineNodesKey(j,6)
                  splineNodesData(k,1) = xEmbe
                  xEmbe = xEmbe + xEmbeInc
                End Do
                embeMin = splineNodesData(splineNodesKey(j,4),1)
                embeMax = splineNodesData(splineNodesKey(j,6),1)
                Call outputEmbeRescale(i,j,rhoMin,rhoMax,embeMin,embeMax)
              End If
            End Do             
          End If
        End If
      End Do
    End If
  End Subroutine scaleEmbedding 
!--------------------------------------------------------------------------------------------------- 
  Subroutine scaleEmbeddingOpt()
! Rescale the embedding function to match the density function  
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i, j, k, functionType, atomA
    Integer(kind=StandardInteger) :: eamStart, eamLength, eamEnd
    Real(kind=DoubleReal) :: rhoMin,rhoMax,embeMin,embeMax,xEmbe,xEmbeInc
    If(embeRescale.eq.1)Then
! Loop through EAM functions
      Do i=1,size(eamKeyOpt,1)
        If(eamKeyOpt(i,1).gt.0)Then        
          functionType = eamKeyOpt(i,3)
          eamStart = eamKeyOpt(i,4)
          eamLength = eamKeyOpt(i,5)
          eamEnd = eamKeyOpt(i,6)
          If(functionType.eq.2)Then   ! Dens
! Find range of density values
            atomA = eamKeyOpt(i,1)
            rhoMin = 2.1D20
            rhoMax = -2.1D20
            Do j=eamStart,eamEnd
              If(eamDataOpt(j,2).lt.rhoMin)Then
                rhoMin = eamDataOpt(j,2)
              End If
              If(eamDataOpt(j,2).gt.rhoMax)Then
                rhoMax = eamDataOpt(j,2)
              End If
            End Do
! Loop through functions to find matching embedding
            Do j=1,size(eamKeyOpt,1)
              If(eamKeyOpt(j,3).eq.3.and.eamKeyOpt(j,1).eq.atomA)Then  ! If embe function and the right atom species
! Assume density starts at zero
                rhoMin = 0.0D0
                xEmbe = 0.0D0
                xEmbeInc = rhoMax/(1.0D0*(splineNodesKeyOpt(j,5)-1.0D0))
! Loop through nodes
                Do k = splineNodesKeyOpt(j,4),splineNodesKeyOpt(j,6)
                  splineNodesDataOpt(k,1) = xEmbe
                  xEmbe = xEmbe + xEmbeInc
                End Do
                embeMin = splineNodesDataOpt(splineNodesKeyOpt(j,4),1)
                embeMax = splineNodesDataOpt(splineNodesKeyOpt(j,6),1)
                Call outputEmbeRescale(i,j,rhoMin,rhoMax,embeMin,embeMax)
              End If
            End Do             
          End If
        End If
      End Do
    End If
  End Subroutine scaleEmbeddingOpt
!---------------------------------------------------------------------------------------------------  
  Subroutine forceEmbeddingFit()
! Force EAM Embedding Function to fit an analytic form
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, n, functionCounter     
    Real(kind=DoubleReal), Dimension(1:1000,1:4) :: dataPoints
! Loop through EAM functions
    functionCounter = 0
    If(forceEmbeFitOpt.eq.1)Then
      dataPoints = 0.0D0    
      Do i=1,size(splineNodesKey,1)
        If(splineNodesKey(i,1).gt.0)Then        
          functionCounter = functionCounter + 1  
        End If
        If(splineNodesKey(i,3).eq.3.or.splineNodesKey(i,3).eq.6&
          .or.splineNodesKey(i,3).eq.7)Then
          n = 0
          Do j=splineNodesKey(i,4),splineNodesKey(i,6)
            n = n + 1
            dataPoints(n,1) = splineNodesData(j,1)
            dataPoints(n,2) = splineNodesData(j,2)
          End Do
          Call forceEmbeddingFitProcess(dataPoints, n)   
          n = 0     
          Do j=splineNodesKey(i,4),splineNodesKey(i,6)
            n = n + 1
            splineNodesData(j,1) = dataPoints(n,1) 
            splineNodesData(j,2) = dataPoints(n,2) 
            splineNodesData(j,3) = dataPoints(n,3) 
            splineNodesData(j,4) = dataPoints(n,4) 
          End Do
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do 
    End If
  End Subroutine forceEmbeddingFit   
!---------------------------------------------------------------------------------------------------  
  Subroutine forceEmbeddingFitProcess(dataPointsIn, dataPointsCounter)
! Force EAM Embedding Function to fit an analytic form
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, dataPointsCounter      
    Real(kind=DoubleReal), Dimension(1:1000,1:4) :: dataPointsIn
    Real(kind=DoubleReal), Dimension(1:dataPointsCounter,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:dataPointsCounter,1:2) :: pointsInterp
    Real(kind=DoubleReal), Dimension(1:7) :: coefficients
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Transfer data
    Do i=1,dataPointsCounter
      dataPoints(i,1) = dataPointsIn(i,1)
      dataPoints(i,2) = dataPointsIn(i,2)
    End Do
! Fit to polynomial
    coefficients = PolyFit(dataPoints,6)
! Replace points using fit
    Do i=1,dataPointsCounter    
      dataPointsIn(i,2) = CalcPolynomial(coefficients, dataPointsIn(i,1))
      pointsInterp(i,1) = dataPointsIn(i,1)
      pointsInterp(i,2) = dataPointsIn(i,2)
    End Do
! Interpolate to give f'(x) and f''(x)   
    Do i=1,dataPointsCounter   
      yArray = PointInterp(pointsInterp,pointsInterp(i,1),&
      4,2,1,dataPointsCounter)
      dataPointsIn(i,3) = yArray(2)
      dataPointsIn(i,4) = yArray(3)
      !If(mpiProcessID.eq.0)Then
      !print *,dataPointsIn(i,1),dataPointsIn(i,2),dataPointsIn(i,3),dataPointsIn(i,4)
      !End If
    End Do
    
  End Subroutine forceEmbeddingFitProcess  
!---------------------------------------------------------------------------------------------------   
  Subroutine simulatedAnnealing(sTemp, tLoops, varLoops)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Real(kind=DoubleReal) :: sTemp
    Integer(kind=StandardInteger) :: i, j, accepted, tLoops, varLoops 
    Real(kind=DoubleReal) :: tLoop, badAcceptance, randNumber    
    Character(len=64) :: fileName
! Init variables
    badAcceptance = 0.0D0    
! Run first evaluation to get rss of optimum nodes
    splineNodesKey = splineNodesKeyOpt                    ! Set working nodes as optimum
    splineNodesData = splineNodesDataOpt
    Call forceEmbeddingFit()                              ! Force embedding function to fit 6th order poly
    Call setEamSpline()                                   ! Make spline EAM from nodes
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then  ! Apply ZBL hard core if required    
      Call eamZblHardCore()
    End If    
    Call evaluate()                                       ! Run evaluation     
    optimumRSS = totalRSS                                 ! Set optimum RSS  
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
      Print *,"SA Start RSS (Opt): ",optimumRSS
    End If
! Temperature loop
    Do i=1,tLoops                                         ! Loop through decreasing temperatures
      tLoop = sTemp / (2.0D0**(i-1))                      ! Loop temperature
      Do j=1,varLoops                                     ! Loop through variations at this temperature        
        accepted = 0                                  
        splineNodesKey = splineNodesKeyOpt                ! Set working nodes as optimum
        splineNodesData = splineNodesDataOpt
        Call saVary()                                     ! Vary nodes
        If(i.eq.1.and.j.eq.1)Then
          Call outputSplineNodes("sa.nodes", 1)
        Else
          Call outputSplineNodes("sa.nodes", 0)
        End If        
        Call forceEmbeddingFit()                          ! Force embedding function to fit 6th order poly
        Call setEamSpline()                               ! Make spline EAM from nodes
        If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then  ! Apply ZBL hard core if required    
          Call eamZblHardCore()
        End If    
        Call evaluate()                                   ! Run evaluate        
        If(totalRSS.lt.optimumRSS)Then                    ! If better, accept new points
          splineNodesKeyOpt = splineNodesKey
          splineNodesDataOpt = splineNodesData
          optimumRSS = totalRSS
          accepted = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
            Print *,"SA ",tLoop,j,totalRSS,accepted,(optimumRSS/startRSS),"*"
          End If 
          ! Save best potential so far            
          fileName = "potTempOptSA.pot"  
          Call saveEamFile(fileName)     
        Else
          badAcceptance = &
          exp(-1.0D0*abs((totalRSS-optimumRSS)/optimumRSS)/tLoop)
          Call RANDOM_NUMBER(randNumber)
          If(randNumber.le.badAcceptance)Then             ! Accept bad by chance 
            splineNodesKeyOpt = splineNodesKey
            splineNodesDataOpt = splineNodesData
            optimumRSS = totalRSS
            accepted = 1
          End If
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
            Print *,"SA ",tLoop,j,totalRSS,accepted,(optimumRSS/startRSS)
          End If 
        End If
      End Do                                              ! End loop through variations
    End Do                                                ! End Temperature loop
  End Subroutine simulatedAnnealing  
!---------------------------------------------------------------------------------------------------   
  Subroutine saVary()
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: j, nodeToVary  
    Integer(kind=StandardInteger), Dimension(1:splineTotalNodes) :: splineNodeList
    Real(kind=DoubleReal) :: sigma, maxVariation
    Real(kind=DoubleReal) :: unvariedNode, variedNode
! Init variables     
    sigma = varyNodeOptions(2)    
! Set optimum spline points as working spline points
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt   
! Make shuffled node list    
    Call makeSplineNodeList(splineNodeList)
! vary every node
    Do j=1,splineTotalNodes
! vary node value      
      nodeToVary = splineNodeList(j)
      If(varyFixedNodes.eq.1.or.&                       ! Update if node not fixed, or force vary fixed is on
      (varyFixedNodes.eq.0.and.&
      splineNodesData(nodeToVary,6).lt.0.5D0))Then  
        unvariedNode = splineNodesData(nodeToVary,2)    !y value only, fix x value
        maxVariation = saMaxVariation * unvariedNode  ! vary 0.1% of original value
        variedNode = RandomVaryPoint(unvariedNode, maxVariation, sigma)   
        splineNodesData(nodeToVary,2) = variedNode
! Update y'(x) and y''(x) using 4 point interp
        Call nodeGrad(nodeToVary)      
      End If
    End Do
! Distribute array  
    Call M_distDouble2D(splineNodesData)  
! Make spline EAM from nodes
    Call setEamSpline()
! Apply ZBL hard core if required      
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      Call eamZblHardCore()
    End If        
  End Subroutine saVary    
  
  
  
  
  
  
  
  
  
  
  
  
!--------------------------------------------------------------------------------------------------- 
!  Test subroutines
!---------------------------------------------------------------------------------------------------   
  Subroutine scaleEmbeddingTest()
! Rescale the embedding function to match the density function  
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i, j, k, functionType, atomA
    Integer(kind=StandardInteger) :: eamStart, eamLength, eamEnd
    Real(kind=DoubleReal) :: rhoMin,rhoMax,embeMin,embeMax,xEmbe,xEmbeInc
    If(embeRescale.eq.1)Then
! Loop through EAM functions
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then        
          functionType = eamKey(i,3)
          eamStart = eamKey(i,4)
          eamLength = eamKey(i,5)
          eamEnd = eamKey(i,6)
          If(functionType.eq.2)Then   ! Dens
! Find range of density values
            atomA = eamKey(i,1)
            rhoMin = 2.1D20
            rhoMax = -2.1D20
            Do j=eamStart,eamEnd
              If(eamData(j,2).lt.rhoMin)Then
                rhoMin = eamData(j,2)
              End If
              If(eamData(j,2).gt.rhoMax)Then
                rhoMax = eamData(j,2)
              End If
            End Do
            !print *,"Density range ",rhoMin,rhoMax
            Do j=1,size(eamKey,1)
              If(eamKey(j,3).eq.3.and.eamKey(j,1).eq.atomA)Then  ! If embe function and the right atom species
! Assume density starts at zero
                rhoMin = 0.0D0
                xEmbe = 0.0D0
                xEmbeInc = rhoMax/(1.0D0*(eamKey(j,5)-1.0D0))
                !print *,"embe inc",xEmbeInc
                Do k = eamKey(j,4),eamKey(j,6)
                  eamData(k,1) = xEmbe
                  xEmbe = xEmbe + xEmbeInc
                End Do
                embeMin = eamData(eamKey(j,4),1)
                embeMax = eamData(eamKey(j,6),1)
                Call outputEmbeRescale(i,j,rhoMin,rhoMax,embeMin,embeMax)
              End If
            End Do             
          End If
        End If
      End Do
    End If
  End Subroutine scaleEmbeddingTest
  
  
  
End Module optimise