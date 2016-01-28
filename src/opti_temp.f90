
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
    
    Call setEamNodes()                         ! readEAM.f90   set spline nodes from eamKey/eamData   
    optEmbeddingFit = 0                        ! Force embedding functional to fit form - on    
    optDensityFit = 0                          ! Spline poly for density - off   
    Call setEamSpline(.true.)                  ! readEAM.f90   make EAM potential by splining nodes
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