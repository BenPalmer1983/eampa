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
  Use eamGen
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
    Implicit None   ! Force declaration of all variables
! Private variables
    Type(saConfig) :: saConfigLive
    Integer(kind=StandardInteger) :: i
! output
    If(TerminalPrint())Then
      print *,""
      print *,""
      print *,"                           Optimise"
      print *,"----------------------------------------------------------------------"
    End If 
! Load EAM stored in eamKeyInput/eamDataInput
    Call loadInputEAM()
! Run initial efs and bp calculations    
    Call runCalcs()
    If(TerminalPrint())Then
      print *,"Input EAM RSS value:"
      print *,totalRSS
    End If
    Call outputBpT()
    Call saveEamFile("opt_001_input.pot")
    
! Spline functions
    Call setEamNodes()         ! Choose the spline nodes from the full eam data for the functions
    Call completeEAMSpline(splineNodesKey, splineNodesData) ! Fill in y'(x) and y''(x)
    splineNodesKeyOpt = splineNodesKey
    splineNodesDataOpt = splineNodesData
    Call setEamSpline()        ! Spline the nodes to make new eam functions
    If(optForceZBL)Then
      Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
    End If    
    Call runCalcs()            ! evalEAM() + evalBulkProperties()
    startRSS = totalRSS
    optimumRSS = totalRSS
    If(TerminalPrint())Then
      print *,"Input EAM RSS value (splined):"
      print *,totalRSS
    End If   
    Call outputBpT()
    Call saveEamFile("opt_002_starting.pot")
! Output to terminal - start of SA
    If(TerminalPrint())Then
      print *,"Start Simulated Annealing"
      print *,"Settings:"
      print *,"Temperature:              ",saConfigIn%temp
      print *,"Temperature loops:        ",saConfigIn%tempLoops
      print *,"Variation loops:          ",saConfigIn%varLoops
      print *,"Maximum variation factor: ",saConfigIn%maxVar    
      print *,"Refinement Loops:         ",saConfigIn%refinementLoops    
    End If    
! Run SA subroutine    
    If(saConfigIn%refinementLoops.gt.0)Then
      Do i=1,saConfigIn%refinementLoops
! copy from input SA settings
        saConfigLive = saConfigIn   
! Update temperature        
        saConfigLive%temp = saConfigLive%temp / (2**((i-1)*saConfigIn%tempLoops))
! SA optimisation
        Call saOpt(saConfigLive)
! Reduce        
        saConfigLive%maxVar = saConfigLive%maxVar/2.0D0
      End Do
    End If

    
    Call saveEamFile("opt_003_postSA.pot")
    
    
    Call LM_Opt()
    
    
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt
    Call setEamSpline()        ! Spline the nodes to make new eam functions
    If(optForceZBL)Then
      Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
    End If    
    Call runCalcs()
    If(TerminalPrint())Then
      print *,"Optimised:"
      print *,totalRSS,optimumRSS
    End If   
    Call outputBpT()
    Call saveEamFile("opt_004_optimised.pot")
    

    
  End Subroutine optiEAM
  
   
  
! --------------------------------------------------------------------------------------------------- 
  Subroutine runCalcs()
! Calculate stress/energy/force of configuration/s   
    Implicit None   ! Force declaration of all variables
! Private variables 
    quietOverride = .true.     ! Turn off verbose
    Call evalEAM()
! Bulk Properties 
    Call evalBulkProperties()
    quietOverride = .false. 
! Save data to file    
    Call outputOptLine()
  End Subroutine runCalcs 
  
 
! ---------------------------------------------------------------------------------------------------  
!
! Simulated Annealing Functions 
!  
! --------------------------------------------------------------------------------------------------- 
  Subroutine saOpt(saConfigLive)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i
    Type(saConfig) :: saConfigLive
! Loop through and decrease temperature
    Do i=1,saConfigLive%tempLoops
      Call saOpt_VarLoop(saConfigLive, i)  
    End Do
  End Subroutine saOpt
! --------------------------------------------------------------------------------------------------- 
  Subroutine saOpt_VarLoop(saConfigLive, loop)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Real(kind=DoubleReal) :: temperature
    Integer(kind=StandardInteger) :: i, loop
    Real(kind=DoubleReal) :: aProb, randDouble
    Type(saConfig) :: saConfigLive
    Logical :: accept, bad
    If(TerminalPrint())Then
      print *,""
      print *,"Start Simulated Annealing Loop ",loop
      print *,""
    End If    
! Temperature    
    temperature = (saConfigLive%temp)/(2.0D0**(loop-1))
! Loop and vary nodes
    Do i=1,saConfigLive%varLoops
! Load optimum nodes    
      splineNodesKey = splineNodesKeyOpt
      splineNodesData = splineNodesDataOpt
! Vary the spline nodes on the root process
      If(mpiProcessID.eq.0)Then
        Call varyEAMNodes(saConfigLive%maxVar, splineNodesKey, splineNodesData)
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
! Acceptance - SA - Root process only
      bad = .false.
      accept = .false.
      If(mpiProcessID.eq.0)Then
        If(totalRSS.lt.optimumRSS)Then  ! improvement
          accept = .true.
        Else  
          aProb = exp((-1.0D0*(totalRSS-optimumRSS))/temperature)
          Call RANDOM_NUMBER(randDouble)
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
      End If
! Print out      
      If(TerminalPrint())Then
        If(accept)Then
          If(bad)Then
            print *,i,": [",optimumRSS,"]  ",totalRSS,"*(Bad)"
          Else
            print *,i,": [",optimumRSS,"]  ",totalRSS,"*(Good)"
          End If          
        Else 
          print *,i,": [",optimumRSS,"]  ",totalRSS
        End If        
      End If       
    End Do  
! Test optimum
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt
    Call setEamSpline()        ! Spline the nodes to make new eam functions
    If(optForceZBL)Then
      Call eamPairZbl()        ! Force overwriting core with ZBL - exp better fit than spline at small r
    End If    
    Call runCalcs()
    If(TerminalPrint())Then
      print *,"Optimised loop ",loop
      print *,totalRSS,optimumRSS
    End If 
! Output Bulk Properties       
    Call outputBpT()    
  End Subroutine saOpt_VarLoop
  
  
! ---------------------------------------------------------------------------------------------------  
!
! Levenberg-Marquardt Algorithm
!  
! --------------------------------------------------------------------------------------------------- 
  Subroutine LM_Opt()
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j, totalNodes
! total nodes
    totalNodes = 0
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then    
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
          totalNodes = totalNodes + 1
        End Do
      Else
        Exit
      End If
    End Do 

    If(mpiProcessID.eq.0)Then
      print *,"Total Nodes: ",totalNodes
    End If

    Call LM_OptProcess(totalNodes)
    
    
    
    
    
    
    
  End Subroutine LM_Opt
  
  
  
  
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
    Real(kind=DoubleReal) :: refRSS, rss, lastRSS
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
        splineNodesResponse = FillSplineResponse(splineNodesResponse,splineNodesKey(i,4),splineNodesKey(i,6))        
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

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

End Module opti
