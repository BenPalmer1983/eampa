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
    !Type(saConfig) :: saConfigLive
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
  
  
    
    
    
! Load EAM stored in eamKeyInput/eamDataInput
    Call loadInputEAM()
! Run initial efs and bp calculations    
    Call runCalcs(0)
    If(TerminalPrint())Then
      print *,"Input RSS ",totalRSS
    End If  
    Call saveEamFile("opt_001_input.pot")
    
! Set EAM nodes

    Call setEamNodesOpti()
    Call setEamSplineOpti()
    
    Call runCalcs(1) 
    
    If(TerminalPrint())Then
      print *,"Spline converted input RSS ",totalRSS
    End If  
    Call saveEamFile("opt_002_input-splined.pot")
    
    
    
    LM_nodes = 0
    Do i=1,eamFunctionCount
      If(splineNodesKey(i,3).eq.1)Then ! pair
        LM_nodes = LM_nodes + splineNodesKey(i,5) - 2      
      End If
      If(splineNodesKey(i,3).eq.2.or.splineNodesKey(i,3).eq.4.or.splineNodesKey(i,3).eq.5)Then ! dens/dden/sden
        LM_nodes = LM_nodes + splineNodesKey(i,5) - 1      
      End If
      If(splineNodesKey(i,3).eq.3.or.splineNodesKey(i,3).eq.6.or.splineNodesKey(i,3).eq.7)Then ! embe/demb/semb
        LM_nodes = LM_nodes + splineNodesKey(i,5)    
      End If
    End Do
    
    print *,splineTotalNodes, LM_nodes
    
    
    Call LM_Opt(crCount,LM_nodes)
    
    
    
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
  Subroutine runCalcs(splineNodesIn)
! Calculate stress/energy/force of configuration/s   
    Implicit None   ! Force declaration of all variables
! Private variables 
    Integer(kind=StandardInteger), Optional :: splineNodesIn
    Integer(kind=StandardInteger) :: splineNodes
! Optional Arguments
    If(Present(splineNodesIn))Then
      splineNodes = splineNodesIn
    End If    
    If(splineNodes.eq.1)Then
      Call setEamSplineOpti() 
    End If        
! Turn off verbose    
    quietOverride = .true.     
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
  Subroutine LM_Opt(sampleCount,totalNodes)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: sampleCount,totalNodes
    Integer(kind=StandardInteger) :: f, i, n, k, m
    Real(kind=DoubleReal) :: varyAmount
    Real(kind=DoubleReal) :: refRSS, lastRSS
    Logical :: unfixedP, calcJR
! Matrix
    Real(kind=DoubleReal), Dimension(1:sampleCount,1:totalNodes) :: J, J_Last
    Real(kind=DoubleReal), Dimension(1:sampleCount) :: refVals, derivRefVals, R, R_Last
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:sampleCount) :: JT     ! Transpose Jacobian
    Real(kind=DoubleReal), Dimension(1:totalNodes,1:totalNodes) :: JTJ, JTJ_Diag  
    Real(kind=DoubleReal), Dimension(1:totalNodes) :: JTR, P, P_B, xP, xP_Last
    Integer(kind=StandardInteger), Dimension(1:totalNodes,1:2) :: xP_map
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
! -----------------
! Start LMA Loop
!
    calcJR = .true.
    Do n=1,4
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
          varyAmount = abs(0.0001D0*xP(i))
          If(varyAmount.eq.0.0D0)Then
            varyAmount = 0.000001D0
          End If
          splineNodesData(k,2) = xP(i)+varyAmount ! Perturb node
          Call CompleteNodeData(splineNodesData, splineNodesKey(f,4), splineNodesKey(f,6))
          Call runCalcs(1)
          splineNodesData(k,2) = xP(i)            ! reset node
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
      JTJ = InvertMatrix(JTJ) ! store inverse (recycle JTJ var)      
      JTR = matmul(JT,R)
      JTR = -1.0D0*JTR ! Recycle JTR var
      P = matmul(JTJ,JTR) 
      
      JTJ = matmul(JT,J)
      JTJ_Diag = lambda*DiagMatrix(JTJ) ! Dampening Matrix
      JTJ = MatAdd(JTJ,JTJ_Diag) ! Recycle JTJ
      P_B = SolveLinearSet(JTJ,JTR)
      
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
        xP = xP + P(i)
        k = xP_map(i,1)
        splineNodesData(k,2) = splineNodesData(k,2)+P(i)
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
    
    
  End Subroutine LM_Opt
  
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
