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
  Use calcEval  

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
    Integer(kind=StandardInteger) :: i, j, loops
    Character(len=255) :: eamTempPotPath
    Character(len=64) :: fileName
    Character(len=32) :: fileNumber
    Integer(kind=StandardInteger), Dimension(1:splineTotalNodes) :: splineNodeList
    Real(kind=DoubleReal) :: timeStartOpt, timeEndOpt
! Start Time
    Call cpu_time(timeStartOpt)
! Print out to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Start Optimise EAM"
    End If
! Make directory for temp potentials to debug (if required)
    eamTempPotPath = Trim(outputDirectory)//"/pots"
    Call makeDir(eamTempPotPath)
! Load original input potential    
    eamKey = eamKeyInput
    eamData = eamDataInput
! Get starting nodes    
    Call setEamNodes()
! Save starting nodes as optimum    
    splineNodesKeyOpt = splineNodesKey
    splineNodesDataOpt = splineNodesData   
! Force ZBL if required
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      Call eamZblHardCore()
    End If    
!-------------------------------------  
! Evaluate Input EAM    
!-------------------------------------
    fileName = "potInput.pot"
    Call saveEamFile(fileName)
    Call evaluate()    
    optimumRSS = totalRSS
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"Input RSS: ",totalRSS,"*"
    End If  
!-------------------------------------      
! Spline EAM - Evaluate Splined EAM
!-------------------------------------
    Call setEamSpline()
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      Call eamZblHardCore()
    End If  
    fileName = "potSpline.pot"  
    Call saveEamFile(fileName) 
    Call evaluate()         
    If(totalRSS.lt.optimumRSS)Then
      splineNodesKeyOpt = splineNodesKey
      splineNodesDataOpt = splineNodesData   
      optimumRSS = totalRSS
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
        Print *,"Spline RSS: ",totalRSS,"*"
      End If 
    Else
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
        Print *,"Spline RSS: ",totalRSS
      End If 
    End If  
!-------------------------------------      
! Loop and vary position of spline nodes
!-------------------------------------   
    Call makeSplineNodeList(splineNodeList)
    j = 0
    Do loops=1,10
      Do i=1,size(splineNodeList,1)
        j = j + 1
! Set optimum spline points as working spline points
        splineNodesKey = splineNodesKeyOpt
        splineNodesData = splineNodesDataOpt
! Vary point    
        Call varyNodePoint(splineNodeList(i))
! Make spline EAM from nodes
        Call setEamSpline()
! Apply ZBL hard core if required      
        If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
          Call eamZblHardCore()
        End If  
! Save potential
        fileNumber = BlankString(fileNumber)
        fileNumber = intToString(i)
        fileName = "pots/pot"//Trim(Adjustl(fileNumber))
        fileName = Trim(Adjustl(fileName))//".pot"  
        Call saveEamFile(fileName) 
! Run evaluation      
        Call evaluate()  
        If(totalRSS.lt.optimumRSS)Then
          splineNodesKeyOpt = splineNodesKey
          splineNodesDataOpt = splineNodesData   
          optimumRSS = totalRSS
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
            Print *,"Variation RSS",j,i,totalRSS,"*"
          End If 
        Else
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then    
            Print *,"Variation RSS",j,i,totalRSS
          End If 
        End If   
      End Do
    End Do
! Test optimum
    splineNodesKey = splineNodesKeyOpt
    splineNodesData = splineNodesDataOpt   
! Make spline EAM from nodes
    Call setEamSpline()
! Apply ZBL hard core if required      
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      Call eamZblHardCore()
    End If  
    fileName = BlankString(fileName)
    fileName = "optimum.pot"
    Call saveEamFile(fileName)
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
    
  End Subroutine runOptimise
  
  
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
! Store value      
      splineNodesData(nodeToVary,2) = variedNode
    End If  
! Distribute array  
    Call M_distDouble2D(splineNodesData)  
  End Subroutine varyNodePoint
  
  
  

End Module optimise