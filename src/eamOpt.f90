Module eamOpt

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

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
  Use initialise
  Use loadData
  Use globals
  Use output
  Use eamGen
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines






  Contains
  
  Subroutine setEamNodesOpti()
! Set node positions using the currently loaded eam functions in eamKey/eamData
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
    Integer(kind=StandardInteger) :: eamStart, eamLength, eamEnd
    Integer(kind=StandardInteger) :: nodes, nodeKey, functionType
    Real(kind=DoubleReal) :: x, xStart, xEnd
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Integer(kind=StandardInteger) :: zA, zB
! Init variables
    splineNodesKey = -1         ! reset key array
    splineNodesData = 0.0D0     ! reset node array
    functionCounter = 0
    nodeKey = 0
! Loop through EAM functions
    Do i=1,size(eamKey,1)
      If(eamKey(i,1).gt.0)Then
        functionCounter = functionCounter + 1
        nodes = splineNodeCount(eamKey(i,3))
        FunctionType = eamKey(i,3)
        eamStart = eamKey(i,4)
        eamLength = eamKey(i,5)
        eamEnd = eamKey(i,6)
! xStart xEnd
        xStart = eamData(eamStart,1)       
        xEnd = eamData(eamEnd,1) 
        If(eamKey(i,3).eq.1)Then ! Function type = pair
          xStart = zblHardCore(1)
          zA = elementsCharge(eamKey(i,1))
          zB = elementsCharge(eamKey(i,2))          
        End If
! Loop through nodes for each function
        Do j=1,nodes
          nodeKey = nodeKey + 1
!If first node, save node key data
          If(j.eq.1)Then  
            splineNodesKey(i,1) = eamKey(i,1)           ! Type A
            splineNodesKey(i,2) = eamKey(i,2)           ! Type B
            splineNodesKey(i,3) = eamKey(i,3)           ! Function type
            splineNodesKey(i,4) = nodeKey               ! Start
            splineNodesKey(i,5) = nodes                 ! Length
            splineNodesKey(i,6) = nodeKey + nodes - 1   ! End
          End If
! x point
          x = xStart+1.0D0*(j-1)*((xEnd-xStart)/(nodes-1))
! last node and Pair/Dens type 0,0,0
          If(j.eq.nodes.and.&
          (eamKey(i,3).eq.1.or.eamKey(i,3).eq.2.or.eamKey(i,3).eq.4.or.eamKey(i,3).eq.5)&
          )Then
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = 0.0D0
            splineNodesData(nodeKey,3) = 0.0D0
            splineNodesData(nodeKey,4) = 0.0D0
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0
          ElseIf(j.eq.1.and.eamKey(i,3).eq.1)Then  ! Pair - get ZBL y(x), y'(x) and y''(x)
            yArray = ZblFull (x, zA, zB)            
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = yArray(1)
            splineNodesData(nodeKey,3) = yArray(2)
            splineNodesData(nodeKey,4) = yArray(3)
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0        
          Else
            yArray = PointInterp(eamData,x,eamInterpPoints,2,eamStart,eamLength)
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = yArray(1)
            splineNodesData(nodeKey,3) = yArray(2)
            splineNodesData(nodeKey,4) = yArray(3)
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0
          End If
        End Do
      End If
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do
! Store total number of nodes
    splineTotalNodes = nodeKey
  End Subroutine setEamNodesOpti
  
  
  
End Module eamOpt  