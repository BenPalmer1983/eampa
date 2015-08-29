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
  Public :: setEamNodesOpti
  Public :: setEamSplineOpti






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
  
  
  Subroutine setEamSplineOpti()
! Spline between nodes and store functions in eamKey/eamData
! Force ZBL Core for pair potentials
! exp(a+bx) style spline from  zbl to first node
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
    Integer(kind=StandardInteger) :: nodes, nodeStart, nodeLength, nodeEnd, eamPoint
    Real(kind=DoubleReal), Dimension(1:1001,1:4) :: splineDataPoints
    Real(kind=DoubleReal) :: x, changeX
    Integer(kind=StandardInteger) :: zA, zB, nodesZ, nodesA
    Integer(kind=StandardInteger) :: pointsZBL, pointsSpline, pointsExpSpline, pointsPolySpline
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Integer(kind=StandardInteger), Dimension(1:1000) :: splineType
! Init variables
    eamKey = 0                 ! clear eam key
    eamData = 0.0D0            ! clear eam data
    splineDataPoints = 0.0D0   ! init temp spline points array
! Loop through EAM functions
    functionCounter = 0
    Do i=1,size(splineNodesKey,1)    
      If(splineNodesKey(i,1).gt.0)Then ! Check that eam function is stored    
! Count function      
        functionCounter = functionCounter + 1
! Update Key Data        
        eamKey(functionCounter,1) = splineNodesKey(i,1) 
        eamKey(functionCounter,2) = splineNodesKey(i,2) 
        eamKey(functionCounter,3) = splineNodesKey(i,3) 
        eamKey(functionCounter,4) = 1+(1001*(functionCounter-1))
        eamKey(functionCounter,5) = 1001
        eamKey(functionCounter,6) = eamKey(functionCounter,4)+eamKey(functionCounter,5)-1
! Nodes used in spline
        nodes = splineNodeCount(splineNodesKey(i,3))
! Pair Potentials
        If(splineNodesKey(i,3).eq.1)Then       
! ZBL
          zA = elementsCharge(splineNodesKey(i,1))
          zB = elementsCharge(splineNodesKey(i,2))   
! Set node start/end points
          nodeStart = splineNodesKey(i,4)
          nodeLength = splineNodesKey(i,5)
          nodeEnd = splineNodesKey(i,6)
! what part of function will be spline - function runs from 0.0 to splineNodesData(nodeEnd,1)
          changeX = splineNodesData(nodeEnd,1)-splineNodesData(nodeStart,1)
! determine 1 to nodes ZBL, nodes ZBL onwards
          pointsSpline = ceiling((changeX/splineNodesData(nodeEnd,1))*1001) ! data points in spline section
          pointsZBL = 1001 - pointsSpline                                   ! data points in zbl section          
! ZBL section
          Do j=1,pointsZBL
            x = (j-1)*(splineNodesData(nodeEnd,1)/(1.0D0*1001))
            yArray = ZblFull (x, zA, zB)
            eamPoint = eamKey(functionCounter,4)+j-1
            eamData(eamPoint,1) = x
            eamData(eamPoint,2) = yArray(1)
            eamData(eamPoint,3) = yArray(2)
            eamData(eamPoint,4) = yArray(3)
          End Do         
          splineType = 1    ! default to poly
          splineType(1) = 2 ! first segment exp(poly) [3rd order]
! Spline section
          splineDataPoints = SplineNodesV(splineNodesData,pointsSpline,nodeStart,nodeEnd,1001,splineType)
          Do j=1,pointsSpline          
            eamPoint = eamKey(functionCounter,4)+j+pointsZBL-1
            eamData(eamPoint,1) = splineDataPoints(j,1)
            eamData(eamPoint,2) = splineDataPoints(j,2)
            eamData(eamPoint,3) = splineDataPoints(j,3)
            eamData(eamPoint,4) = splineDataPoints(j,4)
          End Do
! Dens + Embe          
        Else        
! Set node start/end points
          nodeStart = splineNodesKey(i,4)
          nodeLength = splineNodesKey(i,5)
          nodeEnd = splineNodesKey(i,6)
! spline between nodes
          splineDataPoints = SplineNodes(splineNodesData,1001,nodeStart,nodeEnd)  
! copy data into eamData array
          Do j=1,1001
            eamPoint = eamKey(functionCounter,4)+j-1
            eamData(eamPoint,1) = splineDataPoints(j,1)
            eamData(eamPoint,2) = splineDataPoints(j,2)
            eamData(eamPoint,3) = splineDataPoints(j,3)
            eamData(eamPoint,4) = splineDataPoints(j,4)
          End Do
        End If  
      End If
      If(functionCounter.eq.eamFunctionCount)Then
        Exit  ! Exit, all functions cycled through
      End If      
    !Do i=1,eamFunctionCount
      !i = 1
      !Do j=eamKey(i,4),eamKey(i,6)
      !  If(mpiProcessID.eq.0)Then
      !    print *,i,j,eamData(j,1),eamData(j,2),eamData(j,3),eamData(j,4)
      !  End If
      !End Do
    !End Do
    End Do     
  End Subroutine setEamSplineOpti
  
  
  
End Module eamOpt  