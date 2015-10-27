Module makePotential
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
  Use globals
  Use initialise
  Use loadData
  Use output
  Use readEAM
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: makeEAMRun
  Public :: makeEAM
  Contains
! ---------------------------------------------------------------------------------------------------


  Subroutine makeEAMRun(pairFactorIn, densityFactorIn, embedFactorIn)
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! In    
    Real(kind=DoubleReal), Optional :: pairFactorIn, densityFactorIn, embedFactorIn
    Real(kind=DoubleReal) :: pairFactor, densityFactor, embedFactor
! Private variables    
    Character(len=16), Dimension(1:1) :: labelList
! Optional arguments
    pairFactor = 1.0D0
    densityFactor = 1.0D0
    embedFactor = 1.0D0
    If(Present(pairFactorIn))Then
      pairFactor = pairFactorIn
    End If
    If(Present(densityFactorIn))Then
      densityFactor = densityFactorIn
    End If
    If(Present(embedFactorIn))Then
      embedFactor = embedFactorIn
    End If
    
    !pointsPerFunction = 1001
    labelList(1) = "PD              "
    !labelList(2) = "FE              "
    !labelList(3) = "CR              "
  
  
    print *,"Make EAM Potential"
    
    
    Call makeEAM(labelList, 1001, pairFactor, densityFactor, embedFactor) 
    
  End Subroutine makeEAMRun


  Subroutine makeEAM(labelList, pointsPerFunction, pairFactor, densityFactor, embedFactor)  
! Run simulated annealing optimisation 
    Implicit None   ! Force declaration of all variables
! In
    Character(len=16), Dimension(:) :: labelList
    Integer(kind=StandardInteger) :: pointsPerFunction
    Real(kind=DoubleReal) :: pairFactor, densityFactor, embedFactor
! Private variables        
    Integer(kind=StandardInteger) :: i, j, key, n, keyStart, keyEnd, dataPoint
    Real(kind=DoubleReal) :: zblCutoffRadius, radiusMax, densityMax, x
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
! Initial settings - could be added to input file in future
    zblCutoffRadius = 1.5D0
    radiusMax = 6.5D0
    densityMax = 1.0D0
! Init keys
    key = 0
    keyStart = 0
    dataPoint = 0
! eamKey   1 atomA, 2 atomB, 3 function/al type, 4 func start, 5 func length, 6 func end
! eamData  1 x, 2 y(x), 3 y'(x), 4 y''(x)    
! Clear eam array
    eamKey = 0
    eamData = 0.0D0       
! Pair    
    Do i=1,size(labelList,1)
      Do j=i,size(labelList,1)
        key = key + 1     
        keyStart = keyEnd + 1        
        Do n=1,pointsPerFunction
          x = radiusMax*((n-1.0D0)/(pointsPerFunction-1.0D0))
! Use a ZBL with a sine wave as a starting point          
          yArray = ZblFull (x, 46, 46)
          yArray(1) = pairFactor*(yArray(1) + ((5.0D-3)*sin(2.0D0*x)))
          yArray(2) = pairFactor*(yArray(2) + (2.0D0*(5.0D-3)*cos(2.0D0*x)))
          yArray(3) = pairFactor*(yArray(3) - (4.0D0*(5.0D-3)*sin(2.0D0*x)))
          dataPoint = dataPoint + 1
! Store data point
          eamData(dataPoint,1) = x
          eamData(dataPoint,2) = yArray(1)
          eamData(dataPoint,3) = yArray(2)
          eamData(dataPoint,4) = yArray(3)
        End Do
        keyEnd = dataPoint
! Store key        
        eamKey(key,1) = i
        eamKey(key,2) = j
        eamKey(key,3) = 1
        eamKey(key,4) = keyStart
        eamKey(key,5) = pointsPerFunction
        eamKey(key,6) = keyEnd
      End Do
    End Do    
! Dens
    Do i=1,size(labelList,1)
      key = key + 1
      keyStart = keyEnd + 1  
      Do n=1,pointsPerFunction
        x = radiusMax*((n-1.0D0)/(pointsPerFunction-1.0D0))
        yArray(1) =  2.5D0*x**2*exp(-1.0D0*x**2)&
                    +2.5D1*x**2*exp(-1.0D1*x**2)        
        yArray(2) =  5.0D0*x**1*exp(-1.0D0*x**2)&
                    +2.5D0*x**3*(-1.0D0*2)*exp(-1.0D0*x**2)&
                    +5.0D1*x**1*exp(-1.0D1*x**2)&
                    +2.5D1*x**3*(-1.0D1*2)*exp(-1.0D1*x**2)        
        yArray(3) =  5.0D0*exp(-1.0D0*x**2)&
                    +5.0D0*x**1*(-1.0D0*2*x)*exp(-1.0D0*x**2)&                    
                    +7.5D0*x**2*(-1.0D0*2)*exp(-1.0D0*x**2)&
                    +2.5D0*x**3*(-1.0D0*2)*(-1.0D0*x*2)*exp(-1.0D0*x**2)&
                    +5.0D1*exp(-1.0D1*x**2)&
                    +5.0D1*x**1*(-1.0D1*x*2)*exp(-1.0D1*x**2)&
                    +7.5D1*x**2*(-1.0D1*2)*exp(-1.0D1*x**2)&
                    +2.5D1*x**3*(-1.0D1*2)*(-1.0D1*x*2)*exp(-1.0D1*x**2)
        dataPoint = dataPoint + 1
! Store data point
        eamData(dataPoint,1) = x
        eamData(dataPoint,2) = densityFactor*yArray(1)
        eamData(dataPoint,3) = densityFactor*yArray(2)
        eamData(dataPoint,4) = densityFactor*yArray(3)
      End Do
      keyEnd = dataPoint
! Store key        
      eamKey(key,1) = i
      eamKey(key,2) = j
      eamKey(key,3) = 2
      eamKey(key,4) = keyStart
      eamKey(key,5) = pointsPerFunction
      eamKey(key,6) = keyEnd
    End Do
    
    
    
! Embe
    Do i=1,size(labelList,1)
      key = key + 1
      keyStart = keyEnd + 1  
      Do n=1,pointsPerFunction
        x = densityMax*((n-1.0D0)/(pointsPerFunction-1.0D0))
        yArray(1) = 0.5D0-4.0D0*x**0.5D0+x**2.0D0+1.5D0*x**4.0D0
        yArray(2) = -2.0D0*x**(-0.5D0)+2.0D0*x+6.0D0*x**3.0D0
        yArray(3) = 1.0D0*x**(-1.5D0)+2.0D0+18.0D0*x**2.0D0
        dataPoint = dataPoint + 1
! Store data point
        eamData(dataPoint,1) = x
        eamData(dataPoint,2) = embedFactor*yArray(1)
        eamData(dataPoint,3) = embedFactor*yArray(2)
        eamData(dataPoint,4) = embedFactor*yArray(3)
      End Do
      keyEnd = dataPoint
! Store key        
      eamKey(key,1) = i
      eamKey(key,2) = j
      eamKey(key,3) = 3
      eamKey(key,4) = keyStart
      eamKey(key,5) = pointsPerFunction
      eamKey(key,6) = keyEnd
    End Do
    
    
  End Subroutine makeEAM
  
  

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

End Module makePotential  
  