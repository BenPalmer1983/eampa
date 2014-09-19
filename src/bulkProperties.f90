Module bulkProperties

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Read user input file 

!----------------------------------------
! Updated: 12th Aug 2014
!----------------------------------------

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
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: calcEquilibrium
  
Contains
  
  
  
  Subroutine calcEquilibrium()
! Calculate the equilibrium volume/energy
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, i    
    Real(kind=DoubleReal) :: configEnergy  
    Real(kind=DoubleReal), Dimension(1:5,1:2) :: dataPoints    
    Real(kind=DoubleReal) :: eqTimeStart, eqTimeEnd
    Real(kind=DoubleReal), Dimension(1:1024,1:10) :: outputEqVol
    Real(kind=DoubleReal), Dimension(1:4) :: polyCoeffs
! Start time
    Call cpu_time(eqTimeStart)
! Reset calc ev array
    configCalcEV = -2.1D20   
! Blank output array    
    outputEqVol = -2.1D20 
! Loop through configurations
    Do configID=1,configCount
! Assign process ID
      If(processMap(configID,2).eq.mpiProcessID)Then
! Store original neighbour list to reload   
        Call saveConfigNL(configID) 
        Call calcEnergy(configID, configEnergy, 0)
        dataPoints(3,1) = configVolume(configID)
        dataPoints(3,2) = configEnergy
        Do i=-2,2
          If(i.ne.0)Then
            If(i.gt.-2)Then
              Call loadConfigNL(configID)
            End If
            Call isotropicDistortion(configID,0.05D0*i)
            Call calcEnergy(configID, configEnergy, 0)
            dataPoints(3+i,1) = ((1.0D0+0.05D0*i)**3)*configVolume(configID)
            dataPoints(3+i,2) = configEnergy
          End If
        End Do  
        Do i=1,5
          outputEqVol(configID,2*(i-1)+1) = dataPoints(i,1)
          outputEqVol(configID,2*(i-1)+2) = dataPoints(i,2)
        End Do        
        Call loadConfigNL(configID)   
! Find optimum energy/volume
        configCalcEV(configID) = MinPolyFit(dataPoints,3) 
        polyCoeffs = PolyFit(dataPoints,3)
        configCalcEE(configID) = CalcPolynomial(polyCoeffs, configCalcEV(configID))
      End If  
    End Do
! Distribute eqvol array    
    Call M_collDouble1DMap(configCalcEV,processMap,2)
    Call M_collDouble1DMap(configCalcEE,processMap,2)
    Call M_distDouble1D(configCalcEV)
    Call M_distDouble1D(configCalcEE)
! End time
    Call cpu_time(eqTimeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("Eq Vol Calc",eqTimeEnd-eqTimeStart)    
  End Subroutine calcEquilibrium   
  
  
  
  
  
  
  
!------------------------------------------------
! Distortion + NL subroutines
!------------------------------------------------ 
  
  
  Subroutine isotropicDistortion(configID, distortionAmount)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: distortionAmount
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion    
! set distortion matrix
    distortion(1,1) = 1.0D0+distortionAmount
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = 1.0D0+distortionAmount
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = 1.0D0+distortionAmount
    Call applyDistortionNL(configID, distortion, 1)  
  End Subroutine isotropicDistortion
  
  Subroutine applyDistortionNL(configID, distortion, typeCalc)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, typeCalc, distOption
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Integer(kind=StandardInteger) :: i, j, arrayCheck, keyS, keyE
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
    distOption = 2 ! Full nl modification for distortion
! Check type of distortion (just energy, isotropic distortion)
    If(typeCalc.eq.1)Then  ! Energy only calculation
      ! Check distortion array
      If(distortion(1,1).eq.distortion(2,2).and.distortion(1,1).eq.distortion(3,3))Then
        arrayCheck = 0
        Do i=1,3 
          Do j=1,3 
            If(j.ne.i)Then
              If(distortion(i,j).ne.0.0D0)Then
                arrayCheck = 1
              End If
            End If
          End Do
        End Do
        If(arrayCheck.eq.0)Then
          distOption = 1
        End If
      End If
    End If
! Distort neighbour list  
    If(distOption.eq.1)Then
! Isotropic distortion, energy calculation only - just change nl separation
      Do i=keyS,keyE
        neighbourListR(i) = distortion(1,1)*neighbourListR(i)
      End Do
    Else
! Non isotropic/energy-force-stress calculations         
    End If
  End Subroutine applyDistortionNL
  
  Subroutine saveConfigNL(configID)
! Save NL for one config to the temp memory array
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, keyS, keyE, i, j
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
! copy data
    Do j=1,size(neighbourListI,2)
      Do i=keyS,keyE
        neighbourListIT(i,j) = neighbourListI(i,j)
      End Do
    End Do  
    Do i=keyS,keyE
      neighbourListRT(i) = neighbourListR(i)
    End Do     
    Do j=1,size(neighbourListCoords,2)
      Do i=keyS,keyE
        neighbourListCoordsT(i,j) = neighbourListCoords(i,j)
      End Do
    End Do 
  End Subroutine saveConfigNL
  
  Subroutine loadConfigNL(configID)
! Load NL for one config from the temp memory array
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, keyS, keyE, i, j
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
! copy data
    Do j=1,size(neighbourListI,2)
      Do i=keyS,keyE
        neighbourListI(i,j) = neighbourListIT(i,j)
      End Do
    End Do  
    Do i=keyS,keyE
      neighbourListR(i) = neighbourListRT(i)
    End Do     
    Do j=1,size(neighbourListCoords,2)
      Do i=keyS,keyE
        neighbourListCoords(i,j) = neighbourListCoordsT(i,j)
      End Do
    End Do 
  End Subroutine loadConfigNL
  
  
End Module bulkProperties  