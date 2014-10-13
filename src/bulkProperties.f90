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
  Public :: calcBP
  Public :: calcEquilibrium
  Public :: calcBM
  Public :: calcSpecificEC
  
Contains
  
    
  Subroutine calcBP(configID)
! Calculate the equilibrium volume/energy
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: i, configID
    Real(kind=DoubleReal) :: eqVol, eqEnergy, eqLattice, bm
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPoints  
    Real(kind=DoubleReal), Dimension(1:5) :: polyCoeffs  
! Title print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,""
      print *,""      
      print *,"================================================================="
      print *,"                   Analysis of EAM Potential                     "      
      print *,"================================================================="
      print *,"" 
    End If      
! Calculate equilibrium   
    Call calcSpecificEquilibrium(configID, eqVol, eqEnergy, eqLattice, dataPoints, polyCoeffs)
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Equilibrium Volume/Energy/Lattice Parameter" 
      print *,"-----------------------------------------------"
      print *,"Data points:" 
      Do i=1,7
        print *,i,dataPoints(i,1),dataPoints(i,2)
      End Do
      print *,"Fitting coefficients:" 
      Do i=1,5
        print "(A3,I1,A1,E12.6)"," x^",(i-1)," ",polyCoeffs(i)
      End Do  
      print *,"In Vol:  ",configVolume(configID)
      print *,"Eq Vol:  ",eqVol
      print *,"Eq Ene:  ",eqEnergy
      print *,"Eq LatF: ",eqLattice
      print *,"Eq Lat:  ",(eqLattice*configurationsR(configID,1))
      print *,"" 
    End If
! Calculate bulk modulus
    Call calcSpecificBM(configID, eqVol, eqLattice, bm, dataPoints, polyCoeffs)
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Bulk Modulus" 
      print *,"-----------------------------------------------"
      print *,"Data points:" 
      Do i=1,5
        print *,i,dataPoints(i,1),dataPoints(i,2)
      End Do
      print *,"Fitting coefficients:" 
      Do i=1,5
        print "(A3,I1,A1,E12.6)"," x^",(i-1)," ",polyCoeffs(i)
      End Do  
      print *,"BM:  ",bm
      print *,"" 
    End If    
    
  End Subroutine calcBP  
    
    
  
  Subroutine calcEquilibrium()
! Calculate the equilibrium volume/energy
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: eqTimeStart, eqTimeEnd
    Real(kind=DoubleReal) :: eqVol, eqEnergy, eqLattice
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPoints  
    Real(kind=DoubleReal), Dimension(1:5) :: polyCoeffs  
! Start time
    Call cpu_time(eqTimeStart)
! Reset calc ev array
    configCalcEV = -2.1D20   
    configCalcEE = -2.1D20  
    configCalcEL = -2.1D20  
! Loop through configurations
    Do configID=1,configCount
! Assign process ID
      If(processMap(configID,2).eq.mpiProcessID)Then
        Call calcSpecificEquilibrium(configID, eqVol, eqEnergy, eqLattice, dataPoints, polyCoeffs)    
        configCalcEV(configID) = eqVol
        configCalcEE(configID) = eqEnergy
        configCalcEL(configID) = eqLattice
      End If  
    End Do
! Distribute eqvol array    
    Call M_collDouble1DMap(configCalcEV,processMap,2)
    Call M_collDouble1DMap(configCalcEE,processMap,2)
    Call M_collDouble1DMap(configCalcEL,processMap,2)
    Call M_distDouble1D(configCalcEV)
    Call M_distDouble1D(configCalcEE)
    Call M_distDouble1D(configCalcEL)
    configID = 1
! End time
    Call cpu_time(eqTimeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("Eq Vol Calc",eqTimeEnd-eqTimeStart)    
  End Subroutine calcEquilibrium   
  
  Subroutine calcSpecificEquilibrium(configID, eqVol, eqEnergy, eqLattice, dataPointsV, polyCoeffs)
! Calculate the equilibrium volume/energy
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, i    
    Real(kind=DoubleReal) :: configEnergy  
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPointsV, dataPointsL  
    Real(kind=DoubleReal), Dimension(1:5) :: polyCoeffs  
    Real(kind=DoubleReal) :: eqVol, eqEnergy, eqLattice, eqLatticeVal, strain
    Real(kind=DoubleReal) :: aLat, calcVol, aLatIn, calcVolIn
! Init variables    
    aLatIn = configurationsR(configID,1)
    calcVolIn  = configVolume(configID)
! Store original neighbour list to reload   
    Call saveConfigNL(configID) 
    Call calcEnergy(configID, configEnergy, 0)
    dataPointsV(4,1) = calcVolIn
    dataPointsV(4,2) = configEnergy
    dataPointsL(4,1) = aLatIn
    dataPointsL(4,2) = configEnergy
! Coarse first attempt
    Do i=-3,3
      If(i.ne.0)Then
        If(i.gt.-3)Then
          Call loadConfigNL(configID)
        End If
        strain = 0.02D0*i
        alat = aLatIn*(1.0D0+strain)
        Call calcVolume(configID,aLat,calcVol)
        Call isotropicDistortion(configID,strain)
        Call calcEnergy(configID, configEnergy, 0)
        dataPointsV(4+i,1) = calcVol
        dataPointsV(4+i,2) = configEnergy  
        dataPointsL(4+i,1) = alat     
        dataPointsL(4+i,2) = configEnergy  
      End If
    End Do
    Call loadConfigNL(configID)
    eqVol = MinPolyFit(dataPointsV,4) 
    eqLatticeVal = MinPolyFit(dataPointsL,4)
    eqLattice = eqLatticeVal/aLatIn
    polyCoeffs = PolyFit(dataPointsV,4)
    eqEnergy = CalcPolynomial(polyCoeffs, eqVol) 
    Call outputEquilibriumPoints(dataPointsL,dataPointsV)  
  End Subroutine calcSpecificEquilibrium 
  
  Subroutine calcBM()
! Calculate the equilibrium volume/energy
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: bmTimeStart, bmTimeEnd
    Real(kind=DoubleReal) :: eqVol, eqLattice
    Real(kind=DoubleReal) :: bm
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPoints  
    Real(kind=DoubleReal), Dimension(1:5) :: polyCoeffs  
! Start time
    Call cpu_time(bmTimeStart)
! Reset calc ev array
    configCalcBM = -2.1D20   
! Loop through configurations
    Do configID=1,configCount
! Assign process ID
      If(processMap(configID,3).eq.mpiProcessID)Then
        eqVol = configCalcEV(configID)
        eqLattice = configCalcEL(configID)      
        Call calcSpecificBM(configID, eqVol, eqLattice, bm, dataPoints, polyCoeffs)        
        configCalcBM(configID) = bm        
      End If  
    End Do
! Distribute eqvol array    
    Call M_collDouble1DMap(configCalcBM,processMap,3)
    Call M_distDouble1D(configCalcBM)
    configID = 1
! End time
    Call cpu_time(bmTimeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("BM Calc",bmTimeEnd-bmTimeStart)    
  End Subroutine calcBM   
  
  
  Subroutine calcSpecificBM(configID, eqVol, eqLattice, bm, dataPoints, polyCoeffs)
! Calculate the equilibrium volume/energy
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, i   
    Real(kind=DoubleReal) :: configEnergy  
    Real(kind=DoubleReal) :: bm, eqVol, eqLattice    
    Real(kind=DoubleReal), Dimension(1:5,1:2) :: dataPoints  
    Real(kind=DoubleReal), Dimension(1:5) :: polyCoeffs  
    Real(kind=DoubleReal), Dimension(1:4) :: polyCoeffsD 
    Real(kind=DoubleReal), Dimension(1:3) :: polyCoeffsDD     
    Real(kind=DoubleReal) :: aLat, calcVol, aLatIn, calcVolIn, strain
! Init variables    
    aLatIn = configurationsR(configID,1)
    calcVolIn  = configVolume(configID)
! Store original neighbour list to reload   
    Call saveConfigNL(configID)     
! 7 Strains
    Do i=-2,2
      If(i.gt.-2)Then
        Call loadConfigNL(configID)
      End If      
      strain = (eqLattice-1.0D0)+0.01D0*i
      alat = aLatIn*(1.0D0+strain)
      Call calcVolume(configID,aLat,calcVol)
      Call isotropicDistortion(configID,strain)
      Call calcEnergy(configID, configEnergy, 0)
      dataPoints(3+i,1) = calcVol
      dataPoints(3+i,2) = configEnergy    
    End Do
    Call loadConfigNL(configID)
    polyCoeffs = PolyFit(dataPoints,4)
    polyCoeffsD = DerivativePolynomial(polyCoeffs)
    polyCoeffsDD = DerivativePolynomial(polyCoeffsD)
    bm = CalcPolynomial (polyCoeffsDD, eqVol, 0)    
    bm = eqVol * bm
    bm = UnitConvert(bm,"EVAN3","GPA")  
  End Subroutine calcSpecificBM 
  
  
    
  Subroutine calcSpecificEC(configID, eqALat, eqVol, bm, elasticConstants)
! Calculate the elastic constants
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger) :: i, j   
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:13,1:2) :: dataPointsD
    Real(kind=DoubleReal), Dimension(1:3) :: polyCoeffs  
    Real(kind=DoubleReal), Dimension(1:10) :: elasticConstants
    Real(kind=DoubleReal) :: configEnergy, eqVol, bm, scVol 
    Real(kind=DoubleReal) :: aLatTemp, eqALat, aLatFactor    
    Real(kind=DoubleReal) :: tetraQuadTerm, orthoQuadTerm, monoQuadTerm
    Real(kind=DoubleReal) :: cAA, cAB, cDD, cAA_AB
! Init variables    
    scVol = eqVol * configurationCoordsKeyG(configID,2)
    elasticConstants = -2.1D20
! Store original neighbour list to reload   
    Call saveConfigNL(configID)         
! Store original lattice parameter to reload
    aLatTemp = configurationsR(configID,1)
    aLatFactor = eqALat/aLatTemp
! Apply strains - Tetragonal
    j = 0
    Do i=-3,3
      If(i.gt.-3)Then
        Call loadConfigNL(configID)
      End If 
      j = j + 1
      Call tetragonalDistortion(configID, (i*0.001D0), aLatFactor)
      Call calcEnergy(configID, configEnergy, 0)
      dataPoints(j,1) = i*0.001D0
      dataPoints(j,2) = configEnergy      
      !If(mpiProcessID.eq.0)Then
      !  print *,j,dataPoints(j,1),dataPoints(j,2)
      !End If
    End Do  
    Call loadConfigNL(configID)
! Fit polynomial
    polyCoeffs = PolyFit(dataPoints,2)
    tetraQuadTerm = polyCoeffs(3)
! Apply strains - Orthorhombic
    j = 0
    Do i=0,6
      If(i.gt.0)Then
        Call loadConfigNL(configID)
      End If 
      j = j + 1
      Call orthogonalDistortion(configID, (i*0.001D0), aLatFactor)
      Call calcEnergy(configID, configEnergy, 0)
      dataPoints(j,1) = i*0.001D0
      dataPoints(j,2) = configEnergy    
    End Do  
    Call loadConfigNL(configID)
! Fit polynomial
    Do i=1,7
      If(i.eq.1)Then
        dataPointsD(7,1) = dataPoints(1,1)
        dataPointsD(7,2) = dataPoints(1,2)
      Else  
        dataPointsD(6+i,1) = dataPoints(i,1)
        dataPointsD(8-i,1) = -1.0D0*dataPoints(i,1)
        dataPointsD(6+i,2) = dataPoints(i,2)
        dataPointsD(8-i,2) = dataPoints(i,2)
      End If
    End Do
    !Do i=1,13
    !  If(mpiProcessID.eq.0)Then
    !    print *,i,dataPointsD(i,1),dataPointsD(i,2)
    !  End If
    !End Do
    polyCoeffs = PolyFit(dataPointsD,2)
    orthoQuadTerm = polyCoeffs(3)
! Apply strains - Monoclinic
    j = 0
    Do i=0,6
      If(i.gt.0)Then
        Call loadConfigNL(configID)
      End If 
      j = j + 1
      Call monoclinicDistortion(configID, (i*0.001D0), aLatFactor)
      Call calcEnergy(configID, configEnergy, 0)
      dataPoints(j,1) = i*0.001D0
      dataPoints(j,2) = configEnergy    
    End Do  
    Call loadConfigNL(configID)
! Fit polynomial
    Do i=1,7
      If(i.eq.1)Then
        dataPointsD(7,1) = dataPoints(1,1)
        dataPointsD(7,2) = dataPoints(1,2)
      Else  
        dataPointsD(6+i,1) = dataPoints(i,1)
        dataPointsD(8-i,1) = -1.0D0*dataPoints(i,1)
        dataPointsD(6+i,2) = dataPoints(i,2)
        dataPointsD(8-i,2) = dataPoints(i,2)
      End If
    End Do
    !Do i=1,13
      !If(mpiProcessID.eq.0)Then
      !  print *,i,dataPointsD(i,1),dataPointsD(i,2)
      !End If
    !End Do
    polyCoeffs = PolyFit(dataPoints,2)
    monoQuadTerm = polyCoeffs(3)
! Calculate elastic constants
    cAA_AB = (tetraQuadTerm/(3*scVol)+orthoQuadTerm/(scVol))/2
    cAA_AB = UnitConvert(cAA_AB,"EVAN3","GPA") 
    cAB = (3.0D0*bm-cAA_AB)/3.0D0
    cAA = 3.0D0*bm-2.0D0*cAB
    cDD = (2*monoQuadTerm)/scVol 
    cDD = UnitConvert(cDD,"EVAN3","GPA")                          !C44   
! Store results
    elasticConstants(1) = cAA 
    elasticConstants(2) = cAB 
    elasticConstants(3) = cDD 
  End Subroutine calcSpecificEC 
  
  
  
  
  
!------------------------------------------------
! Distortion + NL subroutines
!------------------------------------------------   
  
  Subroutine isotropicDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion    
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input    
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn 
    End If
! set distortion matrix
    distortion(1,1) = factor*(1.0D0+distortionAmount)
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = factor*(1.0D0+distortionAmount)
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*(1.0D0+distortionAmount)
    Call applyDistortionNL(configID, distortion, 1)  
  End Subroutine isotropicDistortion
  
  Subroutine tetragonalDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion    
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input    
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn 
    End If
! set distortion matrix
    distortion(1,1) = factor*(1.0D0+distortionAmount)
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = factor*(1.0D0+distortionAmount)
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*((1.0D0+distortionAmount)**(-2.0D0))
    Call applyDistortionNL(configID, distortion, 1)  
  End Subroutine tetragonalDistortion
  
  Subroutine orthogonalDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion    
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input    
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn 
    End If
! set distortion matrix
    distortion(1,1) = factor*(1.0D0+distortionAmount)
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = factor*(1.0D0-distortionAmount)
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*(1.0D0+(distortionAmount**2)/(1.0D0-distortionAmount**(2.0D0)))
    Call applyDistortionNL(configID, distortion, 1)  
  End Subroutine orthogonalDistortion
  
  Subroutine monoclinicDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID 
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion    
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input    
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn 
    End If
! set distortion matrix
    distortion(1,1) = factor*1.0D0
    distortion(1,2) = factor*0.5D0*distortionAmount
    distortion(1,3) = 0.0D0
    distortion(2,1) = factor*0.5D0*distortionAmount
    distortion(2,2) = factor*1.0D0
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*(1.0D0+(distortionAmount**2)/(4.0D0-distortionAmount**(2.0D0)))
    Call applyDistortionNL(configID, distortion, 1)  
  End Subroutine monoclinicDistortion
  
  
  
  Subroutine applyDistortionNL(configID, distortion, typeCalc)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, typeCalc, distOption
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Integer(kind=StandardInteger) :: i, j, arrayCheck, keyS, keyE
    Real(kind=DoubleReal) :: xCA, yCA, zCA, xCB, yCB, zCB
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
      Do i=keyS,keyE
! Apply distortion vector - Point A
        xCA = neighbourListCoords(i,1)*distortion(1,1)+&
              neighbourListCoords(i,2)*distortion(1,2)+&
              neighbourListCoords(i,3)*distortion(1,3)
        yCA = neighbourListCoords(i,1)*distortion(2,1)+&
              neighbourListCoords(i,2)*distortion(2,2)+&
              neighbourListCoords(i,3)*distortion(2,3)
        zCA = neighbourListCoords(i,1)*distortion(3,1)+&
              neighbourListCoords(i,2)*distortion(3,2)+&
              neighbourListCoords(i,3)*distortion(3,3)
        neighbourListCoords(i,1) = xCA     
        neighbourListCoords(i,2) = yCA      
        neighbourListCoords(i,3) = zCA       
! Apply distortion vector - Point A
        xCB = neighbourListCoords(i,4)*distortion(1,1)+&
              neighbourListCoords(i,5)*distortion(1,2)+&
              neighbourListCoords(i,6)*distortion(1,3)
        yCB = neighbourListCoords(i,4)*distortion(2,1)+&
              neighbourListCoords(i,5)*distortion(2,2)+&
              neighbourListCoords(i,6)*distortion(2,3)
        zCB = neighbourListCoords(i,4)*distortion(3,1)+&
              neighbourListCoords(i,5)*distortion(3,2)+&
              neighbourListCoords(i,6)*distortion(3,3)
        neighbourListCoords(i,4) = xCB     
        neighbourListCoords(i,5) = yCB      
        neighbourListCoords(i,6) = zCB      
! Store differences        
        neighbourListCoords(i,7) = neighbourListCoords(i,4) - neighbourListCoords(i,1)    
        neighbourListCoords(i,8) = neighbourListCoords(i,5) - neighbourListCoords(i,2)   
        neighbourListCoords(i,9) = neighbourListCoords(i,6) - neighbourListCoords(i,3)   
! Store separation
        neighbourListR(i) = (neighbourListCoords(i,7)**2 + &
                            neighbourListCoords(i,8)**2 + &
                            neighbourListCoords(i,9)**2)**0.5
! Store directon magnitudes
        neighbourListCoords(i,10) = neighbourListCoords(i,7)/neighbourListR(i)
        neighbourListCoords(i,11) = neighbourListCoords(i,8)/neighbourListR(i)
        neighbourListCoords(i,12) = neighbourListCoords(i,9)/neighbourListR(i)
      End Do
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
  
  
  Subroutine calcVolume(configID,aLat,volume)
! Load NL for one config from the temp memory array
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, xCopy, yCopy, zCopy, j
    Real(kind=DoubleReal) :: aLat, volume
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: configUnitVector 
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: configVolVector 
! Init variables
    volume = 0.0D0
    xCopy = configurationsI(configID,1) 
    yCopy = configurationsI(configID,2) 
    zCopy = configurationsI(configID,3) 
! Get config unit vector
    configUnitVector(1,1) = configurationsR(configID,2) 
    configUnitVector(1,2) = configurationsR(configID,3) 
    configUnitVector(1,3) = configurationsR(configID,4) 
    configUnitVector(2,1) = configurationsR(configID,5) 
    configUnitVector(2,2) = configurationsR(configID,6) 
    configUnitVector(2,3) = configurationsR(configID,7) 
    configUnitVector(3,1) = configurationsR(configID,8) 
    configUnitVector(3,2) = configurationsR(configID,9) 
    configUnitVector(3,3) = configurationsR(configID,10) 
! Apply global unit vector      
    configUnitVector = matmul(globalConfigUnitVector,configUnitVector)
! Make volume vector
    Do j=1,3 
      configVolVector(j,1) = aLat*xCopy*configUnitVector(j,1)
      configVolVector(j,2) = aLat*yCopy*configUnitVector(j,2)
      configVolVector(j,3) = aLat*zCopy*configUnitVector(j,3)
    End Do  
    volume = TripleProductSq(configVolVector)
  End Subroutine calcVolume
  
  
End Module bulkProperties  