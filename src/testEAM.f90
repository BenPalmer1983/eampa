Module testEAM

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
  Use readConfig     ! read config file  
  Use makeConfig     ! read config file  
  Use neighbourList
  Use calcEAM  
  Use calcEval  
  Use bulkProperties 

! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private    
! Public Subroutines
  Public :: runTestEAM
  
Contains  
  
  Subroutine runTestEAM()
! Tests the EAM potential
    Implicit None   ! Force declaration of all variables
! Private variables    
! Print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"---------------------------------------------------------------------"
      print *,"Run Test EAM Subroutine"
      print *,"---------------------------------------------------------------------"
    End If
!--------------------------------------------
! Lattice parameters + Bulk Modulus
!--------------------------------------------    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"FCC aLat" 
    End If
! FCC
    Call makeConfigFile(1,2.50D0)      ! Make FCC configuration and calculate alat
    Call readConfigFile()              ! Read in configuration file
    Call clearNeighbourList()          ! Clear existing neighbour list
    Call makeNeighbourList()           ! Make new neighbour list
! Estimate FCC lattice parameter
    Call estimateALat(fccALat, fccEMin, fccVolMin, fccBM) 
    Call outputALat("FCC:", fccALat, fccEMin, fccVolMin, fccBM)
    Call outputALatT("FCC:", fccALat, fccEMin, fccVolMin, fccBM)
! BCC    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"BCC aLat" 
    End If
! Make BCC configuration and calculate alat
    Call makeConfigFile(2,2.50D0)
! Read in configuration file
    Call readConfigFile()
    Call clearNeighbourList() 
    Call makeNeighbourList()    
! Estimate BCC lattice parameter
    Call estimateALat(bccALat, bccEMin, bccVolMin, bccBM)    
    Call outputALat("BCC:", bccALat, bccEMin, bccVolMin, bccBM)    
    Call outputALatT("BCC:", bccALat, bccEMin, bccVolMin, bccBM)   
!--------------------------------------------
! Elastic Constants
!--------------------------------------------   
! FCC
    Call makeConfigFile(1,fccALat)     ! Make with optimised lattice parameter
    Call readConfigFile()              ! Read in configuration file
    Call clearNeighbourList()          ! Clear existing neighbour list
    Call makeNeighbourList()           ! Make new neighbour list    
! Calculate elastic constants    
    Call calcSpecificEC(1, fccALat, fccVolMin, fccBM, fccEC, -1)
! BCC
    Call makeConfigFile(1,bccALat)     ! Make with optimised lattice parameter
    Call readConfigFile()              ! Read in configuration file
    Call clearNeighbourList()          ! Clear existing neighbour list
    Call makeNeighbourList()           ! Make new neighbour list    
! Calculate elastic constants    
    Call calcSpecificEC(1, bccALat, bccVolMin, bccBM, bccEC, -1)


!--------------------------------------------
! Output
!--------------------------------------------   
    Call outputTestingSummary() 
    Call outputTestingSummaryT()

    
  End Subroutine runTestEAM
  
  
  Subroutine estimateALat(aLatMin, eMin, volumeMin, bm)
! Assumes just 1 config - estimates alat
! ONLY CALCULATES ON CONFIGID 1
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: i, j, minPoint, sPoint, ePoint, configID
    Real(kind=DoubleReal) :: aLatIn, aLatTest, configEnergy, volume, bm
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, eMin, aLatNew, aLatMin, volumeMin
    Real(kind=DoubleReal), Dimension(1:65,1:5) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:12,1:2) :: dataPointsFit
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
    Real(kind=DoubleReal), Dimension(1:5) :: coefficientsV
    Real(kind=DoubleReal), Dimension(1:5) :: coefficientsBM
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsDBM
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsDDBM    
! Init variables    
    configID = 1
    dataPoints = -2.1D20
    eMin = 2.1D20
    aLatMin = 0.0d0
! Starting aLat
    aLatIn = configurationsR(configID,1)
! Save starting neighbour list    
    Call saveConfigNL(configID)    
! First attempt
!---------------------------------------------
! loop through lattice parameters
    Do i=1,27
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then   ! Split between multiple processes
! Load original neighbour list    
        Call loadConfigNL(configID)
! Set testing alat
        aLatTest = 2.25D0 + (i * 0.25D0)
        Call changeALat(configID, aLatIn, aLatTest)
! Calculate volume      
        Call calcVol(configID, aLatTest, volume)
! Calculate energy
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)
! Store data points      
        dataPoints(i,1) = aLatTest
        dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,3) = pairEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,4) = embeddingEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,5) = volume
      End If
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPoints)    
    Call M_distDouble2D(dataPoints)
! Store minimum
    Do i=1,27
      If(dataPoints(i,2).lt.eMin)Then
        eMin = dataPoints(i,2)
        aLatMin = dataPoints(i,1)
        minPoint = i
      End If
    End Do
! Output points
    Call outputALatTest(dataPoints)
! Make fitting points array
    sPoint = minPoint - 5
    ePoint = minPoint + 6
    If(sPoint.lt.1)Then
      sPoint = 1
      ePoint = 12
    End If
    If(ePoint.gt.27)Then
      sPoint = 16
      ePoint = 27
    End If
! Lattice parameter
    dataPointsFit = 0.0D0
    j = 0
    Do i=sPoint,ePoint
      j = j + 1
      dataPointsFit(j,1) = dataPoints(i,1)
      dataPointsFit(j,2) = dataPoints(i,2)
    End Do    
! Fit points and find minimum
    aLatMin = MinPolyFit(dataPointsFit,3)   
! Output minimum to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Min alat attempt 1: ",aLatMin
    End If 
! Second attempt
!---------------------------------------------
    eMin = 2.1D20 
    dataPoints = -2.1D20   
    aLatNew = aLatMin
! loop through lattice parameters
    Do i=1,27
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then   ! Split between multiple processes
        j = i-13
! Load original neighbour list    
        Call loadConfigNL(1)
! Set testing alat
        aLatTest = aLatNew + (j * 0.03D0)
        Call changeALat(configID, aLatIn, aLatTest)
! Calculate volume      
        Call calcVol(configID, aLatTest, volume)
! Calculate energy
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)
! Store data points      
        dataPoints(i,1) = aLatTest
        dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,3) = pairEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,4) = embeddingEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,5) = volume/(1.0D0*configurationCoordsKeyG(configID,2))
      End If
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPoints)    
    Call M_distDouble2D(dataPoints)
! Store minimum
    Do i=1,27
      If(dataPoints(i,2).lt.eMin)Then
        eMin = dataPoints(i,2)
        aLatMin = dataPoints(i,1)
        volumeMin = dataPoints(i,5)
      End If      
    End Do    
! Output points
    Call outputALatTest(dataPoints)
! Make fitting points array
    sPoint = minPoint - 5
    ePoint = minPoint + 6
    If(sPoint.lt.1)Then
      sPoint = 1
      ePoint = 12
    End If
    If(ePoint.gt.27)Then
      sPoint = 16
      ePoint = 27
    End If
! Lattice parameter
    dataPointsFit = 0.0D0
    j = 0
    Do i=sPoint,ePoint
      j = j + 1
      dataPointsFit(j,1) = dataPoints(i,1)
      dataPointsFit(j,2) = dataPoints(i,2)
    End Do    
! Fit points and find minimum
    aLatMin = MinPolyFit(dataPointsFit,3)       
    coefficients = PolyFit(dataPointsFit,3)
    eMin = CalcPolynomial(coefficients,aLatMin)    
! Volume
    dataPointsFit = 0.0D0
    j = 0
    Do i=sPoint,ePoint
      j = j + 1
      dataPointsFit(j,1) = dataPoints(i,5)
      dataPointsFit(j,2) = dataPoints(i,2)
    End Do       
    volumeMin = MinPolyFit(dataPointsFit,4)   
    coefficientsV = PolyFit(dataPointsFit,4)
! Output minimum to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Min alat attempt 2: ",aLatMin,eMin,volumeMin
    End If     
! Calculate bulk modulus from second fit data points 
    coefficientsBM = PolyFit(dataPointsFit,4)
    coefficientsDBM = DerivativePolynomial(coefficientsBM)
    coefficientsDDBM = DerivativePolynomial(coefficientsDBM)
    bm = CalcPolynomial (coefficientsDDBM, volumeMin, 0)    
    bm = volumeMin * bm
    bm = UnitConvert(bm,"EVAN3","GPA")     
  End Subroutine estimateALat
  

  
  
  
!---------------------------------------------------------------------------------------------------
    
    
  Subroutine changeALat(configID, aLatOld, aLatNew)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: aLatOld, aLatNew
    Integer(kind=StandardInteger) :: i, keyS, keyE
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
! Isotropic distortion, energy calculation only - just change nl separation
    Do i=keyS,keyE
      neighbourListR(i) = neighbourListR(i)*(aLatNew/(1.0D0*aLatOld))
    End Do
  End Subroutine changeALat  
    
    
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
  
  Subroutine calcVol(configID, aLat, volume)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
    Real(kind=DoubleReal) :: aLat, volume
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: superCellVectors,configUnitVector
    !Integer(kind=StandardInteger) :: i, keyS, keyE
! Init variables
    xCopy = configurationsI(configID,1) 
    yCopy = configurationsI(configID,2) 
    zCopy = configurationsI(configID,3) 
! Config unit vector    
    configUnitVector(1,1) = configurationsR(configID,21) 
    configUnitVector(1,2) = configurationsR(configID,22) 
    configUnitVector(1,3) = configurationsR(configID,23) 
    configUnitVector(2,1) = configurationsR(configID,24) 
    configUnitVector(2,2) = configurationsR(configID,25) 
    configUnitVector(2,3) = configurationsR(configID,26) 
    configUnitVector(3,1) = configurationsR(configID,27) 
    configUnitVector(3,2) = configurationsR(configID,28) 
    configUnitVector(3,3) = configurationsR(configID,29) 
! set the supercell vector    
    superCellVectors = 0.0D0
    superCellVectors(1,1) = 1.0D0*aLat*xCopy
    superCellVectors(2,2) = 1.0D0*aLat*yCopy
    superCellVectors(3,3) = 1.0D0*aLat*zCopy
! apply config unit vector to supercell    
    superCellVectors = MatMul(configUnitVector,superCellVectors)
! calculate volume    
    volume = TripleProductSq(superCellVectors)
  End Subroutine calcVol  
  
  
  
  
  
  
  
  
  
  
End Module testEAM