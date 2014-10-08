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
    Real(kind=DoubleReal) :: aLatMin, eMin
! Print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"---------------------------------------------------------------------"
      print *,"Run Test EAM Subroutine"
      print *,"---------------------------------------------------------------------"
    End If
!--------------------------------------------
! Lattice parameters
!--------------------------------------------    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"FCC aLat" 
    End If
! FCC
! Make FCC configuration and calculate alat
    Call makeConfigFile(1)
! Read in configuration file
    Call readConfigFile()
    Call clearNeighbourList() 
    Call makeNeighbourList()    
! Estimate FCC lattice parameter
    Call estimateALat(aLatMin, eMin)    
! Save
    fccALat = aLatMin
    fccEMin = eMin
    Call outputALat("FCC:", fccALat, fccEMin)
! BCC    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"BCC aLat" 
    End If
! Make BCC configuration and calculate alat
    Call makeConfigFile(2)
! Read in configuration file
    Call readConfigFile()
    Call clearNeighbourList() 
    Call makeNeighbourList()    
! Estimate BCC lattice parameter
    Call estimateALat(aLatMin, eMin)    
! Save
    bccALat = aLatMin  
    bccEMin = eMin 
    Call outputALat("BCC:", bccALat, bccEMin)     
!--------------------------------------------
! Bulk modulus
!--------------------------------------------   
! Update lattice parameter for FCC   
    configurationsR(1,1) = fccALat
    Call makeConfigFile(1)  ! Make FCC
    Call readConfigFile()
    Call clearNeighbourList() 
    Call makeNeighbourList()    

    
! Update lattice parameter for BCC   
    configurationsR(1,1) = bccALat
    Call makeConfigFile(2)  ! Make BCC
    Call readConfigFile()
    Call clearNeighbourList() 
    Call makeNeighbourList()    

    
    
    
    
  
! Clear neighbour list
! Make neighbour list   

    !If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
    !  print *,"Total Spline Nodes: ",splineTotalNodes
    !End If
! Start Time
    !Call cpu_time(timeStartOpt)

    
    
  End Subroutine runTestEAM
  
  
  Subroutine estimateALat(aLatMin, eMin)
! Assumes just 1 config - estimates alat
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: i, j, minPoint, sPoint, ePoint
    Real(kind=DoubleReal) :: aLatIn, aLatTest, configEnergy
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, eMin, aLatNew, aLatMin
    Real(kind=DoubleReal), Dimension(1:65,1:4) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:12,1:2) :: dataPointsFit
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients
! Init variables    
    dataPoints = -2.1D20
    eMin = 2.1D20
    aLatMin = 0.0d0
! Starting aLat
    aLatIn = configurationsR(1,1)
! Save starting neighbour list    
    Call saveConfigNL(1)    
! First attempt
! loop through lattice parameters
    Do i=1,65
! Load original neighbour list    
      Call loadConfigNL(1)
! Set testing alat
      aLatTest = 2.4D0 + (i * 0.1D0)
      Call changeALat(1, aLatIn, aLatTest)
! Calculate energy
      Call calcEnergy(1, configEnergy, 0, pairEnergy, embeddingEnergy)
! Store data points      
      dataPoints(i,1) = aLatTest
      dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(1,2))
      dataPoints(i,3) = pairEnergy/(1.0D0*configurationCoordsKeyG(1,2))
      dataPoints(i,4) = embeddingEnergy/(1.0D0*configurationCoordsKeyG(1,2))
! Store minimum
      If(dataPoints(i,2).lt.eMin)Then
        eMin = dataPoints(i,2)
        aLatMin = dataPoints(i,1)
      End If
    End Do
! Output points
    Call outputALatTest(dataPoints)
! Output minimum to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Min alat attempt 1: ",aLatMin
    End If 
! Second attempt
    eMin = 2.1D20 
    dataPoints = -2.1D20   
    aLatNew = aLatMin
! loop through lattice parameters
    Do i=1,30
      j = i-15
! Load original neighbour list    
      Call loadConfigNL(1)
! Set testing alat
      aLatTest = aLatNew + (j * 0.01D0)
      Call changeALat(1, aLatIn, aLatTest)
! Calculate energy
      Call calcEnergy(1, configEnergy, 0, pairEnergy, embeddingEnergy)
! Store data points      
      dataPoints(i,1) = aLatTest
      dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(1,2))
      dataPoints(i,3) = pairEnergy/(1.0D0*configurationCoordsKeyG(1,2))
      dataPoints(i,4) = embeddingEnergy/(1.0D0*configurationCoordsKeyG(1,2))
! Store minimum
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
    If(ePoint.gt.30)Then
      sPoint = 19
      ePoint = 30
    End If
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
! Output minimum to terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Min alat attempt 2: ",aLatMin,eMin
    End If     
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
  
  
End Module testEAM