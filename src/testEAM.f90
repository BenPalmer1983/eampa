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
  Use bulkProperties 

! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private    
! 
  Real(kind=DoubleReal) :: rssMurnFCC
  Real(kind=DoubleReal) :: rssBirchMurnFCC
  Real(kind=DoubleReal) :: rssMurnBCC
  Real(kind=DoubleReal) :: rssBirchMurnBCC
  
  !Integer(kind=StandardInteger) :: neighbourListCountTT, configID
  !Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: neighbourListKeyTT
  !Real(kind=DoubleReal), Dimension(1:1024,1:1) :: neighbourListKeyRTT
  !Integer(kind=StandardInteger), Dimension(1:800000,1:6) :: neighbourListITT 
  !Real(kind=DoubleReal), Dimension(1:800000) :: neighbourListRTT             
  !Real(kind=DoubleReal), Dimension(1:800000,12) :: neighbourListCoordsTT
  
! Public Subroutines
  Public :: runTestAnalysis
  Public :: runTestEAM
  Real(kind=DoubleReal) :: tempRSS
  
Contains 

  Subroutine runTestAnalysis()
! Tests the EAM potential
    Implicit None   ! Force declaration of all variables
! Private variables    
    Call runTestEAM(0) 
! Title print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,""
      print *,""      
      print *,"================================================================="
      print *,"                   Analysis of EAM Potential                     "      
      print *,"================================================================="
      print *,"" 
      print *,"Murnaghan Fit - Alat, V0, E0, B0, B0'"
      print *,"-----------------------------------------------------" 
      print *,"Fit RSS:                 ",rssMurnFCC
      print *,"FCC alat opt:            ",fccALatMurn,"       (",&
              fccReferenceValues(1),")    [",(fccALatMurn-fccReferenceValues(1))**2,"]"
      print *,"FCC vol opt:             ",fccVolMinMurn,"       (",(fccReferenceValues(1)**3/4),")"
      print *,"FCC energy min:          ",fccEMinMurn,"       (",&
              fccReferenceValues(2),")    [",(fccEMinMurn-fccReferenceValues(2))**2,"]"
      print *,"FCC bulk modulus:        ",fccBMMurn,"       (",&
              fccReferenceValues(3),")    [",(fccBMMurn-fccReferenceValues(3))**2,"]"
      print *,"FCC b0':                 ",fccBMPMurn
      tempRSS = (fccALatMurn-fccReferenceValues(1))**2+(fccEMinMurn-fccReferenceValues(2))**2+&
                (fccBMMurn-fccReferenceValues(3))**2
      print *,"RSS:                    ",tempRSS
      print *,"-----------------------------------------------------" 
      print *,"Fit RSS:                 ",rssMurnBCC
      print *,"BCC alat opt:            ",bccALatMurn
      print *,"BCC vol opt:             ",bccVolMinMurn
      print *,"BCC energy min:          ",bccEMinMurn
      print *,"BCC bulk modulus:        ",bccBMMurn
      print *,"BCC b0':                 ",bccBMPMurn
      print *,"" 
      print *,"Birch-Murnaghan - Alat, V0, E0, B0, B0'"
      print *,"-----------------------------------------------------" 
      print *,"Fit RSS:                 ",rssBirchMurnFCC
      print *,"FCC alat opt:            ",fccALatBirchMurn,"       (",&
              fccReferenceValues(1),")    [",(fccALatBirchMurn-fccReferenceValues(1))**2,"]"
      print *,"FCC vol opt:             ",fccVolMinBirchMurn,"       (",&
              (fccReferenceValues(1)**3/4),")"
      print *,"FCC energy min:          ",fccEMinBirchMurn,"       (",&
              fccReferenceValues(2),")    [",(fccEMinBirchMurn-fccReferenceValues(2))**2,"]"
      print *,"FCC bulk modulus:        ",fccBMBirchMurn,"       (",&
              fccReferenceValues(3),")    [",(fccBMBirchMurn-fccReferenceValues(3))**2,"]"
      print *,"FCC b0':                 ",fccBMPBirchMurn
      tempRSS = (fccALatBirchMurn-fccReferenceValues(1))**2+(fccEMinBirchMurn-fccReferenceValues(2))**2+&
                (fccBMBirchMurn-fccReferenceValues(3))**2
      print *,"RSS:                    ",tempRSS
      print *,"-----------------------------------------------------" 
      print *,"Fit RSS:                 ",rssBirchMurnBCC
      print *,"BCC alat opt:            ",bccALatBirchMurn
      print *,"BCC vol opt:             ",bccVolMinBirchMurn
      print *,"BCC energy min:          ",bccEMinBirchMurn
      print *,"BCC bulk modulus:        ",bccBMBirchMurn
      print *,"BCC b0':                 ",bccBMPBirchMurn
      print *,"" 
      If(testingFitChoice.eq.1)Then
        print *,"Murn Fit Elastic Constants" 
      End If 
      If(testingFitChoice.eq.2)Then
        print *,"Birch-Murn Fit Elastic Constants" 
      End If
      print *,"-----------------------------------------------------" 
      print *,"FCC C11:                 ",fccEC(1),"       (",fccReferenceValues(4),")"
      print *,"FCC C12:                 ",fccEC(2),"       (",fccReferenceValues(5),")"
      print *,"FCC C44:                 ",fccEC(3),"       (",fccReferenceValues(6),")"
      print *,"-----------------------------------------------------" 
      print *,"BCC C11:                 ",bccEC(1)
      print *,"BCC C12:                 ",bccEC(2)
      print *,"BCC C44:                 ",bccEC(3)
      print *,"-----------------------------------------------------"       
      print *,""   
      print *,"Testing-Reference RSS:    ",testingRSS
      print *,"" 
      print *,"================================================================="
      
      
      
      
    End If   


  End Subroutine runTestAnalysis  
 
  
  Subroutine runTestEAM(printTestingIn)
! Tests the EAM potential
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: i, configID
    Integer(kind=StandardInteger), optional :: printTestingIn
    Integer(kind=StandardInteger) :: printTesting
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsMurn 
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBirchMurn 
    Real(kind=DoubleReal) :: rssMurnFit,rssBirchMurnFit
! Optional variables
    printTesting = 0
    If(present(printTestingIn))Then
      printTesting = printTestingIn
    End If
! Print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      print *,"---------------------------------------------------------------------"
      print *,"Run Test EAM Subroutine"
      print *,"---------------------------------------------------------------------"
    End If
! Set ID slot to use for temporary configs/neighbour lists
    configID = 901
!--------------------------------------------
! Lattice parameters + Bulk Modulus
!--------------------------------------------    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      print *,"FCC aLat" 
    End If
! FCC
    Call makeConfigFile(1,2.50D0)      ! Make FCC configuration and calculate alat
    Call readConfigFile(0,configID)           ! Read in configuration file into config slot 901
    Call makeNeighbourList(configID,configID,printTesting)           ! Make new neighbour list
! Estimate FCC lattice parameter
    !Call estimateALat(configID,fccALat, fccEMin, fccVolMin, fccBM, &
    !                 coefficientsMurn, coefficientsBirchMurn, &
    !                 rssMurnFit, rssBirchMurnFit) 
    Call estimateALat(configID,&
                     coefficientsMurn, coefficientsBirchMurn, &
                     rssMurnFit, rssBirchMurnFit) 
! Murnaghan fit    
    fccEMinMurn = coefficientsMurn(1)
    fccVolMinMurn = coefficientsMurn(2)
    fccALatMurn = (fccVolMinMurn*4)**(1.0D0/3.0D0)
    fccBMMurn = coefficientsMurn(3)*160.2176568    
    fccBMPMurn = coefficientsMurn(4)
    rssMurnFCC = rssMurnFit
! Birch-Murnaghan fit    
    fccEMinBirchMurn = coefficientsBirchMurn(1)
    fccVolMinBirchMurn = coefficientsBirchMurn(2)
    fccALatBirchMurn = (fccVolMinBirchMurn*4)**(1.0D0/3.0D0)
    fccBMBirchMurn = coefficientsBirchMurn(3)*160.2176568  
    fccBMPBirchMurn = coefficientsBirchMurn(4)
    rssBirchMurnFCC = rssBirchMurnFit  
! BCC    
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      print *,"BCC aLat" 
    End If
! Make BCC configuration and calculate alat
    Call makeConfigFile(2,2.50D0)
! Read in configuration file
    Call readConfigFile(0,configID)
    Call makeNeighbourList(configID,configID)  
! Estimate BCC lattice parameter
    !Call estimateALat(configID,bccALat, bccEMin, bccVolMin, bccBM, &
    !                 coefficientsMurn, coefficientsBirchMurn, &
    !                 rssMurnFit, rssBirchMurnFit) 
    Call estimateALat(configID, &
                     coefficientsMurn, coefficientsBirchMurn, &
                     rssMurnFit, rssBirchMurnFit)     
! Murnaghan fit  
    bccEMinMurn = coefficientsMurn(1)
    bccVolMinMurn = coefficientsMurn(2)
    bccALatMurn = (bccVolMinMurn*2)**(1.0D0/3.0D0)
    bccBMMurn = coefficientsMurn(3)*160.2176568   
    bccBMPMurn = coefficientsMurn(4)
    rssMurnBCC = rssBirchMurnFit   
! Birch-Murnaghan fit    
    bccEMinBirchMurn = coefficientsBirchMurn(1)
    bccVolMinBirchMurn = coefficientsBirchMurn(2)
    bccALatBirchMurn = (bccVolMinBirchMurn*2)**(1.0D0/3.0D0)
    bccBMBirchMurn = coefficientsBirchMurn(3)*160.2176568   
    bccBMPBirchMurn = coefficientsBirchMurn(4)
    rssBirchMurnBCC = rssBirchMurnFit  
    
    If(testingFitChoice.eq.1)Then 
! Save Murn as the values to use    
      fccALat = fccALatMurn
      fccVolMin = fccVolMinMurn
      fccEMin = fccEMinMurn
      fccBM = fccBMMurn
      fccBMP = fccBMPMurn
      bccALat = bccALatMurn
      bccVolMin = bccVolMinMurn
      bccEMin = bccEMinMurn
      bccBM = bccBMMurn
      bccBMP = bccBMPMurn    
    End If
    If(testingFitChoice.eq.2)Then 
! Save BirchMurn as the values to use    
      fccALat = fccALatBirchMurn
      fccVolMin = fccVolMinBirchMurn
      fccEMin = fccEMinBirchMurn
      fccBM = fccBMBirchMurn
      fccBMP = fccBMPBirchMurn
      bccALat = bccALatBirchMurn
      bccVolMin = bccVolMinBirchMurn
      bccEMin = bccEMinBirchMurn
      bccBM = bccBMBirchMurn
      bccBMP = bccBMPBirchMurn
    End If
!--------------------------------------------
! Elastic Constants
!--------------------------------------------   
! FCC
    Call makeConfigFile(1,fccALat)     ! Make with optimised lattice parameter
    Call readConfigFile(0,configID)              ! Read in configuration file
    Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Calculate elastic constants    
    Call calcSpecificEC(configID, fccALat, fccVolMin, fccBM, fccEC, -1)
! BCC
    Call makeConfigFile(1,bccALat)     ! Make with optimised lattice parameter
    Call readConfigFile(0,configID)              ! Read in configuration file
    Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Calculate elastic constants    
    Call calcSpecificEC(configID, bccALat, bccVolMin, bccBM, bccEC, -1)
!--------------------------------------------
! Testing RSS
!--------------------------------------------
    fccCalcValues = -2.1D20
    bccCalcValues = -2.1D20
! Fill calculated arrays
    fccCalcValues(1) = fccALat * rssWeighting(4)
    fccCalcValues(2) = fccEMin * rssWeighting(5)
    fccCalcValues(3) = fccBM * rssWeighting(6)
    fccCalcValues(4) = fccEC(1) * rssWeighting(7)
    fccCalcValues(5) = fccEC(2) * rssWeighting(7)
    fccCalcValues(6) = fccEC(3) * rssWeighting(7)
!    bccCalcValues
    testingRSS = 0.0D0
    Do i=1,size(fccReferenceValues,1)
      If(fccReferenceValues(i).gt.-2.0D20.and.fccCalcValues(i).gt.-2.0D20)Then
        testingRSS = testingRSS + (fccReferenceValues(i)-fccCalcValues(i))**2
      End If
    End Do   
!--------------------------------------------
! Output
!--------------------------------------------   
    Call outputTestingSummary() 
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      Call outputTestingSummaryT() 
    End If    
  End Subroutine runTestEAM
  
  
  
  !Subroutine estimateALat(configID, aLatMin, eMin, volumeMin, bm, &
  !                       coefficientsMurn, coefficientsBirchMurn,&
  !                       rssMurnFit,rssBirchMurnFit)
  Subroutine estimateALat(configID, &
                         coefficientsMurn, coefficientsBirchMurn,&
                         rssMurnFit,rssBirchMurnFit)
! Assumes just 1 config - estimates alat
! ONLY CALCULATES ON CONFIGID 1
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: i, j, minPoint, sPoint, ePoint, configID
    Real(kind=DoubleReal) :: aLatIn, aLatTest, configEnergy, volume
    Real(kind=DoubleReal) :: rssMurnFit,rssBirchMurnFit
    Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, eMin, aLatNew, aLatMin, volumeMin
    Real(kind=DoubleReal), Dimension(1:65,1:5) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:12,1:2) :: dataPointsFit
    Real(kind=DoubleReal), Dimension(1:6,1:2) :: dataPointsFitS
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsMurn  
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBirchMurn
! Init variables    
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
        minPoint = i
      End If
    End Do
! Output points
    Call outputALatTest(dataPoints)
! Make fitting points array
    sPoint = minPoint - 2
    ePoint = minPoint + 3
    If(sPoint.lt.1)Then
      sPoint = 1
      ePoint = 6
    End If
    If(ePoint.gt.27)Then
      sPoint = 22
      ePoint = 27
    End If
! Lattice parameter
    dataPointsFitS = 0.0D0
    j = 0
    Do i=sPoint,ePoint
      j = j + 1
      dataPointsFitS(j,1) = dataPoints(i,1)
      dataPointsFitS(j,2) = dataPoints(i,2)
    End Do    
! Fit points and find minimum
    aLatMin = MinPolyFit(dataPointsFitS,3)   
! Output minimum to terminal
    !If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
    !  print *,"Min alat attempt 1: ",aLatMin
    !End If 
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
        Call loadConfigNL(configID)
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
! Volume
    dataPointsFit = 0.0D0
    j = 0
    Do i=sPoint,ePoint
      j = j + 1
      dataPointsFit(j,1) = dataPoints(i,5)
      dataPointsFit(j,2) = dataPoints(i,2)
    End Do       
! Murnaghan Fit
    coefficientsMurn = MurnFit(dataPointsFit)
    rssMurnFit = MurnRSS(dataPointsFit,coefficientsMurn)
! Birch-Murnaghan Fit
    coefficientsBirchMurn = BirchMurnFit(dataPointsFit, 0.01D0, 2000, 5, coefficientsMurn)
    rssBirchMurnFit = BirchMurnRSS(dataPointsFit,coefficientsBirchMurn)
    
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