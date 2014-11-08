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
  Use relax  
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
  Real(kind=DoubleReal) :: tempRSS  
! Public Subroutines
  Public :: runTestAnalysis
  Public :: runTestEAM
  Public :: calcEOSFit  
Contains 
!---------------------------------------------------------------------------------------------------
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
      print *,"EoS Fit RSS: ",eosFitRSS      
      print *,"-----------------------------------------------------"    
      print *,""   
      print *,"Testing-Reference RSS:    ",testingRSS
      print *,"" 
      print *,"================================================================="
    End If   
! Save aLat vs Energy Data points    
    If(mpiProcessID.eq.0)Then 
      Open(UNIT=25,FILE=Trim(outputDirectory)//"/alat-vs-energy.dat")
      write(25,"(A30)") "=============================="
      write(25,"(A30)") "ALat vs Energy                "
      write(25,"(A30)") "=============================="
      Close(25)      
    End If    
    
    If(mpiProcessID.eq.0)Then 
      Open(UNIT=25,FILE=Trim(outputDirectory)//"/alat-vs-energy.dat",&
      status="old",position="append",action="write") 
      write(25,"(A3)") " "
      write(25,"(A3)") "FCC"
      write(25,"(A3)") " "
      Close(25)      
    End If    
    Call aLatEnergy(901, 1, fccALatMurn)
    If(mpiProcessID.eq.0)Then 
      Open(UNIT=25,FILE=Trim(outputDirectory)//"/alat-vs-energy.dat",&
      status="old",position="append",action="write") 
      write(25,"(A3)") " "
      write(25,"(A3)") "BCC"
      write(25,"(A3)") " "
      Close(25)      
    End If 
    Call aLatEnergy(901, 2, bccALatMurn)
  End Subroutine runTestAnalysis   
!---------------------------------------------------------------------------------------------------
  Subroutine runTestEAM(printTestingIn)
! Tests the EAM potential
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger), optional :: printTestingIn
    Integer(kind=StandardInteger) :: printTesting
    !Real(kind=DoubleReal), Dimension(1:4) :: coefficientsMurn 
    !Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBirchMurn 
    !Real(kind=DoubleReal) :: rssMurnFit,rssBirchMurnFit
    Integer(kind=StandardInteger) :: aLatSwitch, ecSwitch, eosSwitch, vacSwitch
! Optional variables
    printTesting = 0
    If(present(printTestingIn))Then
      printTesting = printTestingIn
    End If
! Init variables
    aLatSwitch = 1   ! aLat and BM calc switch
    ecSwitch = 1     ! elastic constant switch
    eosSwitch = 1    ! Equation of State
    vacSwitch = 1
! Print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      print *,"---------------------------------------------------------------------"
      print *,"Run Test EAM Subroutine"
      print *,"---------------------------------------------------------------------"
    End If
! Set ID slot to use for temporary configs/neighbour lists
    configID = 901
!--------------------------------------------
! EoS
!--------------------------------------------
    If(aLatSwitch.eq.1.or.ecSwitch.eq.1)Then    
! EoS Fit FCC  
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
        print *,"FCC aLat" 
      End If   
      Call eosFit(configID,1)
! EoS Fit BCC  
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
        print *,"BCC aLat" 
      End If   
      Call eosFit(configID,2)
    End If
!--------------------------------------------
! Elastic Constants
!--------------------------------------------   
! FCC
    If(ecSwitch.eq.1)Then
      Call makeConfigFile(1,fccALat)     ! Make with optimised lattice parameter
      Call readConfigFile(0,configID)              ! Read in configuration file
      Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Calculate elastic constants    
      Call calcSpecificEC(configID, fccALat, fccVolMin, fccBM, fccEC, -1)
! BCC
      Call makeConfigFile(2,bccALat)     ! Make with optimised lattice parameter
      Call readConfigFile(0,configID)              ! Read in configuration file
      Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Calculate elastic constants    
      Call calcSpecificEC(configID, bccALat, bccVolMin, bccBM, bccEC, -1)
    End If
    
    
!--------------------------------------------
! Testing RSS
!--------------------------------------------
    fccCalcValues = -2.1D20
    bccCalcValues = -2.1D20
! Fill calculated arrays
    fccCalcValues(1) = fccALat
    fccCalcValues(2) = fccEMin
    fccCalcValues(3) = fccBM
    fccCalcValues(4) = fccEC(1)
    fccCalcValues(5) = fccEC(2)
    fccCalcValues(6) = fccEC(3)
! Calculate RSS
    testingALatRSS = (fccReferenceValues(1)-fccALat)**2
    testingEMinRSS = (fccReferenceValues(2)-fccEMin)**2    
    testingBMRSS = (fccReferenceValues(3)-fccBM)**2    
    testingECRSS = (fccReferenceValues(4)-fccEC(1))**2+&
                   (fccReferenceValues(5)-fccEC(2))**2+&
                   (fccReferenceValues(6)-fccEC(3))**2
! Total testing RSS
    testingRSS = testingALatRSS+testingEMinRSS+testingBMRSS+testingECRSS+eosFitRSS
!--------------------------------------------
! Output
!--------------------------------------------   
    !Call outputTestingSummary() 
    !If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      !Call outputTestingSummaryT() 
    !End If    
  End Subroutine runTestEAM
  
!---------------------------------------------------------------------------------------------------
  Subroutine eosFit(configID, structureType)
! Assumes just configid - fits birch murnagham EoS
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: configID, structureType   
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:12,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: dataPointsS
    Real(kind=DoubleReal) :: eMin, aLatMin, aLatTest, aLatStart, aLatNL
    Real(kind=DoubleReal) :: configEnergy, pairEnergy, embeddingEnergy, volume
    Real(kind=DoubleReal) :: configEnergyBM, tempRSSAdd
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBirchMurn
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTest
! Init variables    
    dataPoints = -2.1D20
    eMin = 2.1D20
    aLatMin = 0.0d0    
    If(structureType.eq.1)Then  ! FCC
      aLatStart = 2.5D0
    End If     
    If(structureType.eq.2)Then  ! BCC
      aLatStart = 2.5D0
    End If  
!--------------------------------------------
! Estimate A Lat (between 2.5 and 8.5)
!--------------------------------------------
    aLatNL = aLatStart
    Do i=1,12
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then   ! Split between multiple processes
        Call makeConfigFile(structureType,aLatNL)      ! Make FCC(1)/BCC(2) configuration and calculate alat
        Call readConfigFile(0,configID)           ! Read in configuration file into config slot 901
        Call makeNeighbourList(configID,configID,0)           ! Make new neighbour list   
! Set testing alat
        aLatTest = aLatStart + ((i-1) * 0.5D0)
        Call changeALat(configID, aLatNL, aLatTest)
! Calculate volume      
        Call calcVol(configID, aLatTest, volume)
! Calculate energy
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)
! Store data points      
        dataPoints(i,1) = aLatTest
        dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
! Store aLatNL        
        aLatNL = aLatNL + mpiProcessCount * 0.5D0
      End If
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPoints)    
    Call M_distDouble2D(dataPoints)
    Do i=1,12
      If(dataPoints(i,2).lt.eMin)Then
        aLatMin = dataPoints(i,1) 
        eMin = dataPoints(i,2)
      End If
    End Do    
!--------------------------------------------
! Estimate A Lat (between aLatMin-0.5D0 and aLatMin+0.5D0)
!--------------------------------------------   
    aLatStart = aLatMin-0.5D0
    Call makeConfigFile(structureType,aLatStart)      ! Make FCC(1)/BCC(2) configuration and calculate alat
    Call readConfigFile(0,configID)           ! Read in configuration file into config slot 901
    Call makeNeighbourList(configID,configID,0)           ! Make new neighbour list  
    Call saveConfigNL(configID)    
    Do i=1,4
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then   ! Split between multiple processes
 ! Load original neighbour list    
        Call loadConfigNL(configID)   
 ! Set testing alat
        aLatTest = aLatStart + ((i-1) / 4.0D0)   
        Call changeALat(configID, aLatStart, aLatTest)
 ! Calculate volume      
        Call calcVol(configID, aLatTest, volume)
! Calculate energy
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)  
! Store data points      
        dataPointsS(i,1) = volume/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPointsS(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
      End If
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPointsS)    
    Call M_distDouble2D(dataPointsS)
! Fit data
    coefficientsBirchMurn = BirchMurnFit(dataPointsS, 0.01D0, 100, 2)
! Predicted aLat min
    If(structureType.eq.1)Then  ! FCC
      aLatMin = (coefficientsBirchMurn(2)*4.0)**(1.0D0/3.0D0)
    End If     
    If(structureType.eq.2)Then  ! BCC
      aLatMin = (coefficientsBirchMurn(2)*2.0)**(1.0D0/3.0D0)
    End If  
!--------------------------------------------
! Estimate A Lat (near to aLatMin)
!--------------------------------------------   
    aLatStart = aLatMin-(5*0.05D0)
    Call makeConfigFile(structureType,aLatStart)      ! Make FCC(1)/BCC(2) configuration and calculate alat
    Call readConfigFile(0,configID)           ! Read in configuration file into config slot 901
    Call makeNeighbourList(configID,configID,0)           ! Make new neighbour list  
    Call saveConfigNL(configID)    
    Do i=1,12
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then   ! Split between multiple processes
 ! Load original neighbour list    
        Call loadConfigNL(configID)   
 ! Set testing alat
        aLatTest = aLatStart + ((i-1) * 0.05D0)   
        Call changeALat(configID, aLatStart, aLatTest)
 ! Calculate volume      
        Call calcVol(configID, aLatTest, volume)
! Calculate energy
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)  
! Store data points      
        dataPoints(i,1) = volume/(1.0D0*configurationCoordsKeyG(configID,2))
        dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
      End If
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPoints)    
    Call M_distDouble2D(dataPoints)
! Fit data
    coefficientsBirchMurn = BirchMurnFit(dataPoints, 0.01D0, 500, 3)
! Predicted aLat min
    If(structureType.eq.1)Then  ! FCC
      aLatMin = (coefficientsBirchMurn(2)*4.0)**(1.0D0/3.0D0)
    End If     
    If(structureType.eq.2)Then  ! BCC
      aLatMin = (coefficientsBirchMurn(2)*2.0)**(1.0D0/3.0D0)
    End If  
!--------------------------------------------
! Save EoS Parameters
!-------------------------------------------- 
    If(structureType.eq.1)Then  ! FCC
      fccALat = aLatMin
      fccVolMin = coefficientsBirchMurn(2)
      fccEMin = coefficientsBirchMurn(1)
      fccBM = coefficientsBirchMurn(3)*160.2176568 
      fccBMP = coefficientsBirchMurn(4)
    End If
    If(structureType.eq.2)Then  ! BCC
      bccALat = aLatMin
      bccVolMin = coefficientsBirchMurn(2)
      bccEMin = coefficientsBirchMurn(1)
      bccBM = coefficientsBirchMurn(3)*160.2176568 
      bccBMP = coefficientsBirchMurn(4)
    End If
!--------------------------------------------
! EoS Test vs Experimental RSS
!-------------------------------------------- 
    eosFitRSS = 0.0D0
    If(structureType.eq.1)Then  ! FCC
      coefficientsTest(1) = fccReferenceValues(2)                 ! E0
      coefficientsTest(2) = (fccReferenceValues(1)**3)/4.0D0      ! V0
      coefficientsTest(3) = fccReferenceValues(3)/160.2176568     ! B0
      coefficientsTest(4) = 3.0D0                                 ! B'0
      coefficientsTest = BirchMurnFitBP(dataPoints, coefficientsTest)     ! No reference B'0 so adjust to match calculated
! Loop through points
      Do i=1,12
        configEnergyBM = BirchMurnCalc(dataPoints(i,1),coefficientsTest)
        tempRSSAdd = (configEnergyBM-dataPoints(i,2))**2
        eosFitRSS = eosFitRSS + tempRSSAdd
      End Do   
    End If  
  End Subroutine eosFit
  
  
  
  
  
  
  
!---------------------------------------------------------------------------------------------------
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
  Subroutine aLatEnergy(configID, structureType, aLatOpt)
! Makes array of alat-energy data points
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: configID, i, j       ! probably 901
    Integer(kind=StandardInteger) :: structureType  ! 1 FCC, 2 BCC
    Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
    Real(kind=DoubleReal) :: aLatTest, aLatOpt, volume, energyM, energyBM
    Real(kind=DoubleReal) :: configEnergy, pairEnergy, embeddingEnergy
    Real(kind=DoubleReal), Dimension(1:41,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:41,1:4) :: dataPointsPrint
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsMurn  
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBirchMurn
! Init variables
    If(structureType.eq.1)Then     ! FCC
      xCopy = 4
      yCopy = 4
      zCopy = 4
    End If
    If(structureType.eq.2)Then     ! BCC
      xCopy = 5
      yCopy = 5
      zCopy = 5
    End If
! Open file       
    If(mpiProcessID.eq.0)Then 
      Open(UNIT=25,FILE=Trim(outputDirectory)//"/alat-vs-energy.dat",&
      status="old",position="append",action="write") 
    End If      
! Make config + neighbour list    
    Call makeConfigFile(structureType,2.0D0)     ! Make with optimised lattice parameter
    Call readConfigFile(0,configID)                ! Read in configuration file
    Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Save neighbour list    
    Call saveConfigNL(configID)
! Lattice/Energy calculations (2.0D0 Ang to 10.0D0 Ang)
    dataPointsPrint = 0.0D0
    Do i=1,41
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then   ! Split between multiple processes
        aLatTest = 2.0D0 +((i-1)*0.2D0)
        Call isotropicDistortion(configID, 0.0D0, aLatTest/2.0D0)
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)      
        configEnergy = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        pairEnergy = pairEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        embeddingEnergy = embeddingEnergy/(1.0D0*configurationCoordsKeyG(configID,2))        
        dataPointsPrint(i,1) = aLatTest
        dataPointsPrint(i,2) = configEnergy
        dataPointsPrint(i,3) = pairEnergy
        dataPointsPrint(i,4) = embeddingEnergy
! Reload list        
        Call loadConfigNL(configID)
      End If
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPointsPrint)    
    Call M_distDouble2D(dataPointsPrint)
! output data
    If(mpiProcessID.eq.0)Then
      write(25,"(A8,A16,A16,A16,A16)") &
      "        ","aLat/Ang        ","Energy/eV       ",&        
      "Pair Energy/eV  ","Embe Energy/eV  "
      Do j=1,41
        write(25,"(I8,F16.8,F16.8,F16.8,F16.8)") j,&
        dataPointsPrint(j,1),dataPointsPrint(j,2),&
        dataPointsPrint(j,3),dataPointsPrint(j,4)
      End Do
    End If        
    write(25,"(A1)") " "
! Make config + neighbour list    
    Call makeConfigFile(structureType,(aLatOpt-20*0.005D0))     ! Make with optimised lattice parameter
    Call readConfigFile(0,configID)                ! Read in configuration file
    Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Save neighbour list    
    Call saveConfigNL(configID)    
! Lattice/Energy calculations (2.0D0 Ang to 10.0D0 Ang)
    dataPoints = 0.0D0
    dataPointsPrint = 0.0D0
    j = 0
    Do i=-20,20
      j = j + 1
      If(mpiProcessID.eq.mod(j-1,mpiProcessCount))Then   ! Split between multiple processes
        aLatTest = aLatOpt+(i*0.005D0) 
        Call isotropicDistortion(configID, 0.0D0, aLatTest/(aLatOpt-20*0.005D0))        
        Call calcEnergy(configID, configEnergy, 0, pairEnergy, embeddingEnergy)      
        configEnergy = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        pairEnergy = pairEnergy/(1.0D0*configurationCoordsKeyG(configID,2))
        embeddingEnergy = embeddingEnergy/(1.0D0*configurationCoordsKeyG(configID,2))  
        If(structureType.eq.1)Then
          dataPoints(j,1) = (aLatTest**3)/4.0D0
        End If
        If(structureType.eq.2)Then
          dataPoints(j,1) = (aLatTest**3)/2.0D0
        End If
        dataPoints(j,2) = configEnergy
        dataPointsPrint(j,1) = aLatTest
        dataPointsPrint(j,2) = configEnergy
        dataPointsPrint(j,3) = pairEnergy
        dataPointsPrint(j,4) = embeddingEnergy
! Reload list        
        Call loadConfigNL(configID)
      End If  
    End Do
! Collect and distribute data array   
    Call M_collDouble2D(dataPoints) 
    Call M_collDouble2D(dataPointsPrint)    
    Call M_distDouble2D(dataPoints)  
    Call M_distDouble2D(dataPointsPrint)
! Murnaghan Fit
    coefficientsMurn = MurnFit(dataPoints)
! Birch-Murnaghan Fit
    coefficientsBirchMurn = BirchMurnFit(dataPoints, 0.01D0, 2000, 5, coefficientsMurn)
! output data
    If(mpiProcessID.eq.0)Then
      write(25,"(A8,A16,A16,A16,A16,A16,A16)") &
      "        ","aLat/Ang        ","Energy/eV       ",&        
      "Murn Energy/eV  ","BM Energy/eV    ",&
      "Pair Energy/eV  ","Embe Energy/eV  "
      Do j=1,41
        volume = dataPoints(j,1)
        energyM = MurnCalc(volume, coefficientsMurn)
        energyBM = BirchMurnCalc(volume, coefficientsMurn)
        write(25,"(I8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)") j,&
        dataPointsPrint(j,1),dataPointsPrint(j,2),&
        energyM,energyBM,&
        dataPointsPrint(j,3),dataPointsPrint(j,4)
      End Do
    End If
    write(25,"(A1)") " "
! Close file       
    If(mpiProcessID.eq.0)Then 
      Close(25)
    End If            
  End Subroutine aLatEnergy
!---------------------------------------------------------------------------------------------------
  Subroutine calcEOSFit(configID, structureType, aLatOpt)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID, structureType, i
    Real(kind=DoubleReal) :: aLatTest, aLatLower, aLatOpt, eosFitRSSValues
    Real(kind=DoubleReal) :: configEnergy, configEnergyBM, tempRSSAdd
    Real(kind=DoubleReal), Dimension(1:16,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients ! E0, V0, B0, B'0
! Init values    
    eosFitRSSValues = 0.0D0
! Starting aLat
    aLatLower = aLatOpt - 0.5D0  
! Make config + neighbour list    
    Call makeConfigFile(structureType,aLatLower)     ! Make with optimised lattice parameter
    Call readConfigFile(0,configID)                ! Read in configuration file
    Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Save neighbour list    
    Call saveConfigNL(configID)    
! Calculate points from EAM potential
    dataPoints = 0.0D0
    Do i=1,16
      If(mpiProcessID.eq.mod(i-1,mpiProcessCount))Then       
        aLatTest = aLatLower + (i-1)*(0.5D0/7.0D0)
        Call isotropicDistortion(configID, 0.0D0, aLatTest/aLatLower)
        Call calcEnergy(configID, configEnergy, 0) 
        If(structureType.eq.1)Then  !FCC
          dataPoints(i,1) = (aLatTest**3/4.0D0)
        End If  
        If(structureType.eq.2)Then  !BCC
          dataPoints(i,1) = (aLatTest**3/2.0D0)
        End If  
        dataPoints(i,2) = configEnergy/(1.0D0*configurationCoordsKeyG(configID,2)) 
  ! Reload list        
        Call loadConfigNL(configID)
      End If
    End Do   
! Collect and distribute data array   
    Call M_collDouble2D(dataPoints)   
    Call M_distDouble2D(dataPoints)    
! Compare these points to 
    If(structureType.eq.1)Then  !FCC
      coefficients(1) = fccReferenceValues(2)                 ! E0
      coefficients(2) = (fccReferenceValues(1)**3)/4.0D0      ! V0
      coefficients(3) = fccReferenceValues(3)/160.2176568     ! B0
      coefficients(4) = 3.0D0                                 ! B'0
    End If  
    If(structureType.eq.2)Then  !BCC   
          
    End If  
! Fit B'0 to points
    coefficients = BirchMurnFitBP(dataPoints, coefficients)   
! Calculate rss    
    Open(unit=969,file=trim(outputDirectory)//"/"//"EoSFitting.dat",&
    status="old",position="append",action="write")  
    Write(969,"(A13,F8.4)") "EoS Fitting  ",ProgramTime()  
    Do i=1,16 
      configEnergyBM = BirchMurnCalc(dataPoints(i,1),coefficients)
      tempRSSAdd = (configEnergyBM-dataPoints(i,2))**2
      eosFitRSSValues = eosFitRSSValues + tempRSSAdd
      If(i.lt.10)Then
        Write(969,"(I1,A2,F12.6,A1,F12.6,A1,F12.6,A1,E12.6)") &
        i," ",dataPoints(i,1)," ",dataPoints(i,2)," ",configEnergyBM," ",tempRSSAdd
      Else  
        Write(969,"(I2,A1,F12.6,A1,F12.6,A1,F12.6,A1,E12.6)") &
        i," ",dataPoints(i,1)," ",dataPoints(i,2)," ",configEnergyBM," ",tempRSSAdd
      End If  
    End Do  
    Write(969,"(A13,E12.6)") "EoS Fit RSS: ",eosFitRSSValues
    Write(969,"(A1)") " "
    Write(969,"(A1)") " "
    Close(969)
    eosFitRSS = eosFitRSS + eosFitRSSValues
  End Subroutine calcEOSFit
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
!---------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------
  Subroutine vacancyEnergy()
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables   
    !Integer(kind=StandardInteger) :: configID  
  
  
  
  End Subroutine vacancyEnergy

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !---------------------------------------------------------------------------------------------------
  Subroutine runTestEAMOld(printTestingIn)
! Tests the EAM potential
    Implicit None   ! Force declaration of all variables
! Private variables       
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger), optional :: printTestingIn
    Integer(kind=StandardInteger) :: printTesting
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsMurn 
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBirchMurn 
    Real(kind=DoubleReal) :: rssMurnFit,rssBirchMurnFit
    Integer(kind=StandardInteger) :: aLatSwitch, ecSwitch, eosSwitch, vacSwitch
! Optional variables
    printTesting = 0
    If(present(printTestingIn))Then
      printTesting = printTestingIn
    End If
! Init variables
    aLatSwitch = 1   ! aLat and BM calc switch
    ecSwitch = 1     ! elastic constant switch
    eosSwitch = 1    ! Equation of State
    vacSwitch = 1
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
    If(aLatSwitch.eq.1.or.ecSwitch.eq.1)Then    
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
        print *,"FCC aLat" 
      End If
! FCC
      Call makeConfigFile(1,2.50D0)      ! Make FCC configuration and calculate alat
      Call readConfigFile(0,configID)           ! Read in configuration file into config slot 901
      Call makeNeighbourList(configID,configID,printTesting)           ! Make new neighbour list
! Estimate FCC lattice parameter
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
    End If
!--------------------------------------------
! Elastic Constants
!--------------------------------------------   
! FCC
    If(ecSwitch.eq.1)Then
      Call makeConfigFile(1,fccALat)     ! Make with optimised lattice parameter
      Call readConfigFile(0,configID)              ! Read in configuration file
      Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Calculate elastic constants    
      Call calcSpecificEC(configID, fccALat, fccVolMin, fccBM, fccEC, -1)
! BCC
      Call makeConfigFile(2,bccALat)     ! Make with optimised lattice parameter
      Call readConfigFile(0,configID)              ! Read in configuration file
      Call makeNeighbourList(configID,configID)           ! Make new neighbour list    
! Calculate elastic constants    
      Call calcSpecificEC(configID, bccALat, bccVolMin, bccBM, bccEC, -1)
    End If
!--------------------------------------------
! EoS RSS
!--------------------------------------------    
    If(eosSwitch.eq.1)Then
      eosFitRSS = 0.0D0
      If(eosFitRSSOption.eq.1)Then
        Call calcEOSFit(configID, 1, fccReferenceValues(1))    
      End If
    End If
!--------------------------------------------
! EoS RSS
!--------------------------------------------     
    If(vacSwitch.eq.1)Then
      Call vacancyEnergy()
    End If
!--------------------------------------------
! Testing RSS
!--------------------------------------------
    fccCalcValues = -2.1D20
    bccCalcValues = -2.1D20
! Fill calculated arrays
    fccCalcValues(1) = fccALat
    fccCalcValues(2) = fccEMin
    fccCalcValues(3) = fccBM
    fccCalcValues(4) = fccEC(1)
    fccCalcValues(5) = fccEC(2)
    fccCalcValues(6) = fccEC(3)
! Calculate RSS
    testingALatRSS = (fccReferenceValues(1)-fccALat)**2
    testingEMinRSS = (fccReferenceValues(2)-fccEMin)**2    
    testingBMRSS = (fccReferenceValues(3)-fccBM)**2    
    testingECRSS = (fccReferenceValues(4)-fccEC(1))**2+&
                   (fccReferenceValues(5)-fccEC(2))**2+&
                   (fccReferenceValues(6)-fccEC(3))**2
! Total testing RSS
    testingRSS = testingALatRSS+testingEMinRSS+testingBMRSS+testingECRSS+eosFitRSS
!--------------------------------------------
! Output
!--------------------------------------------   
    Call outputTestingSummary() 
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      Call outputTestingSummaryT() 
    End If    
  End Subroutine runTestEAMOld
  
  
  
End Module testEAM