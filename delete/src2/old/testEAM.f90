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
  Contains
! ---------------------------------------------------------------------------------------------------
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
! ---------------------------------------------------------------------------------------------------
  Subroutine runTestEAM(printTestingIn)
! Tests the EAM potential
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger), optional :: printTestingIn
    Integer(kind=StandardInteger) :: printTesting
    Integer(kind=StandardInteger) :: ecSwitch, eosSwitch, vacSwitch
    Real(kind=DoubleReal) :: testStartTime, testEndTime
    testStartTime = ProgramTime()
! Optional variables
    printTesting = 0
    If(present(printTestingIn))Then
      printTesting = printTestingIn
    End If
! Init variables
    eosSwitch = 1    ! Equation of State
    ecSwitch = 1     ! elastic constant switch
    vacSwitch = 1
    testConfigRSS = 0.0D0
! Print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
      print *,"---------------------------------------------------------------------"
      print *,"Run Test EAM Subroutine"
      print *,"---------------------------------------------------------------------"
    End If
! Set ID slot to use for temporary configs/neighbour lists
    configID = 901
! --------------------------------------------
! EoS
! --------------------------------------------
    If(eosSwitch.eq.1)Then
! EoS Fit FCC
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
        print *,"FCC aLat"
      End If
      Call eosFit(configID,1)                        ! FCC EoS Fit
! Store RSS values
      If(fccReferenceValues(1).gt.-2.0D0)Then
        testConfigRSS(1) = (fccReferenceValues(1)-fccALat)**2   ! FCC ALat
      End If
      If(fccReferenceValues(2).gt.-2.0D0)Then
        testConfigRSS(2) = (fccReferenceValues(2)-fccEMin)**2   ! FCC EMin
      End If
      If(fccReferenceValues(3).gt.-2.0D0)Then
        testConfigRSS(3) = (fccReferenceValues(3)-fccBM)**2     ! FCC BM
      End If
      testConfigRSS(4) = eosFitRSS
! EoS Fit BCC
      If(mpiProcessID.eq.0.and.printToTerminal.eq.1.and.printTesting.eq.1)Then
        print *,"BCC aLat"
      End If
      Call eosFit(configID,2)                        ! BCC EoS Fit
! Store RSS values
      If(bccReferenceValues(1).gt.-2.0D0)Then
        testConfigRSS(8) = (bccReferenceValues(1)-bccALat)**2   ! BCC ALat
      End If
      If(bccReferenceValues(2).gt.-2.0D0)Then
        testConfigRSS(9) = (bccReferenceValues(2)-bccEMin)**2   ! BCC EMin
      End If
      If(bccReferenceValues(3).gt.-2.0D0)Then
        testConfigRSS(10) = (bccReferenceValues(3)-bccBM)**2    ! BCC BM
      End If
      testConfigRSS(11) = eosFitRSS
    End If
    print *,testConfigRSS(1),testConfigRSS(2),testConfigRSS(3)
! --------------------------------------------
! Elastic Constants
! --------------------------------------------
! FCC
    If(ecSwitch.eq.1)Then
      Call makeConfigFile(1,fccALat)     ! Make with optimised lattice parameter
      Call readConfigFile(0,configID)              ! Read in configuration file
      Call makeNeighbourList(configID,configID)           ! Make new neighbour list
! Calculate elastic constants
      Call calcSpecificEC(configID, fccALat, fccVolMin, fccBM, fccEC, -1)
! Store RSS values
      If(fccReferenceValues(4).gt.-2.0D0)Then
        testConfigRSS(5) = (fccReferenceValues(4)-fccEC(1))**2   ! FCC ALat
      End If
      If(fccReferenceValues(5).gt.-2.0D0)Then
        testConfigRSS(6) = (fccReferenceValues(5)-fccEC(2))**2   ! FCC ALat
      End If
      If(fccReferenceValues(6).gt.-2.0D0)Then
        testConfigRSS(7) = (fccReferenceValues(6)-fccEC(3))**2   ! FCC ALat
      End If
! BCC
      Call makeConfigFile(2,bccALat)     ! Make with optimised lattice parameter
      Call readConfigFile(0,configID)              ! Read in configuration file
      Call makeNeighbourList(configID,configID)           ! Make new neighbour list
! Calculate elastic constants
      Call calcSpecificEC(configID, bccALat, bccVolMin, bccBM, bccEC, -1)
    End If
! --------------------------------------------
! Testing RSS
! --------------------------------------------
! Total testing RSS
    testingRSS = 0.0D0
    testingRSS = testingALatRSS+testingEMinRSS+testingBMRSS+testingECRSS+eosFitRSS
! --------------------------------------------
! Output
! --------------------------------------------
    testEndTime = ProgramTime()
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A25,F8.6,A6,E16.8)") &
      "Test config calc.  Time: ",(testEndTime-testStartTime)," RSS: ",testingRSS
    End If
  End Subroutine runTestEAM
! ---------------------------------------------------------------------------------------------------
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
! --------------------------------------------
! Estimate A Lat (between 2.5 and 8.5)
! --------------------------------------------
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
! --------------------------------------------
! Estimate A Lat (between aLatMin-0.5D0 and aLatMin+0.5D0)
! --------------------------------------------
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
! --------------------------------------------
! Estimate A Lat (near to aLatMin)
! --------------------------------------------
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
! --------------------------------------------
! Save EoS Parameters
! --------------------------------------------
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
! --------------------------------------------
! EoS Test vs Experimental RSS
! --------------------------------------------
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
! ---------------------------------------------------------------------------------------------------
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
! ---------------------------------------------------------------------------------------------------
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
! ---------------------------------------------------------------------------------------------------
  Subroutine calcVol(configID, aLat, volume)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
    Real(kind=DoubleReal) :: aLat, volume
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: superCellVectors,configUnitVector
! Integer(kind=StandardInteger) :: i, keyS, keyE
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
! ---------------------------------------------------------------------------------------------------
  Subroutine vacancyEnergy()
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
! Integer(kind=StandardInteger) :: configID
  End Subroutine vacancyEnergy

End Module testEAM
