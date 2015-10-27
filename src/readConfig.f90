Module readConfig

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
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use initialise
  Use loadData
  Use globals
  Use output
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: readConfigFile
  Public :: readBpConfigFile

  Contains
  Subroutine readConfigFile()
    Implicit None   ! Force declaration of all variables
! Print out
    If(TerminalPrint())Then
      print *,"READ CONFIG: ",trim(configFilePath)
    End If
! Private variables
    Call cpu_time(timeStart)
! Load config file into memory
    Call loadFile()
! Read in any listed DFT files
    Call readDFTFiles()
    Call resetConfigVars()
    Call processFile()
    Call expandCoordinates()
    Call outputConfigPoints()
! Synch MPI processes
    Call M_synchProcesses()
    Call cpu_time(timeEnd)
    Call timeAcc(configLoadTime,timeStart,timeEnd)
! Output
    If(TerminalPrint())Then
      print *,"Atom configurations loaded: ",(timeEnd-timeStart),"s"
    End If
  End Subroutine readConfigFile

  Subroutine loadFile()
! Creates a temporary file to contain the config potential
! Converts alpha to upper, and strips out comment lines etc
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: ios, i, j
    Character(len=255) :: fileRow, fileRowCaps
! Output to Terminal
    If(TerminalPrint())Then
      Print *,"Loading user config file ",trim(configFilePath)
    End If
! Read config file into memory
    j = 0
    Open(UNIT=1,FILE=trim(configFilePath))
    Do i=1,maxFileRows
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
      fileRow = Trim(Adjustl(fileRow))
      fileRow = RemoveComments(fileRow)
      fileRow = RemoveQuotes(fileRow)
      If(fileRow(1:1).eq."!".or.fileRow(1:1).eq." ")Then
! skip
      Else
        j = j + 1
        fileRow = Trim(Adjustl(fileRow))
        fileRowCaps = StrToUpper(fileRow)
        If(fileRowCaps(1:5).eq."#PATH")Then
          configInputData(j) = Trim("#PATH "//fileRow(7:255))
        Else
          configInputData(j) = fileRowCaps
        End If
      End If
    End Do
    Close(1)
  End Subroutine loadFile

  Subroutine readDFTFiles()
! Creates a temporary file to contain the config potential
! Converts alpha to upper, and strips out comment lines etc
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i,j,m,n,writeFile,tempLastRow,readFlag
    Integer(kind=StandardInteger) :: replCount
    Real(kind=DoubleReal) :: radiusCutoff,confWeight
    Character(len=255) :: fileRow, fileRowB
    Character(len=255) :: dftFilePath
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD
    Character(len=32) :: dftType
! Init variables
    bufferA = BlankString(bufferA)
    bufferB = BlankString(bufferB)
    bufferC = BlankString(bufferC)
    bufferD = BlankString(bufferD)
    confWeight = 1.0D0
! Loop through file rows
    m = 0
    n = 0
    writeFile = 0
    Do i=1,maxFileRows
      fileRow = configInputData(i)     !read line
! Break out
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
! Check if normal or dft config
      If(fileRow(1:7).eq."#NEWDFT")Then
        writeFile = 2
      ElseIf(fileRow(1:4).eq."#NEW")Then
        writeFile = 1
      End If
! Save to memory
      If(writeFile.eq.1)Then
        m = m + 1
        configInputDataTemp(m) = trim(fileRow)
      End If
      If(writeFile.eq.2)Then
        n = n + 1
        configInputDataDFT(n) = trim(fileRow)
      End If
! reset writeFile flag at end
      If(fileRow(1:4).eq."#END")Then
        writeFile = 0
      End If
    End Do
    tempLastRow = m
    Do i=1,maxFileRows
      fileRow = configInputDataTemp(i)     !read line
! Break out
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
    End Do
    Do i=1,maxFileRows
      fileRow = configInputDataDFT(i)     !read line
! Break out
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
! print *,trim(fileRow)
    End Do
! Clear dft temp array
    configInputDataDFTTemp = BlankStringArray(configInputDataDFTTemp)
! Loop through DFT files
    readFlag = 0
    Do i=1,maxFileRows
      fileRow = configInputDataDFT(i)     !read line
! Break out
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
      If(readFlag.eq.0)Then
        If(fileRow(1:7).eq."#NEWDFT")Then
          readFlag = 1
          configLabelReplace = BlankStringArray(configLabelReplace)
        End If
      Else
        If(fileRow(1:3).eq."#RC")Then
          Read(fileRow,*) bufferA, bufferB, bufferC
          Read(bufferB,*) radiusCutoff
        End If
        If(fileRow(1:5).eq."#TYPE")Then
          Read(fileRow,*) bufferA, bufferB
          dftType = Trim(bufferB)
        End If
        If(fileRow(1:5).eq."#PATH")Then
          Read(fileRow,*) bufferA, bufferB
          dftFilePath = Trim(bufferB)
        End If
        If(fileRow(1:5).eq."#REPL")Then
          Read(fileRow,*) bufferA, bufferB, bufferC
          replCount = 0
          Do j=1,maxFileRows
            replCount = replCount + 1
            fileRowB = configLabelReplace(j,1)
            If(fileRowB(1:1).eq." ")Then
              EXIT
            End If
          End Do
          configLabelReplace(replCount,1) = trim(StrToUpper(bufferB))
          configLabelReplace(replCount,2) = trim(StrToUpper(bufferC))
        End If
        If(fileRow(1:7).eq."#ENDDFT")Then
          readFlag = 0
          Do j=1,maxFileRows
            fileRowB = configLabelReplace(j,1)
            If(fileRowB(1:1).eq." ")Then
              EXIT
            End If
          End Do
! Read in file
          Call readPWSCFFile(dftFilePath, radiusCutoff, confWeight)
        End If
      End If
    End Do
! Merge arrays
    Do i=1,maxFileRows
      fileRow = configInputDataDFTTemp(i)     !read line
! Break out
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
      configInputDataTemp(tempLastRow+i) = trim(fileRow)
    End Do
    configInputData = configInputDataTemp
  End Subroutine readDFTFiles

  Subroutine readPWSCFFile(dftFilePath, dftInRadiusCutoff, confWeight)
! Read in configuration file
    Implicit None  ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: ios, i, j, k, m
    Character(len=255) :: dftFilePath
    Character(len=255) :: fileRow, fileLineBuffer
    Character(len=16) :: labelTemp
    Character(len=255) :: bufferA, bufferB, bufferC
    Integer(kind=StandardInteger) :: readType, lastScf, numberOfAtoms
    Real(kind=DoubleReal) :: aLat, dftInRadiusCutoff
    Real(kind=DoubleReal) :: configTotalEnergy, configEnergyPerAtom
    Real(kind=DoubleReal) :: dftInEqVol, confWeight
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: stress
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: crystalAxes    !Unit vector
    Character(len=8), Dimension(1:4096) :: atomType
    Character(len=8) :: atomTypeTemp
    Real(kind=DoubleReal), Dimension(1:4096,1:3) :: atomCoords
    Real(kind=DoubleReal), Dimension(1:4096,1:3) :: atomForcess
    Real(kind=DoubleReal) :: energyOffset
! Init variables
    readType = 0    ! SCF calculation
    lastScf = 1
    numberOfAtoms = 0
    aLat = 0.0D0
    configTotalEnergy = 0.0D0
    configEnergyPerAtom = 0.0D0
    stress = 0.0D0
    crystalAxes = 0.0D0
    atomCoords = 0.0D0
    atomForcess = 0.0D0
! Check type of DFT calculation
    Open(UNIT=101,FILE=Trim(dftFilePath))
    Do i=1,maxFileRows
! Read in line
      Read(101,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If(ios /= 0)Then
        EXIT
      End If
      If(fileRow(1:54).eq."     A final scf calculation at the relaxed structure.")Then
        readType = 1 ! vc-relax calculation
        lastScf = 0
      End If
    End Do
    Close(101)
! Read in data
    Open(UNIT=101,FILE=Trim(dftFilePath))
    Do i=1,maxFileRows
! Read in line
      Read(101,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If(ios /= 0)Then
        EXIT
      End If
      If(fileRow(1:23).eq."Begin final coordinates")Then
        lastScf = 1
      End If
      If(fileRow(1:33).eq."     number of atoms/cell      = ")Then
        fileLineBuffer = fileRow(34:100)
        read(fileLineBuffer,*) numberOfAtoms
      End If
      If(lastScf.eq.1)Then
! Total energy
        If(fileRow(1:17).eq."!    total energy")Then
          fileLineBuffer = fileRow(33:100)
          read(fileLineBuffer,*) bufferA, bufferB
          read(bufferA,*) configTotalEnergy
          configTotalEnergy = UnitConvert(configTotalEnergy,"RY","EV")
!  (dftInOptEnergy-dftInCohEnergy)
        End If
! Stresses
        If(fileRow(1:38).eq."          total   stress  (Ry/bohr**3)")Then
          Read(101,"(A255)",IOSTAT=ios) fileRow
          read(fileRow,*) bufferA, bufferB, bufferC
          read(bufferA,*) stress(1,1)
          read(bufferB,*) stress(1,2)
          read(bufferC,*) stress(1,3)
          Read(101,"(A255)",IOSTAT=ios) fileRow
          read(fileRow,*) bufferA, bufferB, bufferC
          read(bufferA,*) stress(2,1)
          read(bufferB,*) stress(2,2)
          read(bufferC,*) stress(2,3)
          Read(101,"(A255)",IOSTAT=ios) fileRow
          read(fileRow,*) bufferA, bufferB, bufferC
          read(bufferA,*) stress(3,1)
          read(bufferB,*) stress(3,2)
          read(bufferC,*) stress(3,3)
          Do j=1,3
            Do k=1,3
              stress(j,k) = UnitConvert(stress(j,k),"RYBOH3","GPA")
            End Do
          End Do
        End If
! Lattice Parameter
        If(fileRow(1:32).eq."     lattice parameter (alat)  =")Then
          fileLineBuffer = fileRow(33:100)
          read(fileLineBuffer,*) bufferA, bufferB
          read(bufferA,*) aLat
          aLat = UnitConvert(aLat,"Bohr","A")
        End If
! Crystal Axes/Unit Vector - scf
        If(readType.eq.0.and.fileRow(1:18).eq."     crystal axes:")Then
          Read(101,"(A255)",IOSTAT=ios) fileRow
          fileLineBuffer = fileRow(24:100)
          read(fileLineBuffer,*) bufferA, bufferB, bufferC
          read(bufferA,*) crystalAxes(1,1)
          read(bufferB,*) crystalAxes(1,2)
          read(bufferC,*) crystalAxes(1,3)
          Read(101,"(A255)",IOSTAT=ios) fileRow
          fileLineBuffer = fileRow(24:100)
          read(fileLineBuffer,*) bufferA, bufferB, bufferC
          read(bufferA,*) crystalAxes(2,1)
          read(bufferB,*) crystalAxes(2,2)
          read(bufferC,*) crystalAxes(2,3)
          Read(101,"(A255)",IOSTAT=ios) fileRow
          fileLineBuffer = fileRow(24:100)
          read(fileLineBuffer,*) bufferA, bufferB, bufferC
          read(bufferA,*) crystalAxes(3,1)
          read(bufferB,*) crystalAxes(3,2)
          read(bufferC,*) crystalAxes(3,3)
        End If
! Crystal Axes/Unit Vector - vc-relax
        If(readType.eq.1.and.fileRow(1:15).eq."CELL_PARAMETERS")Then
          Read(101,"(A255)",IOSTAT=ios) fileRow
          fileLineBuffer = fileRow(1:100)
          read(fileLineBuffer,*) bufferA, bufferB, bufferC
          read(bufferA,*) crystalAxes(1,1)
          read(bufferB,*) crystalAxes(1,2)
          read(bufferC,*) crystalAxes(1,3)
          Read(101,"(A255)",IOSTAT=ios) fileRow
          fileLineBuffer = fileRow(1:100)
          read(fileLineBuffer,*) bufferA, bufferB, bufferC
          read(bufferA,*) crystalAxes(2,1)
          read(bufferB,*) crystalAxes(2,2)
          read(bufferC,*) crystalAxes(2,3)
          Read(101,"(A255)",IOSTAT=ios) fileRow
          fileLineBuffer = fileRow(1:100)
          read(fileLineBuffer,*) bufferA, bufferB, bufferC
          read(bufferA,*) crystalAxes(3,1)
          read(bufferB,*) crystalAxes(3,2)
          read(bufferC,*) crystalAxes(3,3)
        End If
! Symbol/Coords - scf
        If(readType.eq.0.and.fileRow(1:22).eq."     site n.     atom ")Then
          Do j=1,numberOfAtoms
            Read(101,"(A255)",IOSTAT=ios) fileRow
            atomType(j) = TrimSpaces(fileRow(11:23))
            Do k=1,maxFileRows
              labelTemp = configLabelReplace(k,1)
              If(labelTemp(1:1).eq." ")Then
                EXIT
              End If
              If(Adjustl(StrToUpper(labelTemp)).eq.&
                Adjustl(StrToUpper(atomType(j))))Then
                atomType(j) = Adjustl(StrToUpper(labelTemp))
              End If
            End Do
            fileLineBuffer = fileRow(39:75)
            read(fileLineBuffer,*) bufferA, bufferB, bufferC
            read(bufferA,*) atomCoords(j,1)
            read(bufferB,*) atomCoords(j,2)
            read(bufferC,*) atomCoords(j,3)
          End Do
        End If
! Symbol/Coords - vcrelax
        If(fileRow(1:16).eq."ATOMIC_POSITIONS")Then
          Do j=1,numberOfAtoms
            Read(101,"(A255)",IOSTAT=ios) fileRow
            atomType(j) = TrimSpaces(fileRow(1:7))
            Do k=1,maxFileRows
              labelTemp = configLabelReplace(k,1)
              If(labelTemp(1:1).eq." ")Then
                EXIT
              End If
              If(Adjustl(StrToUpper(labelTemp)).eq.&
                Adjustl(StrToUpper(atomType(j))))Then
                atomType(j) = Adjustl(StrToUpper(labelTemp))
              End If
            End Do
            fileLineBuffer = fileRow(8:48)
            read(fileLineBuffer,*) bufferA, bufferB, bufferC
            read(bufferA,*) atomCoords(j,1)
            read(bufferB,*) atomCoords(j,2)
            read(bufferC,*) atomCoords(j,3)
          End Do
        End If
! atom forces
        If(fileRow(1:27).eq."     Forces acting on atoms")Then
          Read(101,"(A255)",IOSTAT=ios) fileRow !read blank line
          Do j=1,numberOfAtoms
            Read(101,"(A255)",IOSTAT=ios) fileRow
            fileLineBuffer = fileRow(33:76)
            read(fileLineBuffer,*) bufferA, bufferB, bufferC
            read(bufferA,*) atomForcess(j,1)
            read(bufferB,*) atomForcess(j,2)
            read(bufferC,*) atomForcess(j,3)
          End Do
          Do j=1,numberOfAtoms
            Do k=1,3
              atomForcess(j,k) = UnitConvert(1.0D0*atomForcess(j,k),"RYBOHR","EVANG")
            End Do
          End Do
        End If
      End If
    End Do
    Close(101)
! Sort out energy
    energyOffset = 0.0D0
    Do i=1,numberOfAtoms
      atomTypeTemp = StrToUpper(Trim(atomType(i)))
      Do j=1,size(dftElement)
        If(dftElement(j).eq."  ")Then
          Exit
        End If
        If(dftElement(j).eq.atomTypeTemp(1:2))Then
          energyOffset = energyOffset - (dftOptEnergy(j)-dftCohEnergy(j))
        End If
      End Do
    End Do
    configTotalEnergy = configTotalEnergy + energyOffset
    configEnergyPerAtom = configTotalEnergy / (1.0D0 * numberOfAtoms)
! Check next free slot
    m = 0
    Do i=1,maxFileRows
      m = m + 1
      fileRow = configInputDataDFTTemp(i)     !read line
! Break out
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
    End Do
! Save to file
    fileRow = BlankString(fileRow)
    write(fileRow,"(A27,I8,A9,E20.10)") "#NEW  !added config, atoms ",numberOfAtoms,&
    ", energy ",configTotalEnergy
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4,E20.10)") "#LP ",aLat
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A3,F12.7,F12.7,F12.7)") "#X ",crystalAxes(1,1),crystalAxes(1,2),crystalAxes(1,3)
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A3,F12.7,F12.7,F12.7)") "#Y ",crystalAxes(2,1),crystalAxes(2,2),crystalAxes(2,3)
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A3,F12.7,F12.7,F12.7)") "#Z ",crystalAxes(3,1),crystalAxes(3,2),crystalAxes(3,3)
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4,F12.7,F12.7,F12.7)") "#SX ",stress(1,1),stress(1,2),stress(1,3)
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4,F12.7,F12.7,F12.7)") "#SY ",stress(2,1),stress(2,2),stress(2,3)
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4,F12.7,F12.7,F12.7)") "#SZ ",stress(3,1),stress(3,2),stress(3,3)
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A9)") "#CC 1 1 1"
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4,F10.7)") "#RC ",dftInRadiusCutoff
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4,F10.7)") "#CW ",confWeight
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A5,F10.7,A3)") "#EPA ",configEnergyPerAtom," EV"
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    If(dftInEqVol.gt.-2.0D20)Then
      fileRow = BlankString(fileRow)
      write(fileRow,"(A4,F16.7,A5)") "#EV ",dftInEqVol," ANG3"
      configInputDataDFTTemp(m) = trim(fileRow)
      m = m + 1
    End If
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4)") "#F Y"
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
!
    Do i=1,numberOfAtoms
      fileRow = BlankString(fileRow)
      write(fileRow,"(A3,A2,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10)") &
      atomType(i),"  ",atomCoords(i,1)," ",atomCoords(i,2)," ",atomCoords(i,3),&
      " ",atomForcess(i,1)," ",atomForcess(i,2)," ",atomForcess(i,3)
      configInputDataDFTTemp(m) = trim(fileRow)
      m = m + 1
    End Do
!
    fileRow = BlankString(fileRow)
    write(fileRow,"(A4)") "#END"
    configInputDataDFTTemp(m) = trim(fileRow)
    m = m + 1
  End Subroutine readPWSCFFile

  Subroutine resetConfigVars()
! Reset arrays used to store coords etc
    Implicit None   ! Force declaration of all variables
! Reset
    configCount = 0
    configurationsI = 0
    configurationsR = 0.0D0
    coordCount = 0
    configurationCoordsKey = 0
    configurationCoordsI = 0
    configurationCoordsR = 0.0D0
    configurationForcesR = -2.1D20
    coordCountG = 0
    configurationCoordsKeyG = 0
    configurationCoordsIG = 0
    configurationCoordsRG = 0.0D0
    configVolume = 0.0D0
    configVolumeOpt = -2.1D0
  End Subroutine resetConfigVars

  Subroutine processFile()
!
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i,coordStart,coordLength
    Integer(kind=StandardInteger) :: configID
    Character(len=255) :: fileRow
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF, bufferG
! Read config file into memory
    configID = 0
    configCount = 0
    coordStart = 1
    coordLength = 0
    coordCount = 0
    Do i=1,maxFileRows
      fileRow = configInputData(i)
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
      If(fileRow(1:4).eq."#NEW")Then
        configID = configID + 1
        configCount = configCount + 1
      End If
! Reals/Doubles
      If(fileRow(1:3).eq."#LP")Then
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) configurationsR(configID,1)
      End If
      If(fileRow(1:2).eq."#X")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configurationsR(configID,2)
        Read(bufferC,*) configurationsR(configID,3)
        Read(bufferD,*) configurationsR(configID,4)
      End If
      If(fileRow(1:2).eq."#Y")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configurationsR(configID,5)
        Read(bufferC,*) configurationsR(configID,6)
        Read(bufferD,*) configurationsR(configID,7)
      End If
      If(fileRow(1:2).eq."#Z")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configurationsR(configID,8)
        Read(bufferC,*) configurationsR(configID,9)
        Read(bufferD,*) configurationsR(configID,10)
      End If
      If(fileRow(1:3).eq."#RC")Then
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) configurationsR(configID,11)
      End If
! Integers
      If(fileRow(1:3).eq."#CC")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configurationsI(configID,1)
        Read(bufferC,*) configurationsI(configID,2)
        Read(bufferD,*) configurationsI(configID,3)
      End If
      If(fileRow(1:2).eq."#F")Then
        Read(fileRow,*) bufferA, bufferB
        If(bufferB(1:1).eq."Y")Then
          configurationsI(configID,4) = 1
        Else
          configurationsI(configID,4) = 0
        End If
      End If
! Co-ordinates
      If(fileRow(1:1).ne."#")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Call AddUniqueElement(bufferA(1:2))
        If(QueryUniqueElement(bufferA(1:2)).gt.0)Then
          coordCount = coordCount + 1
          coordLength = coordLength + 1
          configurationCoordsI(coordCount,1) = QueryUniqueElement(bufferA(1:2))
          Read(bufferB,*) configurationCoordsR(coordCount,1)
          Read(bufferC,*) configurationCoordsR(coordCount,2)
          Read(bufferD,*) configurationCoordsR(coordCount,3)
! Read in forces
          If(configurationsI(configID,4).eq.1)Then
            Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF, bufferG
            Read(bufferE,*) configurationForcesR(coordCount,1)
            Read(bufferF,*) configurationForcesR(coordCount,2)
            Read(bufferG,*) configurationForcesR(coordCount,3)
          End If
        End If
      End If
! Read in configuration weighting
      If(fileRow(1:3).eq."#CW")Then
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) configWeighting(configID)
      End If
! Read in energy
      If(fileRow(1:4).eq."#EPA")Then
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) configRef(configID,1)
        configRef(configID,1) = UnitConvert(configRef(configID,1),bufferC,"EV")
        configRefEnergies(configID) = configRef(configID,1)  ! Also store in energies array
      End If
! Read in EqVol
      If(fileRow(1:3).eq."#EV")Then
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) configRef(configID,2)
        configRef(configID,2) = UnitConvert(configRef(configID,2),bufferC,"ANG3")
        configRefEV(configID) = configRef(configID,2)  ! Also store in eqvol array
      End If
! Read in Bulk Modulus
      If(fileRow(1:3).eq."#BM")Then
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) configRefBM(configID)
        configRefBM(configID) = UnitConvert(configRefBM(configID),bufferC,"GPA")
      End If
! Read in stresses
      If(fileRow(1:3).eq."#SX")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configRefStresses(configID,1)
        Read(bufferC,*) configRefStresses(configID,2)
        Read(bufferD,*) configRefStresses(configID,3)
      End If
      If(fileRow(1:3).eq."#SY")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configRefStresses(configID,4)
        Read(bufferC,*) configRefStresses(configID,5)
        Read(bufferD,*) configRefStresses(configID,6)
      End If
      If(fileRow(1:3).eq."#SZ")Then
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferB,*) configRefStresses(configID,7)
        Read(bufferC,*) configRefStresses(configID,8)
        Read(bufferD,*) configRefStresses(configID,9)
      End If
! End - store start/length/end coord list
      If(fileRow(1:4).eq."#END")Then
        configurationCoordsKey(configID,1) = coordStart
        configurationCoordsKey(configID,2) = coordLength
        configurationCoordsKey(configID,3) = coordStart+coordLength-1
        coordStart = coordStart + coordLength
        If(TerminalPrint())Then
          print *,"Loaded: ",configID,coordStart,(coordStart+coordLength-1),"(",coordLength,")"
        End If
        coordLength = 0
      End If
    End Do
  End Subroutine processFile

  Subroutine expandCoordinates()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID,i
    Integer(kind=StandardInteger) :: coordStart, coordLength, coordEnd
    Integer(kind=StandardInteger) :: coordG, coordStartG, coordLengthG, coordEndG
    Integer(kind=StandardInteger) :: xLoop, yLoop, zLoop
    Integer(kind=StandardInteger) :: xCopies, yCopies, zCopies
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: crystalUnitCellTemp
    Real(kind=DoubleReal), Dimension(1:3) :: aVect
! Blank arrays
    configurationCoordsKeyG = 0
    configurationCoordsIG = 0
    configurationCoordsRG = 0.0D0
! Initialise counters
    coordG = 0
    coordStartG = 1
    coordLengthG = 0
    coordEndG = 0
    configsAtomTotal = 0
! Loop through configs
    Do configID=1,configCount
! Clear temp unit cell matrix
      crystalUnitCellTemp = 0.0D0
! Store user input
      crystalUnitCellTemp(1,1) = configurationsR(configID,2)   !xx
      crystalUnitCellTemp(1,2) = configurationsR(configID,3)   !xy
      crystalUnitCellTemp(1,3) = configurationsR(configID,4)   !xz
      crystalUnitCellTemp(2,1) = configurationsR(configID,5)   !yx
      crystalUnitCellTemp(2,2) = configurationsR(configID,6)   !yy
      crystalUnitCellTemp(2,3) = configurationsR(configID,7)   !yz
      crystalUnitCellTemp(3,1) = configurationsR(configID,8)   !zx
      crystalUnitCellTemp(3,2) = configurationsR(configID,9)   !zy
      crystalUnitCellTemp(3,3) = configurationsR(configID,10)   !zz
! Multiply by global unit cell vector
      crystalUnitCellTemp = MatMul(globalConfigUnitVector,crystalUnitCellTemp)
! Lattice parameter configurationsR(configID,1)
! copies x dir      configurationsI(configID,1)
! copies y dir      configurationsI(configID,2)
! copies z dir      configurationsI(configID,3)
      xCopies = configurationsI(configID,1)
      yCopies = configurationsI(configID,2)
      zCopies = configurationsI(configID,3)
! x vector
      crystalUnitCellTemp(1,1) = crystalUnitCellTemp(1,1) * &
      configurationsR(configID,1) * configurationsI(configID,1)
      crystalUnitCellTemp(1,2) = crystalUnitCellTemp(1,2) * &
      configurationsR(configID,1) * configurationsI(configID,2)
      crystalUnitCellTemp(1,3) = crystalUnitCellTemp(1,3) * &
      configurationsR(configID,1) * configurationsI(configID,3)
! y vector
      crystalUnitCellTemp(2,1) = crystalUnitCellTemp(2,1) * &
      configurationsR(configID,1) * configurationsI(configID,1)
      crystalUnitCellTemp(2,2) = crystalUnitCellTemp(2,2) * &
      configurationsR(configID,1) * configurationsI(configID,2)
      crystalUnitCellTemp(2,3) = crystalUnitCellTemp(2,3) * &
      configurationsR(configID,1) * configurationsI(configID,3)
! z vector
      crystalUnitCellTemp(3,1) = crystalUnitCellTemp(3,1) * &
      configurationsR(configID,1) * configurationsI(configID,1)
      crystalUnitCellTemp(3,2) = crystalUnitCellTemp(3,2) * &
      configurationsR(configID,1) * configurationsI(configID,2)
      crystalUnitCellTemp(3,3) = crystalUnitCellTemp(3,3) * &
      configurationsR(configID,1) * configurationsI(configID,3)
! store
      crystalUnitCell(configID,1) = crystalUnitCellTemp(1,1)
      crystalUnitCell(configID,2) = crystalUnitCellTemp(1,2)
      crystalUnitCell(configID,3) = crystalUnitCellTemp(1,3)
      crystalUnitCell(configID,4) = crystalUnitCellTemp(2,1)
      crystalUnitCell(configID,5) = crystalUnitCellTemp(2,2)
      crystalUnitCell(configID,6) = crystalUnitCellTemp(2,3)
      crystalUnitCell(configID,7) = crystalUnitCellTemp(3,1)
      crystalUnitCell(configID,8) = crystalUnitCellTemp(3,2)
      crystalUnitCell(configID,9) = crystalUnitCellTemp(3,3)
! Volume
      configVolume(configID) = TripleProductSq(crystalUnitCellTemp)
! Generate co-ordinates
      coordStart = configurationCoordsKey(configID,1)
      coordLength = configurationCoordsKey(configID,2)
      coordEnd = configurationCoordsKey(configID,3)
! If(mpiProcessID.eq.0)Then
! print *,"input",coordStart,coordLength,coordEnd
! End If
      Do xLoop=1,xCopies
        Do yLoop=1,yCopies
          Do zLoop=1,zCopies
            Do i=coordStart,coordEnd
              coordG = coordG + 1
! Atom label
              configurationCoordsIG(coordG,1) = configurationCoordsI(i,1)
              configurationCoordsRG(coordG,1) = &
              (xLoop + configurationCoordsR(i,1) - 1.0D0)/(1.0D0*xCopies)
              configurationCoordsRG(coordG,2) = &
              (yLoop + configurationCoordsR(i,2) - 1.0D0)/(1.0D0*yCopies)
              configurationCoordsRG(coordG,3) = &
              (zLoop + configurationCoordsR(i,3) - 1.0D0)/(1.0D0*zCopies)
! Forces
              If(configurationsI(configID,4).eq.0)Then
                configRefForces(coordG,1) = -2.1D20
                configRefForces(coordG,2) = -2.1D20
                configRefForces(coordG,3) = -2.1D20
              Else
                configRefForces(coordG,1) = configurationForcesR(i,1)
                configRefForces(coordG,2) = configurationForcesR(i,2)
                configRefForces(coordG,3) = configurationForcesR(i,3)
              End If
            End Do
          End Do
        End Do
      End Do
      coordEndG = coordG
      coordLengthG = coordEndG - coordStartG + 1
      configsAtomTotal = configsAtomTotal + coordLengthG
      configurationCoordsKeyG(configID,1) = coordStartG
      configurationCoordsKeyG(configID,2) = coordLengthG
      configurationCoordsKeyG(configID,3) = coordEndG
! reset start coord count
      coordStartG = coordEndG + 1
    End Do
! Transform coordinates
    Do configID=1,configCount
! Load transformation matrix
      crystalUnitCellTemp(1,1) = crystalUnitCell(configID,1)
      crystalUnitCellTemp(1,2) = crystalUnitCell(configID,2)
      crystalUnitCellTemp(1,3) = crystalUnitCell(configID,3)
      crystalUnitCellTemp(2,1) = crystalUnitCell(configID,4)
      crystalUnitCellTemp(2,2) = crystalUnitCell(configID,5)
      crystalUnitCellTemp(2,3) = crystalUnitCell(configID,6)
      crystalUnitCellTemp(3,1) = crystalUnitCell(configID,7)
      crystalUnitCellTemp(3,2) = crystalUnitCell(configID,8)
      crystalUnitCellTemp(3,3) = crystalUnitCell(configID,9)
! Load start-end
      coordStartG = configurationCoordsKeyG(configID,1)
      coordEndG = configurationCoordsKeyG(configID,3)
      Do coordG=coordStartG,coordEndG
! Store fractional values
        configurationCoordsRG(coordG,4) = 1.0D0*configurationCoordsRG(coordG,1)
        configurationCoordsRG(coordG,5) = 1.0D0*configurationCoordsRG(coordG,2)
        configurationCoordsRG(coordG,6) = 1.0D0*configurationCoordsRG(coordG,3)
! Store in arrays
        aVect(1) = 1.0D0*configurationCoordsRG(coordG,1)
        aVect(2) = 1.0D0*configurationCoordsRG(coordG,2)
        aVect(3) = 1.0D0*configurationCoordsRG(coordG,3)
! Transform coords
        aVect = TransformCoords(aVect,crystalUnitCellTemp)
! Set "real" co-ordinates
        configurationCoordsRG(coordG,1) = 1.0D0*aVect(1)
        configurationCoordsRG(coordG,2) = 1.0D0*aVect(2)
        configurationCoordsRG(coordG,3) = 1.0D0*aVect(3)
      End Do
    End Do
  End Subroutine expandCoordinates

  Subroutine AddUniqueElement(element)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: i, k, found
! convert to uppercase
    element = adjustl(StrToUpper(element))  ! Force uppercase, align to the left
! loop through elements array
    k = 0
    found = 0
    Do i=1,size(elements,1)
      k = k + 1
      If(elements(i).eq."ZZ")Then
        exit
      End If
      If(elements(i).eq.element)Then
        found = 1
        exit
      End If
    End Do
! save element If not found
    If(found.eq.0)Then
      elements(k) = element
! save charge
      Do i=0,size(elementSymbol,1)-1
        If(element.eq.elementSymbol(i))Then
          elementsCharge(k) = i
        End If
      End Do
    End If
  End Subroutine AddUniqueElement
  
  
  
  
  
  
  Subroutine readBpConfigFile()
    Implicit None   ! Force declaration of all variables
! Private
    Integer(kind=StandardInteger) :: i, fileRows, rowID, configID
    Integer(kind=StandardInteger) :: intTemp
    Real(kind=DoubleReal) :: dpTemp
    Character(Len=255) :: fileRowUC
    Character(len=32) :: bufferA, bufferB, bufferC
    Logical :: readData    
! Print out
    If(TerminalPrint())Then
      print *,"READ BP CONFIG: ",trim(bpConfigFilePath)
    End If
! Private variables
    Call cpu_time(timeStart)
! Load config file into memory
    Call readFile(trim(bpConfigFilePath), bpConfigInputData, fileRows)
! Print to terminal    
    If(TerminalPrint())Then
      print *,"Rows: ",fileRows
    End If  
! Read in data
    configID = 0
    rowID = 0
    readData = .false.
    Do i=1,fileRows
      fileRowUC = strtoupper(trim(bpConfigInputData(i)))
      If(fileRowUC(1:4).eq."#NEW")Then
        readData = .true.
        configID = configID + 1
      End If
      If(fileRowUC(1:4).eq."#END")Then
        readData = .false.
! ----- Fill in gaps
! Calculate bulk properties V0
        If(bpInArr(configID)%structure.eq."FCC")Then
          bpInArr(configID)%v0 = ((bpInArr(configID)%aLat)**3)/4.0D0
        End If  
        If(bpInArr(configID)%structure.eq."BCC")Then
          bpInArr(configID)%v0 = ((bpInArr(configID)%aLat)**3)/2.0D0
        End If  
        !If(bpInArr(configID)%structure.eq."HCP")Then
          !bpInArr(configID)%v0 = ((bpInArr(configID)%aLat)**3)/2.0D0
        !End If  
      End If
      If(readData)Then
        If(fileRowUC(1:5).eq."#STRU")Then        
          Read(fileRowUC,*) bufferA, bufferB
          bpInArr(configID)%structure = bufferB(1:3)          
        End If
        If(fileRowUC(1:5).eq."#SIZE")Then        
          Read(fileRowUC,*) bufferA, bufferB
          Read(bufferB,*) intTemp
          bpInArr(configID)%size = intTemp          
        End If
        If(fileRowUC(1:4).eq."#ELE")Then        
          Read(fileRowUC,*) bufferA, bufferB
          bpInArr(configID)%element = bufferB(1:2)          
        End If
        If(fileRowUC(1:5).eq."#ALAT")Then        
          Read(fileRowUC,*) bufferA, bufferB, bufferC
          Read(bufferB,*) dpTemp
          bpInArr(configID)%alat = UnitConvert(dpTemp, bufferC, "ANGS")
        End If
        If(fileRowUC(1:3).eq."#E0")Then        
          Read(fileRowUC,*) bufferA, bufferB, bufferC
          Read(bufferB,*) dpTemp
          bpInArr(configID)%e0 = UnitConvert(dpTemp, bufferC, "EV")          
        End If
        If(fileRowUC(1:3).eq."#B0")Then        
          Read(fileRowUC,*) bufferA, bufferB, bufferC
          Read(bufferB,*) dpTemp
          bpInArr(configID)%b0 = UnitConvert(dpTemp, bufferC, "EVAN3")          
        End If
        If(fileRowUC(1:4).eq."#BP0")Then        
          Read(fileRowUC,*) bufferA, bufferB
          Read(bufferB,*) dpTemp
          bpInArr(configID)%bp0 = dpTemp          
        End If
        If(fileRowUC(1:4).eq."#C11")Then        
          Read(fileRowUC,*) bufferA, bufferB, bufferC
          Read(bufferB,*) dpTemp
          bpInArr(configID)%c11 = UnitConvert(dpTemp, bufferC, "EVAN3")         
        End If
        If(fileRowUC(1:4).eq."#C12")Then        
          Read(fileRowUC,*) bufferA, bufferB, bufferC
          Read(bufferB,*) dpTemp
          bpInArr(configID)%c12 = UnitConvert(dpTemp, bufferC, "EVAN3")         
        End If
        If(fileRowUC(1:4).eq."#C44")Then        
          Read(fileRowUC,*) bufferA, bufferB, bufferC
          Read(bufferB,*) dpTemp
          bpInArr(configID)%c44 = UnitConvert(dpTemp, bufferC, "EVAN3")         
        End If
      End If
    End Do

    

! Synch MPI processes
    Call M_synchProcesses()
    Call cpu_time(timeEnd)
    
    
! Output
    If(TerminalPrint())Then
      print *,"Atom bp configurations loaded: ",(timeEnd-timeStart),"s"
    End If
  End Subroutine readBpConfigFile
  
  

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!
  Function QueryUniqueElement (element) RESULT (output)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: output
    Integer(kind=StandardInteger) :: i, k, found
! convert to uppercase
    element = adjustl(StrToUpper(element))
! loop through elements array
    k = 0
    found = 0
    Do i=1,size(elements,1)
      k = k + 1
      If(elements(i).eq.element)Then
        found = 1
        exit
      End If
    End Do
! save element if not found
    If(found.eq.1)Then
      output = k
    Else
      output = 0
    End If
  End Function QueryUniqueElement

  Function QueryFunctionType (functionType) RESULT (output)
    Character(len=4) :: functionType
    Integer(kind=StandardInteger) :: output
    Integer(kind=StandardInteger) :: i
! convert to uppercase
    FunctionType = adjustl(StrToUpper(functionType))
    Do i=1,size(eamFunctionTypes,1)
      If(functionType.eq.eamFunctionTypes(i))Then
        output = i
        Exit
      End If
    End Do
  End Function QueryFunctionType

End Module readConfig
