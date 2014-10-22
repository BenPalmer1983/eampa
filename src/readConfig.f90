Module readConfig

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
  Use initialise
  Use loadData   
  Use globals
  Use output 
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: readConfigFile
  Public :: generateCoords
  
Contains
  Subroutine readConfigFile(resetVarsIn,readIDFromIn)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger), optional :: resetVarsIn, readIDFromIn
    Integer(kind=StandardInteger) :: resetVars, readIDFrom
    Real(kind=DoubleReal) :: timeStartRC, timeEndRC
! Optional variables    
    resetVars = 0                   ! Whether to clear out the config variables/arrays or not
    If(present(resetVarsIn))Then
      resetVars = resetVarsIn
    End If
    readIDFrom = 1                  ! Which ID config to read into
    If(present(readIDFromIn))Then
      readIDFrom = readIDFromIn
    End If
! Start Time
    Call cpu_time(timeStartRC) 
! Clear config arrays
    If(resetVars.eq.1)Then
      Call resetConfigVars()
    End If  
! Prepare the temporary config file
    If(mpiProcessID.eq.0)Then
      Call prepFile()   
      Call readDFTFiles()      
    End If
! Synch MPI processes    
    Call M_synchProcesses()
! Send temp file name to all processes
    Call M_distChar(configFilePathT)
! Read the temporary config file
    Call readFile(readIDFrom) 
! Generate co-ordinates   
    Call generateCoords(readIDFrom)
! Synch MPI processes    
    Call M_synchProcesses()
! Output summary of config to the output file
    If(readIDFrom.eq.1)Then
      Call outputConfigSummaryT()  
    End If  
! Output files
    If(mpiProcessID.eq.0)Then
      Call outputConfigFiles()   
    End If   
! Synch MPI processes    
    Call M_synchProcesses() 
! Remove file
    Call rmFile(configFilePathT)
! End Time
    Call cpu_time(timeEndRC)        
! Store Time    
    Call storeTime(7,timeEndRC-timeStartRC)  
  End Subroutine readConfigFile 
    
  
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
  
  Subroutine prepFile()
! Creates a temporary file to contain the config potential
! Converts alpha to upper, and strips out comment lines etc
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i
    Character(len=255) :: fileRow
    Character(len=255) :: configFilePathTA
    Character(len=8) :: tempName
! Set temp file path 
    Call tempFileName(tempName)
    configFilePathT = Trim(tempDirectory)//"/"//tempName
    configFilePathTA = Trim(configFilePathT)//".temp.conf"
! Output to Terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      !Print *,"Reading user config file ",trim(configFilePath)
    End If  
! Read EAM file and make new temp EAM pot file
    Open(UNIT=1,FILE=Trim(configFilePath)) 
    Open(UNIT=2,FILE=Trim(configFilePathTA)) 
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If 
      fileRow = RemoveComments(fileRow)
      fileRow = RemoveQuotes(fileRow)
      fileRow = Trim(Adjustl(fileRow))
      If(fileRow(1:1).ne."!".or.fileRow(1:1).ne." ")Then
        If(StrToUpper(fileRow(1:5)).eq."#PATH")Then  !Ensure path remains case sensitive
          fileRow = "#PATH"//fileRow(6:255)
          write(2,"(A)") Trim(fileRow)
        Else
          write(2,"(A)") StrToUpper(Trim(fileRow))
        End If
      End If
    End Do
    Close(2)
    Close(1)
  End Subroutine prepFile 
  
  
  Subroutine readDFTFiles()
! Creates a temporary file to contain the config potential
! Converts alpha to upper, and strips out comment lines etc
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i, j, writeFile
    Real(kind=DoubleReal) :: radiusCutoff, eqVol, confWeight
    Character(len=255) :: fileRow
    Character(len=255) :: configFilePathTA, configFilePathTB, configFilePathTC, configFilePathTDFT
    Character(len=255) :: dftFilePath
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD
    Character(len=32) :: dftType
! Init variables
    bufferA = BlankString(bufferA)
    bufferB = BlankString(bufferB)
    bufferC = BlankString(bufferC)
    bufferD = BlankString(bufferD)
    configFilePathTA = Trim(configFilePathT)//".temp.conf"    ! Prepared file
    configFilePathTB = Trim(configFilePathT)//"B.temp.conf"   ! Input configs
    configFilePathTC = Trim(configFilePathT)//"C.temp.conf"   ! DFT input configs
    configFilePathTDFT = Trim(configFilePathT)//"DFT.temp.conf" ! DFT file
    confWeight = 1.0D0
! Separate configs from dft file configs
    Open(UNIT=1,FILE=Trim(configFilePathTA)) 
    Open(UNIT=2,FILE=Trim(configFilePathTB)) 
    Open(UNIT=3,FILE=Trim(configFilePathTC)) 
    writeFile = 0
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If      
      If(fileRow(1:7).eq."#NEWDFT")Then
        writeFile = 3
      Elseif(fileRow(1:4).eq."#NEW")Then
        writeFile = 2
      End If
      If(writeFile.gt.0)Then
        write(writeFile,"(A)") Trim(fileRow)   
      End If
    End Do
    Close(3)
    Close(2)
    Close(1)
! Set temp file name    
    configFilePathT = Trim(configFilePathT)//".tempfinal.conf" 
! Copy already set configs to config temp file 
    Open(UNIT=1,FILE=Trim(configFilePathTB))    
    Open(UNIT=2,FILE=Trim(configFilePathT)) 
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If  
      write(2,"(A)") Trim(fileRow)
    End Do  
    Close(2)
    Close(1)
! Loop through DFT files    
    Open(UNIT=1,FILE=Trim(configFilePathTC))
    Do i=1,maxFileRows   
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If  
      If(fileRow(1:7).eq."#NEWDFT")Then  ! Clean variables
        dftFilePath = BlankString(dftFilePath)
        dftType = BlankString(dftType)
        dftReplaceLabel = BlankString2DArray(dftReplaceLabel)
        eqVol = -2.1D20
      End If
      If(fileRow(1:5).eq."#PATH")Then  ! Path to DFT file     
        Read(fileRow,*) bufferA, dftFilePath
      End If
      If(fileRow(1:5).eq."#TYPE")Then  ! Ab init file type to read in
        Read(fileRow,*) bufferA, dftType
      End If
      If(fileRow(1:3).eq."#EV")Then  ! Equilibrium volume
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) eqVol   
        eqVol = UnitConvert(eqVol, bufferC, "ANG3")
      End If
      If(fileRow(1:3).eq."#CW")Then  ! Config Weighting
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) confWeight   
      End If
      If(fileRow(1:3).eq."#RC")Then  ! Radius cutoff 
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) radiusCutoff
        radiusCutoff = UnitConvert(radiusCutoff, bufferC, "ANGS")
      End If
      If(fileRow(1:5).eq."#REPL")Then  ! Replace label
        Read(fileRow,*) bufferA, bufferB, bufferC
        Do j=1,size(dftReplaceLabel,1)
          If(dftReplaceLabel(j,1).eq."        ")Then
            dftReplaceLabel(j,1) = trim(bufferB)
            dftReplaceLabel(j,2) = trim(bufferC)
            Exit
          End If
        End Do
      End If
      If(fileRow(1:7).eq."#ENDDFT")Then  ! Clean variables
        If(trim(dftType).eq."PWSCF")Then
          Call readPWSCFFile(dftFilePath, configFilePathT, &
          eqVol, radiusCutoff, confWeight)
        End If  
      End If
    End Do
    Close(1)
! Remove unnecessary files    
    Call rmFile(configFilePathTA)
    Call rmFile(configFilePathTB)
    Call rmFile(configFilePathTC)
    Call rmFile(configFilePathTDFT)
    
! Add files to clean
    !Call fileToClean(configFilePathT)
    !Call fileToClean(configFilePathTA)
    !Call fileToClean(configFilePathTB)
    !Call fileToClean(configFilePathTC)
    !Call fileToClean(configFilePathTDFT)
    !Call rmFile(configFilePathT)
    !Call rmFile(configFilePathTA)
    !Call rmFile(configFilePathTB)
    !Call rmFile(configFilePathTC)
    !Call rmFile(configFilePathTDFT)
  End Subroutine readDFTFiles 
  
  
  Subroutine readFile(readIDFrom)
! Reads in the Config potential from the temporary Config file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i, coordStart, coordLength
    Integer(kind=StandardInteger) :: configID, readIDFrom
    Character(len=255) :: fileRow
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF, bufferG
! Initialise variables
    fileRow = BlankString(fileRow)
    bufferA = BlankString(bufferA)
    bufferB = BlankString(bufferB)
    bufferC = BlankString(bufferC)
    bufferD = BlankString(bufferD)  
    bufferE = BlankString(bufferE) 
    bufferF = BlankString(bufferF) 
    bufferG = BlankString(bufferG)  
    configID = readIDFrom-1
    configCountRI = 0
! coord start point
    If(readIDFrom.eq.1)Then
      coordCount = 0
      coordStart = 1  
      coordLength = 0 
    Else
      If(configCount.eq.0)Then
        coordCount = 0
        coordStart = 1  
        coordLength = 0 
      Else
        coordCount = configurationCoordsKey(configCount,3)
        coordStart = configurationCoordsKey(configCount,3)+1
        coordLength = 0 
      End If  
    End If    
! Read Config file and make new temp Config pot file
    Open(UNIT=1,FILE=Trim(configFilePathT)) 
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If 
      fileRow = Trim(Adjustl(fileRow))
      If(fileRow(1:4).eq."#NEW")Then
        configID = configID + 1      
        configCountRI = configCountRI + 1        
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
        coordLength = 0
      End If
    End Do
    Close(1) 
! Number of configurations read in
    If(readIDFrom.eq.1)Then   !If configs read into slot 1 upwards, count as total configs
      configCount = configCountRI
    End If    
  End Subroutine readFile 
  

  
  
  Subroutine generateCoords(readIDFrom)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: readIDFrom
    Integer(kind=StandardInteger) :: i, j, coordStart, coordEnd, coordLast
    Integer(kind=StandardInteger) :: coordStartG, coordLengthG
    Integer(kind=StandardInteger) :: x, y, z, xCopy, yCopy, zCopy
    Real(kind=DoubleReal) :: aLat
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: configUnitVector 
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: configVolVector 
    Real(kind=DoubleReal) :: xC, yC, zC, xCa, yCa, zCa    
! Init variables
    configUnitVector = 0.0D0
    coordStart = 0
    coordEnd = 0
    xCopy = 0
    yCopy = 0
    zCopy = 0
    aLat = 0.0D0
! coord start point
    If(readIDFrom.eq.1)Then
      coordCountG = 0
      coordLengthG = 0
      coordStartG = 1
    Else  
      If(configCount.eq.0)Then
        coordCountG = 0
        coordLengthG = 0
        coordStartG = 1
      Else
        coordLast = 0
        Do i=1,(readIDFrom-1)
          If(coordLast.lt.configurationCoordsKeyG(i,3))Then
            coordLast = configurationCoordsKeyG(i,3)
          End If  
        End Do
        coordCountG = coordLast   
        coordLengthG = 0   
        coordStartG = coordLast+1
      End If  
    End If  
! Loop through configurations
    !print *,readIDFrom,(readIDFrom+configCountRI-1)
    Do i=readIDFrom,(readIDFrom+configCountRI-1)
! Set config variables    
      coordStart = configurationCoordsKey(i,1)
      coordEnd = configurationCoordsKey(i,3)
      xCopy = configurationsI(i,1) 
      yCopy = configurationsI(i,2) 
      zCopy = configurationsI(i,3) 
      aLat = configurationsR(i,1)    
! Get config unit vector
      configUnitVector(1,1) = configurationsR(i,2) 
      configUnitVector(1,2) = configurationsR(i,3) 
      configUnitVector(1,3) = configurationsR(i,4) 
      configUnitVector(2,1) = configurationsR(i,5) 
      configUnitVector(2,2) = configurationsR(i,6) 
      configUnitVector(2,3) = configurationsR(i,7) 
      configUnitVector(3,1) = configurationsR(i,8) 
      configUnitVector(3,2) = configurationsR(i,9) 
      configUnitVector(3,3) = configurationsR(i,10) 
! Apply global unit vector      
      configUnitVector = matmul(globalConfigUnitVector,configUnitVector)
! Store config unit vector      
      configurationsR(i,21) = configUnitVector(1,1)
      configurationsR(i,22) = configUnitVector(1,2)
      configurationsR(i,23) = configUnitVector(1,3)
      configurationsR(i,24) = configUnitVector(2,1)
      configurationsR(i,25) = configUnitVector(2,2)
      configurationsR(i,26) = configUnitVector(2,3)
      configurationsR(i,27) = configUnitVector(3,1)
      configurationsR(i,28) = configUnitVector(3,2)
      configurationsR(i,29) = configUnitVector(3,3)
! Loop through coords to make supercell
      Do x=1,xCopy
        Do y=1,yCopy
          Do z=1,zCopy    
            Do j=coordStart,coordEnd  
! Increment counter            
              coordCountG = coordCountG + 1
              coordLengthG = coordLengthG + 1
! Calculate coordinate in supercell
              xCa = aLat*(x + configurationCoordsR(j,1) - 1)
              yCa = aLat*(y + configurationCoordsR(j,2) - 1)
              zCa = aLat*(z + configurationCoordsR(j,3) - 1)
! Apply unit vector
              xC = xCa*configUnitVector(1,1)+yCa*configUnitVector(1,2)+zCa*configUnitVector(1,3)
              yC = xCa*configUnitVector(2,1)+yCa*configUnitVector(2,2)+zCa*configUnitVector(2,3)
              zC = xCa*configUnitVector(3,1)+yCa*configUnitVector(3,2)+zCa*configUnitVector(3,3)
! Store results
              configurationCoordsIG(coordCountG,1) = configurationCoordsI(j,1)      ! Atom type
              configurationCoordsRG(coordCountG,1) = xC                             ! X coord
              configurationCoordsRG(coordCountG,2) = yC                             ! Y coord
              configurationCoordsRG(coordCountG,3) = zC                             ! Z coord
              configRefForces(coordCountG,1) = configurationForcesR(j,1)            ! Force x dir
              configRefForces(coordCountG,2) = configurationForcesR(j,2)            ! Force y dir
              configRefForces(coordCountG,3) = configurationForcesR(j,3)            ! Force z dir          
            End Do
          End Do
        End Do
      End Do
! Store coord key      
      configurationCoordsKeyG(i,1) = coordStartG
      configurationCoordsKeyG(i,2) = coordLengthG
      configurationCoordsKeyG(i,3) = coordStartG+coordLengthG-1      
      coordStartG = coordStartG + coordLengthG
      coordLengthG = 0
! Store configuration volume   
      Do j=1,3 
        configVolVector(j,1) = aLat*xCopy*configUnitVector(j,1)
        configVolVector(j,2) = aLat*yCopy*configUnitVector(j,2)
        configVolVector(j,3) = aLat*zCopy*configUnitVector(j,3)
      End Do  
      configVolume(i) = TripleProductSq(configVolVector)
    End Do    
! Count total configurations
    configCountT = 0
    Do i=1,1024
      If(configurationCoordsKeyG(i,1).gt.0)Then
        configCountT = configCountT + 1
      End If
    End Do  
  End Subroutine generateCoords 
  
  

 
  Subroutine outputConfigFiles() 
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: i, j, coordStart, coordStartG, coordEnd, coordEndG
    Character(len=8) :: elementLabel
    If(saveConfigFile(1:1).ne." ")Then      
      Open(UNIT=102,FILE=Trim(outputDirectory)//"/"//Trim(saveConfigFile))
      Do i=1,configCount 
        write(102,"(A4)") "#NEW"    
        write(102,"(A4,E20.10)") "#LP ",configurationsR(i,1)
        write(102,"(A3,F12.7,F12.7,F12.7)") "#X ",configurationsR(i,2),&
          configurationsR(i,3),configurationsR(i,4)
        write(102,"(A3,F12.7,F12.7,F12.7)") "#Y ",configurationsR(i,5),&
          configurationsR(i,6),configurationsR(i,7)
        write(102,"(A3,F12.7,F12.7,F12.7)") "#Z ",configurationsR(i,8),&
          configurationsR(i,9),configurationsR(i,10)
        If(configRefStresses(i,1).gt.-2.0D20)Then
          write(102,"(A4,F12.7,F12.7,F12.7)") "#SX ",configRefStresses(i,1),&
          configRefStresses(i,2),configRefStresses(i,3)
          write(102,"(A4,F12.7,F12.7,F12.7)") "#SY ",configRefStresses(i,4),&
          configRefStresses(i,5),configRefStresses(i,6)
          write(102,"(A4,F12.7,F14.7,F12.7)") "#SZ ",configRefStresses(i,7),&
          configRefStresses(i,8),configRefStresses(i,9)
        End If 
        write(102,"(A4,I2,A1,I2,A1,I2)") "#CC ",configurationsI(i,1)," ",&
        configurationsI(i,2)," ",configurationsI(i,3)
        write(102,"(A4,F10.7)") "#RC ",configurationsR(i,11)
        If(configRef(i,1).gt.-2.0D20)Then
          write(102,"(A5,F10.7,A3)") "#EPA ",configRef(i,1)," EV"
        End If
        If(configRefEV(i).gt.-2.0D20)Then
          write(102,"(A4,F12.7,A5)") "#EV ",configRefEV(i)," ANG3"
        End If
        coordStart = configurationCoordsKey(i,1)
        coordEnd = configurationCoordsKey(i,3)     
        If(configurationForcesR(coordStart,1).gt.-2.0D20)Then
          write(102,"(A4)") "#F Y"
          Do j=coordStart,coordEnd
            elementLabel = elements(configurationCoordsI(j,1))
            elementLabel = SpacesRight(elementLabel)
            write(102,"(A3,A2,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10)") &
            elementLabel,"  ",&
            configurationCoordsR(j,1)," ",configurationCoordsR(j,2)," ",&
            configurationCoordsR(j,3)," ",&
            configurationForcesR(j,1)," ",configurationForcesR(j,2)," ",configurationForcesR(j,3)
          End Do
        Else
          write(102,"(A4)") "#F N"
          Do j=coordStart,coordEnd
            elementLabel = elements(configurationCoordsI(j,1))
            elementLabel = SpacesRight(elementLabel)
            write(102,"(A3,A2,F16.10,A1,F16.10,A1,F16.10)") &
            elementLabel,"  ",&
            configurationCoordsR(j,1)," ",configurationCoordsR(j,2)," ",configurationCoordsR(j,3)
          End Do
        End If
        write(102,"(A4)") "#END" 
      End Do
      Close(102)
    End If
    If(saveExpConfigFile(1:1).ne." ")Then
      Open(UNIT=102,FILE=Trim(outputDirectory)//"/"//Trim(saveExpConfigFile))
      Do i=1,configCount 
        write(102,"(A4)") "#NEW"    
        write(102,"(A4,E20.10)") "#LP  ",(configurationsR(i,1)*configurationsI(i,1))
        write(102,"(A3,F12.7,F12.7,F12.7)") "#X ",configurationsR(i,2),&
          configurationsR(i,3),configurationsR(i,4)
        write(102,"(A3,F12.7,F12.7,F12.7)") "#Y ",configurationsR(i,5),&
          configurationsR(i,6),configurationsR(i,7)
        write(102,"(A3,F12.7,F12.7,F12.7)") "#Z ",configurationsR(i,8),&
          configurationsR(i,9),configurationsR(i,10)
        If(configRefStresses(i,1).gt.-2.0D20)Then
          write(102,"(A4,F12.7,F12.7,F12.7)") "#SX ",configRefStresses(i,1),&
          configRefStresses(i,2),configRefStresses(i,3)
          write(102,"(A4,F12.7,F12.7,F12.7)") "#SY ",configRefStresses(i,4),&
          configRefStresses(i,5),configRefStresses(i,6)
          write(102,"(A4,F12.7,F12.7,F12.7)") "#SZ ",configRefStresses(i,7),&
          configRefStresses(i,8),configRefStresses(i,9)
        End If 
        write(102,"(A12)") "#CC    1 1 1"
        !configurationsR(i,1)
        !,configurationsI(i,1)," ",&
        !configurationsI(i,2)," ",configurationsI(i,3)
        write(102,"(A4,F10.7)") "#RC ",configurationsR(i,11)
        If(configRef(i,1).gt.-2.0D20)Then
          write(102,"(A5,F10.7,A3)") "#EPA ",configRef(i,1)," EV"
        End If
        coordStartG = configurationCoordsKeyG(i,1)
        coordEndG = configurationCoordsKeyG(i,3)     
        If(configRefForces(coordStartG,1).gt.-2.0D20)Then
          write(102,"(A4)") "#F Y"
          Do j=coordStartG,coordEndG
            elementLabel = elements(configurationCoordsIG(j,1))
            elementLabel = SpacesRight(elementLabel)
            write(102,"(A3,A2,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10)") &
            elementLabel,"  ",&
            (configurationCoordsRG(j,1)/(configurationsR(i,1)*configurationsI(i,1)))," ",&
            (configurationCoordsRG(j,2)/(configurationsR(i,1)*configurationsI(i,2)))," ",&
            (configurationCoordsRG(j,3)/(configurationsR(i,1)*configurationsI(i,3)))," ",&
            configRefForces(j,1)," ",configRefForces(j,2)," ",configRefForces(j,3)
          End Do
        Else
          write(102,"(A4)") "#F N"
          Do j=coordStartG,coordEndG
            elementLabel = elements(configurationCoordsIG(j,1))
            elementLabel = SpacesRight(elementLabel)
            write(102,"(A3,A2,F16.10,A1,F16.10,A1,F16.10)") &
            elementLabel,"  ",&
            (configurationCoordsRG(j,1)/(configurationsR(i,1)*configurationsI(i,1)))," ",&
            (configurationCoordsRG(j,2)/(configurationsR(i,1)*configurationsI(i,2)))," ",&
            (configurationCoordsRG(j,3)/(configurationsR(i,1)*configurationsI(i,3)))
          End Do
        End If
        write(102,"(A4)") "#END" 
      End Do
      Close(102)
    End If
  
  End Subroutine outputConfigFiles 
  
  
  
  Subroutine AddUniqueElement(element)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: i, k, found
!convert to uppercase
    element = adjustl(StrToUpper(element))  ! Force uppercase, align to the left
!loop through elements array	
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
!save element If not found
    If(found.eq.0)Then
      elements(k) = element
!save charge
      Do i=0,size(elementSymbol,1)-1
        If(element.eq.elementSymbol(i))Then
          elementsCharge(k) = i
        End If
      End Do      
    End If    
  End Subroutine AddUniqueElement 
  
  
  
  
  
  
  Subroutine readPWSCFFile(dftFilePath, confFilePath, &
                          dftInEqVol, dftInRadiusCutoff, confWeight)
! Read in configuration file
! readPWSCFFile(dftFilePath, configFilePathT, optEnergyPerAtom, cohEnergy)
    Implicit None  ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 1000000 
    Integer(kind=StandardInteger) :: ios, i, j, k
    Character(len=255) :: dftFilePath, confFilePath
    Character(len=255) :: fileRow, fileLineBuffer
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
      If (ios /= 0) Then
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
      If (ios /= 0) Then
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
!Total energy
        If(fileRow(1:17).eq."!    total energy")Then
          fileLineBuffer = fileRow(33:100)
          read(fileLineBuffer,*) bufferA, bufferB
          read(bufferA,*) configTotalEnergy
          configTotalEnergy = UnitConvert(configTotalEnergy,"RY","EV")
          !configEnergyPerAtom = (configTotalEnergy/(1.0D0*numberOfAtoms))-&
          !  (dftInOptEnergy-dftInCohEnergy)  
        End If      
!Stresses        
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
!Lattice Parameter        
        If(fileRow(1:32).eq."     lattice parameter (alat)  =")Then
          fileLineBuffer = fileRow(33:100)
          read(fileLineBuffer,*) bufferA, bufferB
          read(bufferA,*) aLat
          aLat = UnitConvert(aLat,"Bohr","A")
        End If
!Crystal Axes/Unit Vector - scf
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
!Crystal Axes/Unit Vector - vc-relax
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
!Symbol/Coords - scf
        If(readType.eq.0.and.fileRow(1:22).eq."     site n.     atom ")Then
          Do j=1,numberOfAtoms
            Read(101,"(A255)",IOSTAT=ios) fileRow    
            atomType(j) = TrimSpaces(fileRow(11:23))
            Do k=1,size(dftReplaceLabel,1)
              If(dftReplaceLabel(k,1).eq."        ")Then
                Exit
              End If
              If(Adjustl(StrToUpper(dftReplaceLabel(k,1))).eq.&
                Adjustl(StrToUpper(atomType(j))))Then
                atomType(j) = StrToUpper(dftReplaceLabel(k,2))
                Exit
              End If
            End Do
            fileLineBuffer = fileRow(39:75)
            read(fileLineBuffer,*) bufferA, bufferB, bufferC
            read(bufferA,*) atomCoords(j,1)
            read(bufferB,*) atomCoords(j,2)
            read(bufferC,*) atomCoords(j,3)
          End Do
        End If
!Symbol/Coords - vcrelax
        If(fileRow(1:16).eq."ATOMIC_POSITIONS")Then
          Do j=1,numberOfAtoms
            Read(101,"(A255)",IOSTAT=ios) fileRow    
            atomType(j) = TrimSpaces(fileRow(1:7))
            Do k=1,size(dftReplaceLabel,1)
              If(dftReplaceLabel(k,1).eq."        ")Then
                Exit
              End If
              If(Adjustl(StrToUpper(dftReplaceLabel(k,1))).eq.&
                Adjustl(StrToUpper(atomType(j))))Then
                atomType(j) = StrToUpper(dftReplaceLabel(k,2))
                Exit
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
! Save to file
    open(unit=101,file=trim(confFilePath),status="old",position="append",action="write")
    write(101,"(A27,I8,A9,E20.10)") "#NEW  !added config, atoms ",numberOfAtoms,&
        ", energy ",configTotalEnergy     
    write(101,"(A4,E20.10)") "#LP ",aLat
    write(101,"(A3,F12.7,F12.7,F12.7)") "#X ",crystalAxes(1,1),crystalAxes(1,2),crystalAxes(1,3)
    write(101,"(A3,F12.7,F12.7,F12.7)") "#Y ",crystalAxes(2,1),crystalAxes(2,2),crystalAxes(2,3)
    write(101,"(A3,F12.7,F12.7,F12.7)") "#Z ",crystalAxes(3,1),crystalAxes(3,2),crystalAxes(3,3)
    write(101,"(A4,F12.7,F12.7,F12.7)") "#SX ",stress(1,1),stress(1,2),stress(1,3)
    write(101,"(A4,F12.7,F12.7,F12.7)") "#SY ",stress(2,1),stress(2,2),stress(2,3)
    write(101,"(A4,F12.7,F12.7,F12.7)") "#SZ ",stress(3,1),stress(3,2),stress(3,3)
    write(101,"(A9)") "#CC 1 1 1"
    write(101,"(A4,F10.7)") "#RC ",dftInRadiusCutoff
    write(101,"(A4,F10.7)") "#CW ",confWeight    
    write(101,"(A5,F10.7,A3)") "#EPA ",configEnergyPerAtom," EV"
    If(dftInEqVol.gt.-2.0D20)Then
      write(101,"(A4,F16.7,A5)") "#EV ",dftInEqVol," ANG3" 
    End If
    write(101,"(A4)") "#F Y"
    Do i=1,numberOfAtoms
      write(101,"(A3,A2,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10)") &
      atomType(i),"  ",atomCoords(i,1)," ",atomCoords(i,2)," ",atomCoords(i,3),&
      " ",atomForcess(i,1)," ",atomForcess(i,2)," ",atomForcess(i,3)
    End Do
    write(101,"(A4)") "#END"    
    Close(101)
    
  End Subroutine readPWSCFFile

  
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!   
  Function QueryUniqueElement (element) RESULT (output)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: output 
    Integer(kind=StandardInteger) :: i, k, found
!convert to uppercase
    element = adjustl(StrToUpper(element))   
!loop through elements array	
    k = 0    
    found = 0
    Do i=1,size(elements,1)
      k = k + 1
      If(elements(i).eq.element)Then
        found = 1
        exit
      End If
    End Do
!save element if not found
    If(found.eq.1)Then
      output = k
    Else
      output = 0
    End If
  End function QueryUniqueElement  
  
  Function QueryFunctionType (functionType) RESULT (output)
    Character(len=4) :: functionType
    Integer(kind=StandardInteger) :: output 
    Integer(kind=StandardInteger) :: i
!convert to uppercase
    functionType = adjustl(StrToUpper(functionType))   
    Do i=1,size(eamFunctionTypes,1)
      If(functionType.eq.eamFunctionTypes(i))Then
        output = i
        Exit
      End If
    End Do    
  End function QueryFunctionType  
  
  
End Module readConfig  