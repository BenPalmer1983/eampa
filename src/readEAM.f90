Module readEAM

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
! Privacy of variables/functions/subroutines
  Private    
! Public Subroutines
  Public :: readEAMFile
  Public :: eamDerivatives
  Public :: setEamNodes
  Public :: setEamSpline
  Public :: eamZblHardCore
  Public :: AddUniqueElement
! Public functions 
  Public :: QueryUniqueElement
  
Contains
  Subroutine readEAMFile()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Real(kind=DoubleReal) :: timeStartEAM, timeEndEAM
! Start Time
    Call cpu_time(timeStartEAM)   
! EAM file type
    If(eamFileType.gt.1)Then
      Call eamConvertFile()
    End If
! Synchronise Processes
    Call M_synchProcesses()
! Prepare the eam file
    If(mpiProcessID.eq.0)Then
      Call prepFile()    
    End If
! Send temp file name to all processes
    Call M_distChar(eamFilePathT)
! Read the temporary eam file
    Call readFile()     
! Calculate y'(x) and y''(x)
    Call eamDerivatives()
! Store input potential functions
    eamKeyInput = eamKey
    eamDataInput = eamData
! Set eam spline nodes    
    If(eamForceSpline.eq.1)Then
      Call setEamNodes()
      Call setEamSpline()
    End If    
! Force ZBL if required
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      Call eamZblHardCore()
    End If
! Save the eam file in the output dir if required    
    Call saveEamFile(eamSaveFile) 
! Output summary of EAM to the output file
    Call outputSummary()    
! Synchronise Processes
    Call M_synchProcesses()
! End Time
    Call cpu_time(timeEndEAM)
! Store Time    
    Call storeTime(4,timeEndEAM-timeStartEAM)  
  End Subroutine readEAMFile 
  
  
  Subroutine prepFile()
! Creates a temporary file to contain the EAM potential
! Converts alpha to upper, and strips out comment lines etc
! Will also convert LAMMPS and DLPOLY style EAM
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i
    Character(len=255) :: fileRow
    Character(len=8) :: tempName
! Set temp file path 
    Call tempFileName(tempName)
    eamFilePathT = Trim(tempDirectory)//"/"//tempName//".temp.pot"
    Call fileToClean(eamFilePathT)
! Output to Terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"Reading user eam file ",trim(eamFilePath)
    End If  
! Read EAM file and make new temp EAM pot file
    Open(UNIT=1,FILE=Trim(eamFilePath)) 
    Open(UNIT=2,FILE=Trim(eamFilePathT)) 
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
      fileRow = RemoveComments(fileRow)
      fileRow = RemoveQuotes(fileRow)
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If 
      fileRow = Trim(Adjustl(fileRow))
      If(fileRow(1:1).ne."!".or.fileRow(1:1).ne." ".or.fileRow(1:1).ne."#")Then
        write(2,"(A)") StrToUpper(Trim(fileRow))
        If(fileRow(1:4).eq."DEND".or.fileRow(1:4).eq."EMBD")Then
          eamType = 2
        End If
      End If
    End Do
    Close(2)
    Close(1)
  End Subroutine prepFile 
  
  
  Subroutine readFile()
! Reads in the eam potential from the ew temporary EAM file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, functionCounter, i
    Integer(kind=StandardInteger) :: eamFunctionKey, pointCounter, functionStart, functionLength
    Integer(kind=StandardInteger) :: elementIdA, elementIdB, elementIdMin, elementIdMax
    Character(len=255) :: fileRow
    Character(len=4) :: functionType
    Character(len=2) :: elementA, elementB
    Real(kind=DoubleReal) :: xVal, yVal
! Initialise variables
    eamFunctionKey = 0
    pointCounter = 0
    functionCounter = 0
    functionStart = 1
    functionLength = 0
! Output to Terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"Reading formatted eam file ",trim(eamFilePathT)
    End If  
!---------------------------
! Step 1 - store elements   
!---------------------------
! Read eam pot file and store unique elements, determine EAM type
    Open(UNIT=1,FILE=eamFilePathT) ! Open file
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If 
! Check If New Function
      If(fileRow(1:4).eq."PAIR")Then 
! Store this function elements and type
        Read(fileRow,*) functionType, elementA, elementB  ! Read in line
! Add elements to the config/eam element list
        Call AddUniqueElement(elementA)
        Call AddUniqueElement(elementB)       
      Else If(fileRow(1:4).eq."DENS".or.fileRow(1:4).eq."EMBE".or.&
      fileRow(1:4).eq."DDEN".or.fileRow(1:4).eq."SDEN"&
      .or.fileRow(1:4).eq."DEMB".or.fileRow(1:4).eq."SEMB")Then  
! Store this function elements and type
        Read(fileRow,*) functionType, elementA            ! Read in line
        Call AddUniqueElement(elementA)      
      End If
    End Do
    Close(1)
!---------------------------
! Step 2 - prep function counts
!---------------------------    
    elementsCount = 0
    Do i=1,size(elements)
      If(elements(i).ne."ZZ")Then
        elementsCount = elementsCount + 1
      End If
    End Do
    If(eamType.eq.1)Then        ! Standard EAM
      eamFunctionCount = (elementsCount * ( elementsCount + 5)) / 2
      eamPairCount = (elementsCount * ( elementsCount + 1)) / 2
      eamDensCount = elementsCount
      eamEmbeCount = elementsCount
    End If
    If(eamType.eq.2)Then        ! Two-Band EAM
      If(elementsCount.eq.1)Then
        eamFunctionCount = 1 + (elementsCount * ( elementsCount + 7)) / 2
        eamPairCount = (elementsCount * ( elementsCount + 1)) / 2
        eamSdenCount = 1
        eamDdenCount = elementsCount
        eamSembCount = elementsCount
        eamDembCount = elementsCount
      Else
        eamFunctionCount = (elementsCount * ( elementsCount + 3))
        eamPairCount = (elementsCount * ( elementsCount + 1)) / 2
        eamSdenCount = (elementsCount * ( elementsCount - 1)) / 2
        eamDdenCount = elementsCount
        eamSembCount = elementsCount
        eamDembCount = elementsCount
      End If
    End If  
!---------------------------
! Step 3 - store potential functions   
!---------------------------
! Read eam pot file and store unique elements, determine EAM type
    Open(UNIT=1,FILE=eamFilePathT) ! Open file
    Do i=1,maxFileRows 
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If (ios /= 0) Then
        EXIT 
      End If 
! Check If New Function
!----------------
! PAIR + SDEN
      If(fileRow(1:4).eq."PAIR".or.fileRow(1:4).eq."SDEN")Then 
! Store last function start/length        
        If(functionCounter.gt.0)Then  
          eamKey(eamFunctionKey,4) = functionStart
          eamKey(eamFunctionKey,5) = functionLength
          eamKey(eamFunctionKey,6) = (functionStart + functionLength - 1)
          functionStart = functionStart + functionLength
          functionLength = 0
        End If
! Read function information
        Read(fileRow,*) functionType, elementA, elementB  ! Read in line
        functionCounter = functionCounter + 1             ! Increment function counter
! Set function key
        elementIdA = QueryUniqueElement(elementA)
        elementIdB = QueryUniqueElement(elementB)
        elementIdMin = min(elementIdA,elementIdB)-1
        elementIdMax = max(elementIdA,elementIdB)-1
        If(fileRow(1:4).eq."PAIR")Then
          eamFunctionKey = 1+elementIdMin+(elementIdMax*(elementIdMax+1))/2       
        End If          
        If(fileRow(1:4).eq."SDEN")Then
          If(elementsCount.eq.1)Then
            eamFunctionKey = eamPairCount + 1
          Else
            eamFunctionKey = eamPairCount+1+elementIdMin+(elementIdMax*(elementIdMax-1))/2
          End If 
        End If   
! Store this function elements and type        
        eamKey(eamFunctionKey,1) = QueryUniqueElement(elementA)  
        eamKey(eamFunctionKey,2) = QueryUniqueElement(elementB)     
        eamKey(eamFunctionKey,3) = QueryFunctionType(functionType)  
!----------------
! DENS, EMBE, DDEN, DEMB, SEMB      
      Else If(fileRow(1:4).eq."DENS".or.fileRow(1:4).eq."EMBE".or.&
      fileRow(1:4).eq."DDEN"&
      .or.fileRow(1:4).eq."DEMB".or.fileRow(1:4).eq."SEMB")Then
! Store last function start/length
        If(functionCounter.gt.0)Then  
          eamKey(eamFunctionKey,4) = functionStart
          eamKey(eamFunctionKey,5) = functionLength
          eamKey(eamFunctionKey,6) = (functionStart + functionLength - 1)
          functionStart = functionStart + functionLength
          functionLength = 0
        End If      
! Store this function elements and type
        Read(fileRow,*) functionType, elementA            ! Read in line
        functionCounter = functionCounter + 1             ! Increment function counter
! Set function key
        elementIdA = QueryUniqueElement(elementA)
        If(fileRow(1:4).eq."DENS")Then
          eamFunctionKey = eamPairCount+elementIdA            
        End If
        If(fileRow(1:4).eq."EMBE")Then
          eamFunctionKey = eamPairCount+eamDensCount+elementIdA            
        End If
        If(fileRow(1:4).eq."DDEN")Then
          eamFunctionKey = eamPairCount+eamSdenCount+elementIdA            
        End If
        If(fileRow(1:4).eq."SEMB")Then
          eamFunctionKey = eamPairCount+eamSdenCount+eamDdenCount+elementIdA            
        End If
        If(fileRow(1:4).eq."DEMB")Then
          eamFunctionKey = eamPairCount+eamSdenCount+eamDdenCount+eamSembCount+elementIdA            
        End If
        eamKey(eamFunctionKey,1) = QueryUniqueElement(elementA)  
        eamKey(eamFunctionKey,2) = 0  
        eamKey(eamFunctionKey,3) = QueryFunctionType(functionType) 
      Else
! Read in function/al x and y
        pointCounter = pointCounter + 1
        functionLength = functionLength + 1
        Read(fileRow,*) xVal, yVal        
        eamData(pointCounter,1) = xVal
        eamData(pointCounter,2) = yVal
      End If
    End Do
! Store last function    
    If(functionCounter.gt.0)Then  !Store last function start/length
      eamKey(eamFunctionKey,4) = functionStart
      eamKey(eamFunctionKey,5) = functionLength
      eamKey(eamFunctionKey,6) = (functionStart + functionLength - 1)
    End If
    Close(1) ! Close file
  End Subroutine readFile 
  
  
  Subroutine eamConvertFile()
! Fills in first and second order derivatives
    Implicit None   ! Force declaration of all variables  
! Private variables    
    Character(len=255) :: eamFilePathC 
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i, j, k, n, points, rowPoints  
    Integer(kind=StandardInteger) :: rowCount
    Character(len=255) :: fileRow
    Character(len=128) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
    Integer(kind=StandardInteger) :: elementTypeCount, elementZ
    Integer(kind=StandardInteger) :: nrho, nr, potentialLinesRho, potentialLinesR
    Real(kind=DoubleReal) :: drho, dr, cutoff
    Real(kind=DoubleReal) :: radius, rho, tempDouble
    Real(kind=DoubleReal) :: rStart, rEnd, rInc
    Real(kind=DoubleReal) :: x, y     
    Character(len=2), Dimension(1:300) :: functionElements   
! Init variables
    functionElements = "ZZ"    
! Set converted eam potential file name
    eamFilePathC = Trim(tempDirectory)//"/"//Trim(eamFilePath)//".pot"
    Call fileToClean(eamFilePathC)
!--------------------  
! LAMMPS
!--------------------    
    If(eamFileType.eq.2.and.mpiProcessID.eq.0)Then    ! LAMMPS file type
!open output potential file (write to)
      Open(UNIT=2,FILE=trim(eamFilePathC)) 
!open LAMMPS input potential file (read from)
      Open(UNIT=1,FILE=trim(eamFilePath)) 
      rowCount = 0
      Do i=1,maxFileRows 
!Read in line
        Read(1,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
        If (ios /= 0) then
          EXIT 
        End If        
        If(fileRow(1:1).ne."#")Then          
          rowCount = rowCount + 1
          If(rowCount.eq.1)Then
            Read(fileRow,*) elementTypeCount    
          End If
          If(rowCount.eq.2)Then
            Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
            Read(bufferA,*) nrho    !number of points at which p(r) is evaluated
            Read(bufferB,*) drho    
            Read(bufferC,*) nr        !number of points at which V(r) and F(p) are evaluated
            Read(bufferD,*) dr    
            Read(bufferE,*) cutoff            
            potentialLinesR = Ceiling(1.0D0*(nrho/5))
            potentialLinesRho = Ceiling((1.0D0*nrho)/5.0D0)
          End If
          If(rowCount.eq.3)Then
!Read in embedding function and density function for each element
            Do j=1,elementTypeCount
              If(j.gt.1)Then
                !Read in next element row
                Read(1,"(A255)",IOSTAT=ios) fileRow    
              End If
!Read element type                
              Read(fileRow,*) elementZ
!Embedding function
              write(2,"(A8,A2)") "EMBE    ",elementSymbol(elementZ)
              functionElements(j) = elementSymbol(elementZ)
              rho = 0.0D0
!Embedding function
              Do k=1,potentialLinesRho
                Read(1,"(A255)",IOSTAT=ios) fileRow
                If(k.lt.potentialLinesRho.or.mod(nrho,5).eq.0)Then
                  Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE     
                  Read(bufferA,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho 
                  Read(bufferB,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferC,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferD,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferE,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                Else  
                  If(mod(nrho,5).eq.1)Then
                    Read(fileRow,*) bufferA     
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If 
                  If(mod(nrho,5).eq.2)Then
                    Read(fileRow,*) bufferA, bufferB
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If
                  If(mod(nrho,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho 
                    Read(bufferC,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If
                  If(mod(nrho,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho 
                    Read(bufferC,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho 
                    Read(bufferD,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If                  
                End If
              End Do
!Density function
              write(2,"(A8,A2)") "DENS    ",elementSymbol(elementZ)
              radius = 0.0D0  
!Density function
              Do k=1,potentialLinesR+1
                Read(1,"(A255)",IOSTAT=ios) fileRow
                If(k.lt.(potentialLinesR+1).or.(k.eq.(potentialLinesR+1).and.mod(nr,5).eq.0))Then
                  Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE     
                  Read(bufferA,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr 
                  Read(bufferB,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferC,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferD,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferE,*) tempDouble 
                  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                Else  
                  If(mod(nr,5).eq.1)Then
                    Read(fileRow,*) bufferA     
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If 
                  If(mod(nr,5).eq.2)Then
                    Read(fileRow,*) bufferA, bufferB
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If
                  If(mod(nr,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr 
                    Read(bufferC,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If
                  If(mod(nr,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                    Read(bufferA,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr 
                    Read(bufferC,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr 
                    Read(bufferD,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If                  
                End If
              End Do
            End Do !Element count            
!Pair potentials
            Do j=1,elementTypeCount
              Do n=j,elementTypeCount
                write(2,"(A8,A2,A4,A2)") "PAIR    ",&
                functionElements(j),"    ",functionElements(n)
                radius = 0.0D0
!Pair function
                Do k=1,potentialLinesR+1
                  Read(1,"(A255)",IOSTAT=ios) fileRow
                  If(k.lt.(potentialLinesR+1).or.&
                  (k.eq.(potentialLinesR+1).and.mod(nr,5).eq.0))Then
                    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE     
                    Read(bufferA,*) tempDouble 
                    If(radius.eq.0.0D0)Then
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",0.0D0
                    Else
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    radius = radius + dr 
                    Read(bufferB,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                    Read(bufferC,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                    Read(bufferD,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                    Read(bufferE,*) tempDouble 
                    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                  Else  
                    If(mod(nr,5).eq.1)Then
                      Read(fileRow,*) bufferA     
                      Read(bufferA,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If 
                    If(mod(nr,5).eq.2)Then
                      Read(fileRow,*) bufferA, bufferB
                      Read(bufferA,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr 
                      Read(bufferB,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    If(mod(nr,5).eq.3)Then
                      Read(fileRow,*) bufferA, bufferB, bufferC
                      Read(bufferA,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr 
                      Read(bufferB,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr 
                      Read(bufferC,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    If(mod(nr,5).eq.3)Then
                      Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                      Read(bufferA,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr 
                      Read(bufferB,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr 
                      Read(bufferC,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr 
                      Read(bufferD,*) tempDouble 
                      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If                  
                  End If
                End Do
              End Do
            End Do
          End If
        End If        
      End Do
      Close(1)
      Close(2)
    End If
!--------------------  
! DLPOLY
!--------------------    
    If(eamFileType.eq.3.and.mpiProcessID.eq.0)Then    ! DLPOLY file type
!open output potential file (write to)
      Open(UNIT=2,FILE=trim(eamFilePathC)) 
!open DLPOLY input potential file (read from)
      Open(UNIT=1,FILE=trim(eamFilePath)) 
      Do i=1,maxFileRows 
!Read in line
        Read(1,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
        If (ios /= 0) then
          EXIT 
        End If        
        If(fileRow(1:1).ne."#")Then   
! PAIR function        
          If(StrToUpper(fileRow(1:4)).eq."PAIR".or.&
          StrToUpper(fileRow(1:4)).eq."DENS".or.&
          StrToUpper(fileRow(1:4)).eq."EMBE")Then  
            If(StrToUpper(fileRow(1:4)).eq."PAIR")Then             
              Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
              Read(bufferD,*) points
              n = Ceiling(points / 4.0D0)
              Read(bufferE,*) rStart 
              Read(bufferF,*) rEnd
              rInc = (rEnd-rStart)/(points-1)   
              x = rStart            
              write(2,"(A8,A2,A4,A2)") "PAIR    ",bufferB,"    ",bufferC
            End If  
            If(StrToUpper(fileRow(1:4)).eq."DENS")Then             
              Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
              Read(bufferC,*) points
              n = Ceiling(points / 4.0D0)
              Read(bufferD,*) rStart 
              Read(bufferE,*) rEnd
              rInc = (rEnd-rStart)/(points-1)   
              x = rStart            
              write(2,"(A8,A2)") "DENS    ",bufferB
            End If 
            If(StrToUpper(fileRow(1:4)).eq."EMBE")Then             
              Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
              Read(bufferC,*) points
              n = Ceiling(points / 4.0D0)
              Read(bufferD,*) rStart 
              Read(bufferE,*) rEnd
              rInc = (rEnd-rStart)/(points-1)   
              x = rStart            
              write(2,"(A8,A2)") "EMBE    ",bufferB
            End If  
            Do j=1,n              
              Read(1,"(A255)",IOSTAT=ios) fileRow
              If(j.lt.n)Then
                Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                Read(bufferA,*) y
                write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                x = x + rInc
                Read(bufferB,*) y
                write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                x = x + rInc
                Read(bufferC,*) y
                write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                x = x + rInc
                Read(bufferD,*) y
                write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                x = x + rInc                
              Else
                rowPoints = points-((n-1)*4)
                If(rowPoints.eq.4)Then
                  Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                  Read(bufferA,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                  x = x + rInc
                  Read(bufferB,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                  x = x + rInc
                  Read(bufferC,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                  x = x + rInc
                  Read(bufferD,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                End If
                If(rowPoints.eq.3)Then
                  Read(fileRow,*) bufferA, bufferB, bufferC
                  Read(bufferA,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                  x = x + rInc
                  Read(bufferB,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                  x = x + rInc
                  Read(bufferC,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                End If
                If(rowPoints.eq.2)Then
                  Read(fileRow,*) bufferA, bufferB
                  Read(bufferA,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                  x = x + rInc
                  Read(bufferB,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                End If
                If(rowPoints.eq.1)Then
                  Read(fileRow,*) bufferA
                  Read(bufferA,*) y
                  write(2,"(E24.16E3,A2,E24.16E3)") x,"  ",y
                End If
              End If
            End Do
          End If
        End If
      End Do
      Close(1)
      Close(2)
    End If
! Update the EAM input file
    eamFilePath = eamFilePathC
  End Subroutine eamConvertFile 
  
  
  
  Subroutine eamDerivatives()
! Fills in first and second order derivatives
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j, k, xPoint, pointOffset
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: pointsInterp
! Calculate y'(x) and y''(x)
    functionCounter = 0
    Do i=1,size(eamKey,1)
      If(eamKey(i,1).gt.0)Then
        functionCounter = functionCounter + 1
        Do j=eamKey(i,4),eamKey(i,6)
          pointOffset = -2
          If((j+pointOffset).lt.eamKey(i,4))Then
            pointOffset = 0
          End If
          If((j+pointOffset+3).gt.eamKey(i,6))Then
            pointOffset = -3
          End If        
          Do k=0,3
            xPoint = j + k + pointOffset
            pointsInterp(k+1,1) = eamData(xPoint,1)
            pointsInterp(k+1,2) = eamData(xPoint,2)
          End Do          
          yArray = PointInterp(pointsInterp,eamData(j,1),4,3)
          eamData(j,3) = yArray(2)
          eamData(j,4) = yArray(3)
        End Do
      End If
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do
  End Subroutine eamDerivatives   
  
  
  Subroutine setEamNodes()
! Set nodes
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: functionCounter, i, j 
    Integer(kind=StandardInteger) :: eamStart, eamLength, eamEnd
    Integer(kind=StandardInteger) :: nodes, nodeKey, functionType, fixEndNode
    Real(kind=DoubleReal) :: x, xStart, xEnd
    Real(kind=DoubleReal), Dimension(1:3) :: yArray    
! Init variables    
    splineNodesKey = -1         ! reset key array
    splineNodesData = 0.0D0     ! reset node array
    functionCounter = 0
    nodeKey = 0
! Loop through EAM functions
    Do i=1,size(eamKey,1)
      If(eamKey(i,1).gt.0)Then        
        functionCounter = functionCounter + 1
        nodes = splineNodeCount(eamKey(i,3))
        functionType = eamKey(i,3)
        eamStart = eamKey(i,4)
        eamLength = eamKey(i,5)
        eamEnd = eamKey(i,6)
        fixEndNode = 0
        If(functionType.eq.1.or.functionType.eq.2.or.&
        functionType.eq.4.or.functionType.eq.5)Then
          fixEndNode = 1
        End If
        Do j=1,nodes
          nodeKey = nodeKey + 1
          If(j.eq.1)Then
            splineNodesData(nodeKey,1) = eamData(eamStart,1)
            splineNodesData(nodeKey,2) = eamData(eamStart,2)
            splineNodesData(nodeKey,3) = eamData(eamStart,3)
            splineNodesData(nodeKey,4) = eamData(eamStart,4)
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0     
            splineNodesKey(i,1) = eamKey(i,1)
            splineNodesKey(i,2) = eamKey(i,2)
            splineNodesKey(i,3) = eamKey(i,3)
            splineNodesKey(i,4) = nodeKey               ! Start
            splineNodesKey(i,5) = nodes                 ! Length
            splineNodesKey(i,6) = nodeKey + nodes - 1   ! End            
          Else If(j.eq.nodes)Then
            splineNodesData(nodeKey,1) = eamData(eamEnd,1)
            splineNodesData(nodeKey,2) = eamData(eamEnd,2)
            splineNodesData(nodeKey,3) = eamData(eamEnd,3)
            splineNodesData(nodeKey,4) = eamData(eamEnd,4)
            splineNodesData(nodeKey,5) = 1.0D0*j            ! Function node counter
            If(fixEndNode.eq.1)Then
              splineNodesData(nodeKey,6) = 1.0D0
            Else              
              splineNodesData(nodeKey,6) = 0.0D0            
            End If
          Else
            xStart = eamData(eamStart,1)
            xEnd = eamData(eamEnd,1)
            x = xStart+1.0D0*(j-1)*((xEnd-xStart)/(nodes-1))
            yArray = PointInterp(eamData,x,eamInterpPoints,1,eamStart,eamLength)
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = yArray(1)
            splineNodesData(nodeKey,3) = yArray(2)
            splineNodesData(nodeKey,4) = yArray(3)
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0     
          End If          
! Fix nodes below ZBL spline/core
          If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
            If(functionType.eq.1.and.splineNodesData(nodeKey,1).le.zblHardCore(2))Then
              splineNodesData(nodeKey,6) = 1.0D0
            End If
          End If          
          !If(mpiProcessID.eq.0)Then
          !print *,nodeKey,splineNodesData(nodeKey,1),splineNodesData(nodeKey,2),&
          !splineNodesData(nodeKey,3),splineNodesData(nodeKey,4),&
          !splineNodesData(nodeKey,5),splineNodesData(nodeKey,6)
          !End If
        End Do
      End If  
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do  
! Store total number of nodes   
    splineTotalNodes = nodeKey
! Output nodes to file    
    If(eamNodesFilePath(1:1).ne." ")Then
      Call outputSplineNodes(eamNodesFilePath)
    End If
  End Subroutine setEamNodes  
  
  
  
  Subroutine setEamSpline()
! Force ZBL Core
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
    Integer(kind=StandardInteger) :: nodes, nodeStart, nodeLength, nodeEnd, eamPoint
    Real(kind=DoubleReal), Dimension(1:1001,1:4) :: splineDataPoints
! Init variables    
    splineDataPoints = 0.0D0
! Loop through EAM functions
    functionCounter = 0
    Do i=1,size(eamKey,1)
      If(eamKey(i,1).gt.0)Then        
        functionCounter = functionCounter + 1        
        nodes = splineNodeCount(eamKey(i,3))
! Set node start/end points
        nodeStart = splineNodesKey(i,4)
        nodeLength = splineNodesKey(i,5)
        nodeEnd = splineNodesKey(i,6)
        splineDataPoints = SplineNodes(splineNodesData,1001,nodeStart,nodeEnd)        
        Do j=1,1001
          eamPoint = eamKey(i,4)+j-1
          eamData(eamPoint,1) = splineDataPoints(j,1)
          eamData(eamPoint,2) = splineDataPoints(j,2)
          eamData(eamPoint,3) = splineDataPoints(j,3)
          eamData(eamPoint,4) = splineDataPoints(j,4) 
        End Do
      End If  
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do        
  End Subroutine setEamSpline  
    
  
  Subroutine setEamSplineOpt()
! Force ZBL Core
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
    Integer(kind=StandardInteger) :: nodes, nodeKey, nodeStart, nodeLength, nodeEnd, eamPoint
    Real(kind=DoubleReal), Dimension(1:1001,1:4) :: splineDataPointsOpt
! Init variables    
    splineDataPointsOpt = 0.0D0
! Loop through EAM functions
    functionCounter = 0
    nodeKey = 0
    Do i=1,size(eamKeyOpt,1)
      If(eamKeyOpt(i,1).gt.0)Then        
        functionCounter = functionCounter + 1        
        nodes = splineNodeCount(eamKeyOpt(i,3))
! Set node start/end points
        nodeStart = splineNodesKeyOpt(i,4)
        nodeLength = splineNodesKeyOpt(i,5)
        nodeEnd = splineNodesKeyOpt(i,6)
        splineDataPointsOpt = SplineNodes(splineNodesDataOpt,1001,nodeStart,nodeEnd)        
        Do j=1,1001
          eamPoint = eamKeyOpt(i,4)+j-1
          eamDataOpt(eamPoint,1) = splineDataPointsOpt(j,1)
          eamDataOpt(eamPoint,2) = splineDataPointsOpt(j,2)
          eamDataOpt(eamPoint,3) = splineDataPointsOpt(j,3)
          eamDataOpt(eamPoint,4) = splineDataPointsOpt(j,4) 
        End Do
      End If  
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do        
  End Subroutine setEamSplineOpt  
  
  
  
  Subroutine eamZblHardCore()
! Force ZBL Core
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j, zA, zB
    Real(kind=DoubleReal) :: xA, xB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:6) :: splineCoeffs    
! Loop through EAM functions
    If(zblHardCore(1).gt.0.0D0.and.eamForceZBL.eq.1)Then
      functionCounter = 0
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
! Pair potentials
          If(eamKey(i,3).eq.1)Then
! ZBL parameters
            zA = elementsCharge(eamKey(i,1))
            zB = elementsCharge(eamKey(i,2))
            xA = zblHardCore(1)
! Point A
            yArray = ZblFull (xA, zA, zB)
            pointA(1) = xA
            pointA(2) = yArray(1)
            pointA(3) = yArray(2)
            pointA(4) = yArray(3)
! Point B
            xB = zblHardCore(2)          
            yArray = PointInterp(eamData,xB,eamInterpPoints,2,eamKey(i,4),eamKey(i,5)) 
            pointB(1) = xB          
            pointB(2) = yArray(1)
            pointB(3) = yArray(2)
            pointB(4) = yArray(3)
! Get spline coefficients
            splineCoeffs = SplineAB(pointA, pointB)
! Alter eamData
            Do j=eamKey(i,4),eamKey(i,6)
              If(eamData(j,1).lt.pointA(1))Then   ! ZBL Core
                xA = eamData(j,1)
                yArray = ZblFull (xA, zA, zB)
                eamData(j,2) = yArray(1)
                eamData(j,3) = yArray(2)
                eamData(j,4) = yArray(3)
              End If
              If(eamData(j,1).ge.pointA(1).and.eamData(j,1).le.pointB(1))Then   ! Spline
                xA = eamData(j,1)
                eamData(j,2) = CalcPolynomial(splineCoeffs, xA, 0)
                eamData(j,3) = CalcPolynomial(splineCoeffs, xA, 1)
                eamData(j,4) = CalcPolynomial(splineCoeffs, xA, 2)
              End If
            End Do
! Store to file
            Call outputZBL(zA, zB, pointA, pointB, splineCoeffs)
          End If
        End If  
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do  
    End If
  End Subroutine eamZblHardCore   
  
  
  
  

 
  
  Subroutine outputSummary()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, functionCounter
! Print out
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Read EAM Potential Functions"
      print *,"EAM Type: ",eamType
      If(eamType.eq.1)Then
        print *,"Pair: ",eamPairCount,", Dens: ",eamDensCount,&
        ", Dens: ",eamEmbeCount,", Total: ",eamFunctionCount
      End If
      print *,"EAM Potential Functions Summary"
      functionCounter = 0
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
          If(eamKey(i,2).gt.0)Then
            print "(I4,A1,A4,A1,I2,A2,A2,A1,I2,A2,A2,A1,I2,A2,I7,A1,I7,A1,I7)",&
            functionCounter," ",eamFunctionTypes(eamKey(i,3)),"(",eamKey(i,3),") ",&
            elements(eamKey(i,1)),"(",eamKey(i,1),") ",&
            elements(eamKey(i,2)),"(",eamKey(i,3),") ",&
            eamKey(i,4)," ",eamKey(i,5)," ",eamKey(i,6)
          Else
            print "(I4,A1,A4,A1,I2,A2,A2,A1,I2,A2,A7,I7,A1,I7,A1,I7)",&
            functionCounter," ",eamFunctionTypes(eamKey(i,3)),"(",eamKey(i,3),") ",&
            elements(eamKey(i,1)),"(",eamKey(i,1),") ",&
            "       ",&
            eamKey(i,4)," ",eamKey(i,5)," ",eamKey(i,6)
          End If
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
    End If 
  End Subroutine outputSummary 
  
  
  
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
  
  
  
End Module readEAM  