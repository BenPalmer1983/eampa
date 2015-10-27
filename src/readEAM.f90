Module readEAM

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
  Use types
  Use plotTypes
  Use plot
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use initialise
  Use loadData
  Use globals
  Use output
  Use eamGen
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: readEAMFile
  Public :: eamDerivatives
  Public :: setEamNodes
  Public :: rescaleEmbedNodes
  Public :: countEamNodes
  Public :: setEamSpline
  Public :: eamZblHardCore
  Public :: eamPairZbl
  Public :: AddUniqueElement
  Public :: loadInputEAM
  Public :: eamCharts
! Public functions
  Public :: QueryUniqueElement

  
  Contains
  Subroutine readEAMFile()
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: timeStartEAM, timeEndEAM
! Start Time
    Call cpu_time(timeStartEAM)
! Prepare the eam file
    Call prepFile()
! Synchronise Processes
    Call M_synchProcesses()
! EAM file type
    If(eamFileType.gt.1)Then
      Call eamConvertFile()
    End If
! Synchronise Processes
    Call M_synchProcesses()
! Read the temporary eam file
    Call readFileE()
! Calculate y'(x) and y''(x)
    Call eamDerivatives()
! Save the prepared eam functions in the output dir
    Call saveEamFile("input_prepared.pot")    
! Store input potential functions
    eamKeyInput = eamKey
    eamDataInput = eamData
! Set eam spline nodes
    If(eamForceSpline)Then
      Call setEamNodes()
      Call saveEamNodes("input_splined.nodes")
      Call setEamSpline()
      Call saveEamFile("input_splined.pot")   
    End If
! Output summary of EAM to the output file
    Call outputSummary()
! Export EAM data to file
    Call eamDataDump()
! Synchronise Processes
    Call M_synchProcesses()
    If(makeEAMCharts)Then
      Call eamCharts()
    End If  
    Call M_synchProcesses()
! End Time
    Call cpu_time(timeEndEAM)
! Store Time
    eamLoadTime = timeEndEAM-timeStartEAM
  End Subroutine readEAMFile

! -----------------------------------------------
! Prepare file and load into memory
! -----------------------------------------------
  Subroutine prepFile()
! Converts alpha to upper, and strips out comment lines etc
! Will also convert LAMMPS and DLPOLY style EAM
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: ios, i, n
    Character(len=255) :: fileRow
! Output to Terminal
    If(TerminalPrint())Then
      Print *,"Reading user eam file ",trim(eamFilePath)
    End If
    Open(UNIT=1,FILE=Trim(eamFilePath))
    n = 0
    eamType = 1  ! Set EAM as default
    Do i=1,maxFileRows
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
      fileRow = RemoveComments(fileRow)
      fileRow = RemoveQuotes(fileRow)
! Break out If end of file
      If(ios /= 0)Then
        EXIT
      End If
      fileRow = Trim(Adjustl(fileRow))
! Read row into
      If(fileRow(1:1).eq."!".or.fileRow(1:1).eq." ".or.fileRow(1:1).eq."#")Then
! Check for EAM type
        If(fileRow(1:1).eq."#")Then
          fileRow = StrToUpper(fileRow)
          If(fileRow(1:4).eq."#EAM")Then
            eamType = 1  ! Standard EAM
          End If
          If(fileRow(1:6).eq."#TBEAM")Then
            eamType = 2  ! Two Band EAM
          End If
          If(fileRow(1:6).eq."#2BEAM")Then
            eamType = 2  ! Two Band EAM
          End If
          If(fileRow(1:6).eq."#3BEAM")Then
            eamType = 3  ! Two Band EAM
          End If
        End If
      Else
        n = n + 1
        eamInputData(n) = fileRow
        If(fileRow(1:4).eq."DEND".or.fileRow(1:4).eq."EMBD")Then
          eamType = 2
        End If
      End If
    End Do
    Close(1)
    If(TerminalPrint())Then
      Print *,"EAM loaded into memory.  Type: ",eamType
    End If
  End Subroutine prepFile

! -----------------------------------------------
! Convert file type
! -----------------------------------------------
  Subroutine eamConvertFile()
! Fills in first and second order derivatives
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(len=255) :: eamFilePathC
    Integer(kind=StandardInteger) :: ios, i, j, k, n, m, w, points, rowPoints
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
    Character(len=255), Dimension(1:65536) :: eamInputDataT
    eamInputDataT = BlankStringArray(eamInputDataT)
! Init variables
    FunctionElements = "ZZ"
! Set converted eam potential file name
    eamFilePathC = Trim(tempDirectory)//"/"//Trim(eamFilePath)//".pot"
! --------------------
! LAMMPS
! --------------------
!    If(eamFileType.eq.2.and.mpiProcessID.eq.0)Then    ! LAMMPS file type
    If(eamFileType.eq.2)Then    ! LAMMPS file type
      rowCount = 0
      m = 0
      w = 0
      Do i=1,maxFileRows
        m = m + 1
        fileRow = eamInputData(m)     !read line
        If(fileRow(1:1).eq." ")Then
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
! Read in embedding function and density function for each element
            Do j=1,elementTypeCount
              If(j.gt.1)Then
! Read in next element row
! Read(1,"(A255)",IOSTAT=ios) fileRow
                m = m + 1
                fileRow = eamInputData(m)     ! read next line
              End If
! Read element type
              Read(fileRow,*) elementZ
! Embedding function
! write(2,"(A8,A2)") "EMBE    ",elementSymbol(elementZ)
              w = w + 1
              write(eamInputDataT(w),"(A8,A2)") "EMBE    ",elementSymbol(elementZ)
              FunctionElements(j) = elementSymbol(elementZ)
              rho = 0.0D0
! Embedding function
              Do k=1,potentialLinesRho
! Read(1,"(A255)",IOSTAT=ios) fileRow
                m = m + 1
                fileRow = eamInputData(m)     ! read next line
                If(k.lt.potentialLinesRho.or.mod(nrho,5).eq.0)Then
                  Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
                  Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferD,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                  Read(bufferE,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  rho = rho + drho
                Else
                  If(mod(nrho,5).eq.1)Then
                    Read(fileRow,*) bufferA
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If
                  If(mod(nrho,5).eq.2)Then
                    Read(fileRow,*) bufferA, bufferB
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If
                  If(mod(nrho,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho
                    Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If
                  If(mod(nrho,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho
                    Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    rho = rho + drho
                    Read(bufferD,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
                  End If
                End If
              End Do
! Density function
! write(2,"(A8,A2)") "DENS    ",elementSymbol(elementZ)
              w = w + 1
              write(eamInputDataT(w),"(A8,A2)") "DENS    ",elementSymbol(elementZ)
              radius = 0.0D0
! Density function
              Do k=1,potentialLinesR+1
! Read(1,"(A255)",IOSTAT=ios) fileRow
                m = m + 1
                fileRow = eamInputData(m)     ! read next line
                If(k.lt.(potentialLinesR+1).or.(k.eq.(potentialLinesR+1).and.mod(nr,5).eq.0))Then
                  Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
                  Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferD,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                  Read(bufferE,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  w = w + 1
                  write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  radius = radius + dr
                Else
                  If(mod(nr,5).eq.1)Then
                    Read(fileRow,*) bufferA
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If
                  If(mod(nr,5).eq.2)Then
                    Read(fileRow,*) bufferA, bufferB
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If
                  If(mod(nr,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr
                    Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If
                  If(mod(nr,5).eq.3)Then
                    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                    Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr
                    Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    radius = radius + dr
                    Read(bufferD,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                    w = w + 1
                    write(eamInputDataT(w),"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
                  End If
                End If
              End Do
            End Do !Element count
! Pair potentials
            Do j=1,elementTypeCount
              Do n=j,elementTypeCount
                write(2,"(A8,A2,A4,A2)") "PAIR    ",&
                FunctionElements(j),"    ",functionElements(n)
                radius = 0.0D0
! Pair function
                Do k=1,potentialLinesR+1
! Read(1,"(A255)",IOSTAT=ios) fileRow
                  m = m + 1
                  fileRow = eamInputData(m)     ! read next line
                  If(k.lt.(potentialLinesR+1).or.&
                    (k.eq.(potentialLinesR+1).and.mod(nr,5).eq.0))Then
                    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
                    Read(bufferA,*) tempDouble
                    If(radius.eq.0.0D0)Then
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",0.0D0
                    Else
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    radius = radius + dr
                    Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                    Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                    Read(bufferD,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                    Read(bufferE,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    radius = radius + dr
                  Else
                    If(mod(nr,5).eq.1)Then
                      Read(fileRow,*) bufferA
                      Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    If(mod(nr,5).eq.2)Then
                      Read(fileRow,*) bufferA, bufferB
                      Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr
                      Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    If(mod(nr,5).eq.3)Then
                      Read(fileRow,*) bufferA, bufferB, bufferC
                      Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr
                      Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr
                      Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                    If(mod(nr,5).eq.3)Then
                      Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
                      Read(bufferA,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr
                      Read(bufferB,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr
                      Read(bufferC,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                      radius = radius + dr
                      Read(bufferD,*) tempDouble
! write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
                    End If
                  End If
                End Do
              End Do
            End Do
          End If
        End If
      End Do
! Save file
      eamInputData = eamInputDataT
    End If
! --------------------
! DLPOLY
! --------------------
    If(eamFileType.eq.3.and.mpiProcessID.eq.0)Then    ! DLPOLY file type
      m = 0
      w = 0
      Do i=1,maxFileRows
        m = m + 1
        fileRow = eamInputData(m)     !read line
! Break out
        If(fileRow(1:1).eq." ")Then
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
! Save file
      eamInputData = eamInputDataT
    End If
! Update the EAM input file
    eamFilePath = eamFilePathC
  End Subroutine eamConvertFile

! -----------------------------------------------
! Read file data into eam arrays
! -----------------------------------------------

  Subroutine readFileE()
! Reads in the eam potential from the ew temporary EAM file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter,i,m
    Integer(kind=StandardInteger) :: eamFunctionKey, pointCounter, functionStart, functionLength
    Integer(kind=StandardInteger) :: elementIdA, elementIdB, elementIdMin, elementIdMax
    Character(len=255) :: fileRow
    Character(len=4) :: functionType
    Character(len=2) :: elementA, elementB
    Real(kind=DoubleReal) :: xVal, yVal
! Initialise variables
    eamFunctionKey = 0
    pointCounter = 0
    FunctionCounter = 0
    FunctionStart = 1
    FunctionLength = 0
! Output to Terminal
    If(TerminalPrint())Then
      Print *,"Reading data to eam arrays"
    End If
! ---------------------------
! Step 1 - store elements
! ---------------------------
    m = 0
    Do i=1,maxFileRows
      m = m + 1
      fileRow = eamInputData(m)     !read line
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
! Check If New Function
      If(fileRow(1:4).eq."PAIR")Then
! Store this function elements and type
        Read(fileRow,*) functionType, elementA, elementB  ! Read in line
! Add elements to the config/eam element list
        Call AddUniqueElement(elementA)
        Call AddUniqueElement(elementB)
      ElseIf(fileRow(1:4).eq."DENS".or.fileRow(1:4).eq."EMBE".or.&
        fileRow(1:4).eq."DDEN".or.fileRow(1:4).eq."SDEN"&
        .or.fileRow(1:4).eq."DEMB".or.fileRow(1:4).eq."SEMB")Then
! Store this function elements and type
        Read(fileRow,*) functionType, elementA            ! Read in line
        Call AddUniqueElement(elementA)
      End If
    End Do
! Output to Terminal
    If(TerminalPrint())Then
      Print *,"Elements loaded from potential:"
      Do i=1,size(elements,1)
        If(elements(i).eq."ZZ")Then
          exit
        End If
        Print *,elements(i)
      End Do
    End If
! ---------------------------
! Step 2 - prep function counts
! ---------------------------
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
! Output to Terminal
    If(TerminalPrint())Then
      Print *,"EAM elements: ",elementsCount
      Print *,"EAM type: ",eamType
      Print *,"Expected functions: ",eamFunctionCount
    End If
! ---------------------------
! Step 3 - store potential functions
! ---------------------------
    m = 0
    Do i=1,maxFileRows
      m = m + 1
      fileRow = eamInputData(m)     !read line
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
! Check If New Function
! ----------------
! PAIR + SDEN
      If(fileRow(1:4).eq."PAIR".or.fileRow(1:4).eq."SDEN")Then
! Store last function start/length
        If(functionCounter.gt.0)Then
          eamKey(eamFunctionKey,4) = functionStart
          eamKey(eamFunctionKey,5) = functionLength
          eamKey(eamFunctionKey,6) = (functionStart + functionLength - 1)
          FunctionStart = functionStart + functionLength
          FunctionLength = 0
        End If
! Read function information
        Read(fileRow,*) functionType, elementA, elementB  ! Read in line
        FunctionCounter = functionCounter + 1             ! Increment function counter
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
        eamKey(eamFunctionKey,1) = QueryUniqueElement(elementA)    ! Element A id
        eamKey(eamFunctionKey,2) = QueryUniqueElement(elementB)    ! Element B id
        eamKey(eamFunctionKey,3) = QueryFunctionType(functionType) ! Type of EAM function
! ----------------
! DENS, EMBE, DDEN, DEMB, SEMB
      ElseIf(fileRow(1:4).eq."DENS".or.fileRow(1:4).eq."EMBE".or.&
        fileRow(1:4).eq."DDEN"&
        .or.fileRow(1:4).eq."DEMB".or.fileRow(1:4).eq."SEMB")Then
! Store last function start/length
        If(functionCounter.gt.0)Then
          eamKey(eamFunctionKey,4) = functionStart
          eamKey(eamFunctionKey,5) = functionLength
          eamKey(eamFunctionKey,6) = (functionStart + functionLength - 1)
          FunctionStart = functionStart + functionLength
          FunctionLength = 0
        End If
! Store this function elements and type
        Read(fileRow,*) functionType, elementA            ! Read in line
        FunctionCounter = functionCounter + 1             ! Increment function counter
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
        FunctionLength = functionLength + 1
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
  End Subroutine readFileE

  Subroutine eamDerivatives()
! Fills in first and second order derivatives
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j, k, xPoint, pointOffset
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:4,1:2) :: pointsInterp
! Calculate y'(x) and y''(x)
    FunctionCounter = 0
    Do i=1,size(eamKey,1)
      If(eamKey(i,1).gt.0)Then
        FunctionCounter = functionCounter + 1
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
! Set node positions using the currently loaded eam functions in eamKey/eamData
! Pair functions: ZBL hard core, exp(ax
! 
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
    Integer(kind=StandardInteger) :: eamStart, eamLength, eamEnd
    Integer(kind=StandardInteger) :: nodes, nodeKey, functionType
    Real(kind=DoubleReal) :: x, xStart, xEnd
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Integer(kind=StandardInteger) :: zA, zB
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
        FunctionType = eamKey(i,3)
        eamStart = eamKey(i,4)
        eamLength = eamKey(i,5)
        eamEnd = eamKey(i,6)
! xStart xEnd
        xStart = eamData(eamStart,1)       
        xEnd = eamData(eamEnd,1) 
        If(eamKey(i,3).eq.1)Then ! Function type = pair
          xStart = zblHardCore(1)
          zA = elementsCharge(eamKey(i,1))
          zB = elementsCharge(eamKey(i,2))          
        End If
! Loop through nodes for each function
        Do j=1,nodes
          nodeKey = nodeKey + 1
!If first node, save node key data
          If(j.eq.1)Then  
            splineNodesKey(i,1) = eamKey(i,1)           ! Type A
            splineNodesKey(i,2) = eamKey(i,2)           ! Type B
            splineNodesKey(i,3) = eamKey(i,3)           ! Function type
            splineNodesKey(i,4) = nodeKey               ! Start
            splineNodesKey(i,5) = nodes                 ! Length
            splineNodesKey(i,6) = nodeKey + nodes - 1   ! End
          End If
! x point
          x = xStart+1.0D0*(j-1)*((xEnd-xStart)/(nodes-1))
! last node and Pair/Dens type 0,0,0
          If(j.eq.nodes.and.&
          (eamKey(i,3).eq.1.or.eamKey(i,3).eq.2.or.eamKey(i,3).eq.4.or.eamKey(i,3).eq.5)&
          )Then
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = 0.0D0
            splineNodesData(nodeKey,3) = 0.0D0
            splineNodesData(nodeKey,4) = 0.0D0
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0
          ElseIf(j.eq.1.and.eamKey(i,3).eq.1)Then  ! Pair - get ZBL y(x), y'(x) and y''(x)
            yArray = ZblFull (x, zA, zB)            
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = yArray(1)
            splineNodesData(nodeKey,3) = yArray(2)
            splineNodesData(nodeKey,4) = yArray(3)
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0        
          Else
            yArray = PointInterp(eamData,x,eamInterpPoints,2,eamStart,eamLength)
            splineNodesData(nodeKey,1) = x
            splineNodesData(nodeKey,2) = yArray(1)
            splineNodesData(nodeKey,3) = yArray(2)
            splineNodesData(nodeKey,4) = yArray(3)
            splineNodesData(nodeKey,5) = 1.0D0*j
            splineNodesData(nodeKey,6) = 0.0D0
          End If
        End Do
      End If
      If(functionCounter.eq.eamFunctionCount )Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do
! Store total number of nodes
    splineTotalNodes = nodeKey
  End Subroutine setEamNodes
  
  
  Subroutine rescaleEmbedNodes(newDensityLimit)
! Set nodes
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: rescFactor, newDensityLimit
    Do i=1,size(splineNodesKey,1)
      print *,i,splineNodesKey(i,1),splineNodesKey(i,3)
      If(splineNodesKey(i,1).gt.0)Then    
        If(splineNodesKey(i,3).eq.3)Then
          rescFactor = newDensityLimit/splineNodesKey(i,6)
          print *,i,newDensityLimit,splineNodesKey(i,6),rescFactor
          Do j=splineNodesKey(i,4), splineNodesKey(i,6)
            splineNodesData(j,1) = rescFactor * splineNodesData(j,1) 
          End Do
        End If
      Else
        Exit
      End If
    End Do 
    
    
  End Subroutine rescaleEmbedNodes

  

  Subroutine countEamNodes(totalNodes)
! Set nodes
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: totalNodes, i, j
! Count nodes    
    totalNodes = 0
    Do i=1,size(splineNodesKey,1)
      If(splineNodesKey(i,1).gt.0)Then    
        Do j=splineNodesKey(i,4), splineNodesKey(i,6)
          totalNodes = totalNodes + 1
        End Do
      Else
        Exit
      End If
    End Do 
  End Subroutine countEamNodes

  Subroutine setEamSpline(calcDerivIn)
! Spline between nodes and store functions in eamKey/eamData
! Force ZBL Core for pair potentials
! exp(a+bx) style spline from  zbl to first node
    Implicit None   ! Force declaration of all variables
! In    
    Logical, Optional :: calcDerivIn
    Logical :: calcDeriv
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
    Integer(kind=StandardInteger) :: nodes, nodeStart, nodeLength, nodeEnd, eamPoint
    Integer(kind=StandardInteger) :: pointsPerFunction
    Real(kind=DoubleReal), Dimension(1:1001,1:4) :: splineDataPoints
    Real(kind=DoubleReal) :: x, changeX
    Integer(kind=StandardInteger) :: zA, zB
    Integer(kind=StandardInteger) :: pointsZBL, pointsSpline
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Integer(kind=StandardInteger), Dimension(1:1000) :: splineType
    Logical, Dimension(1:1000) :: forceCalcDerv, interpNode
! Optional argument    
    calcDeriv = .false.
    If(Present(calcDerivIn))Then
      calcDeriv = calcDerivIn
    End If
! Init variables
    eamKey = 0                 ! clear eam key
    eamData = 0.0D0            ! clear eam data
    splineDataPoints = 0.0D0   ! init temp spline points array
    pointsPerFunction = 1001   ! number of data points per function
    forceCalcDerv = .false.
    interpNode = .true.
! Loop through EAM functions
    functionCounter = 0
    Do i=1,size(splineNodesKey,1)    
      If(splineNodesKey(i,1).gt.0)Then ! Check that eam function is stored    
! Count function      
        functionCounter = functionCounter + 1
! Update Key Data        
        eamKey(functionCounter,1) = splineNodesKey(i,1)    ! Species 1
        eamKey(functionCounter,2) = splineNodesKey(i,2)    ! Species 2
        eamKey(functionCounter,3) = splineNodesKey(i,3)    ! Function type
        eamKey(functionCounter,4) = 1+(pointsPerFunction*(functionCounter-1))              ! function start key
        eamKey(functionCounter,5) = pointsPerFunction                                      ! function length
        eamKey(functionCounter,6) = eamKey(functionCounter,4)+eamKey(functionCounter,5)-1  ! function end key
! Nodes used in spline
        nodes = splineNodeCount(splineNodesKey(i,3))
! ----------------------------------------------------------------------
! Pair Potentials
! ----------------------------------------------------------------------
        If(splineNodesKey(i,3).eq.1)Then       
! ZBL
          zA = elementsCharge(splineNodesKey(i,1))
          zB = elementsCharge(splineNodesKey(i,2))   
! Set node start/end points
          nodeStart = splineNodesKey(i,4)
          nodeLength = splineNodesKey(i,5)
          nodeEnd = splineNodesKey(i,6)
! what part of function will be spline - function runs from 0.0 to splineNodesData(nodeEnd,1)
          changeX = splineNodesData(nodeEnd,1)-splineNodesData(nodeStart,1)
! determine 1 to nodes ZBL, nodes ZBL onwards
          pointsSpline = ceiling((changeX/splineNodesData(nodeEnd,1))*pointsPerFunction) ! data points in spline section
          pointsZBL = pointsPerFunction - pointsSpline                                   ! data points in zbl section          
! ZBL section
          Do j=1,pointsZBL
            x = (j-1)*(splineNodesData(nodeEnd,1)/(1.0D0*pointsPerFunction))
            yArray = ZblFull (x, zA, zB)
            eamPoint = eamKey(functionCounter,4)+j-1
            eamData(eamPoint,1) = x
            eamData(eamPoint,2) = yArray(1)
            eamData(eamPoint,3) = yArray(2)
            eamData(eamPoint,4) = yArray(3)
          End Do      
! set spline types          
          splineType = 1    ! default to poly
          splineType(1) = 2 ! first segment exp(poly) [3rd order]
! leave node 1, 2 and end untouched (exp poly spline between 1-2 and end needs to be 0,0,0)        
          forceCalcDerv = .true.
          forceCalcDerv(1) = .false.
          forceCalcDerv(nodeEnd-nodeStart+1) = .false.
! choose which nodes to
          interpNode = .true.
          interpNode(1) = .false.
! Spline section
          splineDataPoints = SplineNodes(splineNodesData,pointsSpline,nodeStart,nodeEnd,&
          pointsPerFunction,splineType,forceCalcDerv,interpNode)
          Do j=1,pointsSpline          
            eamPoint = eamKey(functionCounter,4)+j+pointsZBL-1
            eamData(eamPoint,1) = splineDataPoints(j,1)
            eamData(eamPoint,2) = splineDataPoints(j,2)
            eamData(eamPoint,3) = splineDataPoints(j,3)
            eamData(eamPoint,4) = splineDataPoints(j,4)
          End Do
! ----------------------------------------------------------------------
! Dens + Embe       
! ----------------------------------------------------------------------   
        Else        
! Set node start/end points
          nodeStart = splineNodesKey(i,4)
          nodeLength = splineNodesKey(i,5)
          nodeEnd = splineNodesKey(i,6)
! spline between nodes       
          splineType = 1    ! default to poly
! Force density fit to analytic type          
          If(eamKey(functionCounter,3).eq.2.or.&
          eamKey(functionCounter,3).eq.4.or.eamKey(functionCounter,3).eq.5)Then
            If(optDensityFit.eq.1)Then
              splineType = 6
            End If
          End If 
! Force embed fit to analytic type          
          If(eamKey(functionCounter,3).eq.3.or.&
          eamKey(functionCounter,3).eq.6.or.eamKey(functionCounter,3).eq.7)Then
            If(optEmbeddingFit.eq.1)Then
              splineType = 7
            End If
          End If  
! force calculation of derivative/second derivative by point interpolation (except start and end nodes)          
          forceCalcDerv = .true.
          forceCalcDerv(1) = .false.
          forceCalcDerv(nodeEnd-nodeStart+1) = .false.
! spline
          splineDataPoints = SplineNodes(splineNodesData,pointsPerFunction,nodeStart,nodeEnd,&
          pointsPerFunction,splineType,forceCalcDerv)
! copy data into eamData array
          Do j=1,pointsPerFunction
            eamPoint = eamKey(functionCounter,4)+j-1
            eamData(eamPoint,1) = splineDataPoints(j,1)
            eamData(eamPoint,2) = splineDataPoints(j,2)
            eamData(eamPoint,3) = splineDataPoints(j,3)
            eamData(eamPoint,4) = splineDataPoints(j,4)
          End Do
        End If  
      End If
      If(functionCounter.eq.eamFunctionCount)Then
        Exit  ! Exit, all functions cycled through
      End If
    End Do    
  End Subroutine setEamSpline



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
    If(zblHardCore(1).gt.0.0D0)Then
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
            !Call outputZBL(zA, zB, pointA, pointB, splineCoeffs)
          End If
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
    End If
  End Subroutine eamZblHardCore
  
  Subroutine eamPairZbl()
! Force ZBL Core
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j, zA, zB
    Real(kind=DoubleReal) :: xA, xB
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:6) :: splineCoeffs
! Loop through EAM functions
    If(zblHardCore(1).gt.0.0D0)Then
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
            !Call outputZBL(zA, zB, pointA, pointB, splineCoeffs)
          End If
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
    End If
  End Subroutine eamPairZbl

  Subroutine outputSummary()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, functionCounter
! Print out
    If(TerminalPrint())Then
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
            elements(eamKey(i,2)),"(",eamKey(i,2),") ",&
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

  Subroutine eamDataDump()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i, j
! Create output file
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"eam.dat"))
! Calculate y'(x) and y''(x)
      functionCounter = 0
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
          write(999,*) "EAM Function ",functionCounter,eamKey(i,1),eamKey(i,2),&
          eamKey(i,3),eamKey(i,4),eamKey(i,5),eamKey(i,6)
          Do j=eamKey(i,4),eamKey(i,6)
            write(999,*) j,eamData(j,1),eamData(j,2),eamData(j,3),eamData(j,4)
          End Do
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
    End If
  End Subroutine eamDataDump
  
  Subroutine eamCharts(prefixIn)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! In    
    Character(*), Optional :: prefixIn
! Private variables
    Integer(kind=StandardInteger) :: functionCounter, i
    !Integer(kind=StandardInteger), Dimension(1:1) :: fittingArr
    Character(len=30) :: fileName
    Type(plotData) :: eamPlotData
    Character(Len=16) :: label
    Character(len=16) :: prefix
! Optional arguments
    prefix = "                "
    If(Present(prefixIn))Then
      prefix = prefixIn
    End If    
! Create output file
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"eam.dat"))
! Calculate y'(x) and y''(x)
      functionCounter = 0
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
          write(fileName,"(I4)") functionCounter
          fileName = adjustl(trim(prefix))//"EAM_"//adjustl(trim(fileName))
          fileName = trim(fileName)//"_"//eamFunctionTypes(eamKey(i,3))
! Chart settings    
          Call plotInit(eamPlotData)
          eamPlotData%tempDirectory = trim(tempDirectory)
          eamPlotData%outputDirectory = trim(outputDirectory)
          eamPlotData%outputName = trim(fileName)
          eamPlotData%title = "Plot for "//eamFunctionTypes(eamKey(i,3))
          eamPlotData%xAxis = "Radius"
          eamPlotData%yAxis = "Energy"  
          eamPlotData%cleanPyFile = .true.
          If(eamKey(i,3).eq.1)Then          
            eamPlotData%yMin = -10.0D0 
            eamPlotData%yMax = 100.0D0 
          Else
            eamPlotData%yMin = 1.1D99 
            eamPlotData%yMax = -1.1D99
          End If                    
          Call plotAdd(eamPlotData, eamData, label, "", eamKey(i,4), eamKey(i,6), 1, 2)
          Call plotMake(eamPlotData)  
          If(eamKey(i,3).eq.1)Then  ! Do a "close up" of the pair potential            
            write(fileName,"(I4)") functionCounter
            fileName = "Input_EAM_"//adjustl(trim(fileName))
            fileName = trim(fileName)//"_"//eamFunctionTypes(eamKey(i,3))//"_C"
! Set chart values
            eamPlotData%outputName = trim(fileName)
            eamPlotData%yMin = -1.0D0 
            eamPlotData%yMax = 10.0D0 
            Call plotMake(eamPlotData)   
          End If
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
    End If
  End Subroutine eamCharts

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
  
  
  

  Subroutine loadInputEAM()
    Implicit None   ! Force declaration of all variables
! Load from input arrays
    eamKey = eamKeyInput
    eamData = eamDataInput
  End Subroutine loadInputEAM
  

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

End Module readEAM
