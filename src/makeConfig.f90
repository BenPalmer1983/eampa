Module makeConfig

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
  Public :: makeConfigFile

  Contains
  Subroutine makeConfigFile(makeOption, aLatIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: makeOption
    Real(kind=DoubleReal) :: aLatIn
    Real(kind=DoubleReal) :: timeStartMC, timeEndMC
! Start Time
    Call cpu_time(timeStartMC)
! Set configuration file
    configFilePath = Trim(tempDirectory)//"/automatedConfig.conf"
! Set global unit vector
    globalConfigUnitVector(1,1) = 1.0D0
    globalConfigUnitVector(1,2) = 0.0D0
    globalConfigUnitVector(1,3) = 0.0D0
    globalConfigUnitVector(2,1) = 0.0D0
    globalConfigUnitVector(2,2) = 1.0D0
    globalConfigUnitVector(2,3) = 0.0D0
    globalConfigUnitVector(3,1) = 0.0D0
    globalConfigUnitVector(3,2) = 0.0D0
    globalConfigUnitVector(3,3) = 1.0D0
! Only make config file on root process
    If(mpiProcessID.eq.0)Then
! Prepare config file
      Open(UNIT=140,FILE=Trim(tempDirectory)//"/automatedConfig.conf")
      write(140,"(A23)") "! Automatic config file"
      Close(140)
! Make FCC for lattice parameter testing
      If(makeOption.eq.1)Then
        Call makeFCC(aLatIn,6.5D0,3,3,3)
      End If
! Make BCC for lattice parameter testing
      If(makeOption.eq.2)Then
        Call makeBCC(aLatIn,6.5D0,4,4,4)
      End If
    End If
! Synch MPI processes
    Call M_synchProcesses()
! End Time
    Call cpu_time(timeEndMC)
! Store Time
    Call storeTime(9,timeEndMC-timeStartMC)
  End Subroutine makeConfigFile
! ---------------------------------------------------------------------------------------------------
  Subroutine makeFCC(aLat,rCut,xC,yC,zC)
    Real(kind=DoubleReal) :: aLat,rCut,xCoord,yCoord,zCoord
    Integer(kind=StandardInteger) :: xC,yC,zC
    Integer(kind=StandardInteger) :: i, n
    Real(kind=DoubleReal), Dimension(1:4,1:3) :: unitCell
! Init variables
    unitCell(1,1) = 0.0D0
    unitCell(1,2) = 0.0D0
    unitCell(1,3) = 0.0D0
    unitCell(2,1) = 0.5D0
    unitCell(2,2) = 0.5D0
    unitCell(2,3) = 0.0D0
    unitCell(3,1) = 0.5D0
    unitCell(3,2) = 0.0D0
    unitCell(3,3) = 0.5D0
    unitCell(4,1) = 0.0D0
    unitCell(4,2) = 0.5D0
    unitCell(4,3) = 0.5D0
! Open config file
    Open(UNIT=140,FILE=Trim(tempDirectory)//"/automatedConfig.conf",&
    status="old",position="append",action="write")
    Call fileToClean(Trim(tempDirectory)//"/automatedConfig.conf")
    Do i=1,elementsCount
      If(elements(i).ne."ZZ")Then
! Write header
        write(140,"(A4)") "#NEW"
        write(140,"(A6,F12.6)")  "#LP   ",aLat
        write(140,"(A39)") "#X    1.0000000   0.0000000   0.0000000"
        write(140,"(A39)") "#Y    0.0000000   1.0000000   0.0000000"
        write(140,"(A39)") "#Z    0.0000000   0.0000000   1.0000000"
        write(140,"(A6,I4,I4,I4)")  "#CC   ",xC,yC,zC
        write(140,"(A6,F12.6)")  "#RC   ",rCut
! Write co-ords
        Do n=1,size(unitCell,1)
          xCoord = 1.0D0*unitCell(n,1)
          yCoord = 1.0D0*unitCell(n,2)
          zCoord = 1.0D0*unitCell(n,3)
          write(140,"(A4,F12.6,F12.6,F12.6)") elements(i),xCoord,yCoord,zCoord
        End Do
        write(140,"(A4)") "#END"
      Else
        Exit
      End If
    End Do
    Close(140)
  End Subroutine makeFCC
! ---------------------------------------------------------------------------------------------------
  Subroutine makeBCC(aLat,rCut,xC,yC,zC)
    Real(kind=DoubleReal) :: aLat,rCut,xCoord,yCoord,zCoord
    Integer(kind=StandardInteger) :: xC,yC,zC
    Integer(kind=StandardInteger) :: i, n
    Real(kind=DoubleReal), Dimension(1:2,1:3) :: unitCell
! Init variables
    unitCell(1,1) = 0.0D0
    unitCell(1,2) = 0.0D0
    unitCell(1,3) = 0.0D0
    unitCell(2,1) = 0.5D0
    unitCell(2,2) = 0.5D0
    unitCell(2,3) = 0.5D0
! Open config file
    Open(UNIT=140,FILE=Trim(tempDirectory)//"/automatedConfig.conf",&
    status="old",position="append",action="write")
    Do i=1,elementsCount
      If(elements(i).ne."ZZ")Then
! Write header
        write(140,"(A4)") "#NEW"
        write(140,"(A6,F12.6)")  "#LP   ",aLat
        write(140,"(A39)") "#X    1.0000000   0.0000000   0.0000000"
        write(140,"(A39)") "#Y    0.0000000   1.0000000   0.0000000"
        write(140,"(A39)") "#Z    0.0000000   0.0000000   1.0000000"
        write(140,"(A6,I4,I4,I4)")  "#CC   ",xC,yC,zC
        write(140,"(A6,F12.6)")  "#RC   ",rCut
! Write co-ords
        Do n=1,size(unitCell,1)
          xCoord = 1.0D0*unitCell(n,1)
          yCoord = 1.0D0*unitCell(n,2)
          zCoord = 1.0D0*unitCell(n,3)
          write(140,"(A4,F12.6,F12.6,F12.6)") elements(i),xCoord,yCoord,zCoord
        End Do
        write(140,"(A4)") "#END"
      Else
        Exit
      End If
    End Do
    Close(140)
  End Subroutine makeBCC

End Module makeConfig
