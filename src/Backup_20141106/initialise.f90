Module initialise

! Setup Modules
  Use kinds
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use globals        ! declare all globals

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h' 


!Privacy of functions/subroutines/variables
  Private
  
!Variables


!Subroutines
  Public :: runInitialise        !Subroutine
  Public :: fileToClean
  
!Functions
  Public :: ProgramTime             !Function
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines
  Subroutine runInitialise()
! Internal subroutine variables
    Integer(kind=StandardInteger) :: error
    Character(len=255) :: inputFileName
    Integer(kind=StandardInteger), Dimension(1:3) :: theTime, theDate
! Initialise variables  
    inputFileName = BlankString(inputFileName)
! Check for input file, exit if no file      
    Call get_command_argument(1,inputFileName)
    If(inputFileName(1:4).eq."    ")Then
      If(mpiProcessID.eq.0)Then
        Print *,"No input file, exiting."
      End If
      Call Exit(0)
    End If
!Set random seed    
    Call SetRandomSeedArray()        
!get the working directory
    Call getcwd(currentWorkingDirectory)  
!Set output and temp/scratch directories
    outputDirectory = Trim(currentWorkingDirectory)//"/output"
    tempDirectory = Trim(currentWorkingDirectory)//"/temp"
! MPI variables (public)
    Call MPI_Comm_size( MPI_COMM_WORLD ,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)  
!Call date subroutines
    Call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    Call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
! Print run message    
    If(mpiProcessID.eq.0)Then
      print *,"Initialising program."
      print *,"MPI Processes: ",mpiProcessCount
      print "(A7,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)", &
      " Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)  
      print "(A30,F8.4,A5,F8.4,A3)","Large array memory allocated: ",&
        largeArraySize,"MB  (",(mpiProcessCount*largeArraySize),"MB)"
    End If
!Call init subroutines
    Call initVars()
    Call makeDirectories()
    Call initDataFiles()
  End Subroutine runInitialise  
  
!Read initVars
  Subroutine initVars()  
!declare private variables  
    Integer(kind=StandardInteger) :: i,j
    Character(len=512) :: blankLine  
    Do j=1,512
      blankLine(j:j) = " "
    End Do  
    Do i=1,100
      fileCleanupList(i) = blankLine
    End Do    
  End Subroutine initVars  
  
    
!Make required directories
  Subroutine makeDirectories()  
!Internal subroutine variables
!output directory
    Call makeDir(outputDirectory)    
    Call makeDir(tempDirectory)  
  End Subroutine makeDirectories  
  

!Run all the input subroutines
  Subroutine initDataFiles()  
!Internal subroutine variables
  Integer(kind=StandardInteger), Dimension(1:3) :: theTime, theDate
!Call date subroutines
    Call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    Call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
!save output file name
    outputFile = trim(outputDirectory)//"/"//"output.dat"
!Create output file
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(outputFile))
      write(999,"(A38)") "======================================"
      write(999,"(A38)") "            Output File               "
      write(999,"(A38)") "      University of Birmingham        "
      write(999,"(A38)") "             Ben Palmer               "
      write(999,"(A38)") "======================================"
      write(999,"(A1)") " "
      write(999,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)  
      write(999,"(A1)") " "
      write(999,"(A15,I4)") "MPI Processes: ",mpiProcessCount
      write(999,"(A1)") " "
      write(999,"(A1)") " "
      write(999,"(A19,F8.4)") "Load globals time: ",(globalsTimeEnd-globalsTimeStart)
      write(999,"(A31,F8.4,A5,F8.4,A3)") "Large array memory allocated: ",&
        largeArraySize,"MB  (",(mpiProcessCount*largeArraySize),"MB)"
      write(999,"(A1)") " "
!close output file
      close(999)
    End If
!save output file name
    outputFileForces = trim(outputDirectory)//"/"//"forces.dat"
!Create output file
    If(mpiProcessID.eq.0)Then
    open(unit=989,file=trim(outputFileForces))
    write(989,"(A38)") "======================================"
    write(989,"(A38)") "            Forces File               "
    write(989,"(A38)") "      University of Birmingham        "
    write(989,"(A38)") "             Ben Palmer               "
    write(989,"(A38)") "======================================"
    write(989,"(A1)") " "
    write(989,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
    "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)  
    write(989,"(A1)") " "
    write(989,"(A1)") " "
!close output file
    close(989)
  End If   
!Create output file
    If(mpiProcessID.eq.0)Then
    open(unit=979,file=trim(trim(outputDirectory)//"/"//"rssLog.dat"))
    write(979,"(A38)") "======================================"
    write(979,"(A38)") "            RSS Log File              "
    write(979,"(A38)") "      University of Birmingham        "
    write(979,"(A38)") "             Ben Palmer               "
    write(979,"(A38)") "======================================"
    write(979,"(A1)") " "
    write(979,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
    "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)  
    write(979,"(A1)") " "
    write(979,"(A1)") " "
!close output file
    close(979)
  End If   
!save output file name
  outputFileForces = trim(outputDirectory)//"/"//"nlSeparation.dat"
!Create output file
    If(mpiProcessID.eq.0)Then
      open(unit=969,file=trim(outputFileForces))
      write(969,"(A38)") "======================================"
      write(969,"(A38)") "     Neighbour List Separations       "
      write(969,"(A38)") "      University of Birmingham        "
      write(969,"(A38)") "             Ben Palmer               "
      write(969,"(A38)") "======================================"
      write(969,"(A1)") " "
      write(969,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)  
      write(969,"(A1)") " "
      write(969,"(A1)") " "
!close output file
      close(969)
    End If 
!Create output file
    If(mpiProcessID.eq.0)Then
      open(unit=969,file=trim(outputDirectory)//"/"//"EoSFitting.dat")
      write(969,"(A38)") "======================================"
      write(969,"(A38)") "         EoS Fitting Results          "
      write(969,"(A38)") "      University of Birmingham        "
      write(969,"(A38)") "             Ben Palmer               "
      write(969,"(A38)") "======================================"
      write(969,"(A1)") " "
      write(969,"(A6,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I4.4)") &
      "Date: ",theTime(1),":",theTime(2)," ",theDate(1),"/",theDate(2),"/",theDate(3)  
      write(969,"(A1)") " "
      write(969,"(A1)") " "
!close output file
      close(969)
    End If 
  End Subroutine initDataFiles
  
  
    
  
!Other subroutines  
!Run all the input subroutines
  Subroutine fileToClean(fileName)  
!Internal subroutine variables
  Integer(kind=StandardInteger) :: i
  Character(*) :: fileName
  Character(len=512) :: testLine
!Add to list
  Do i=1,100
    testLine = fileCleanupList(i)
    If(testLine(1:2).eq."  ")Then
      fileCleanupList(i) = trim(adjustl(fileName))
    Exit
    End If
  End Do  
  
  
  End Subroutine fileToClean
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!

  Function ProgramTime () RESULT (outputTime)
    ! -- Argument and result
    Real(kind=DoubleReal) :: inputTime, outputTime 
    Call cpu_time(inputTime)
    outputTime = inputTime - programStartTime
  End Function ProgramTime    

  Function OutputProgramTime () RESULT (outputTime)
    ! -- Argument and result
    Real(kind=DoubleReal) :: outputTime 
    open(unit=9991,file=trim(outputFile),status="old",position="append",action="write")
    outputTime = ProgramTime()
    write(9991,"(A16,F8.4)") "Program time:   ",outputTime
    close(9991)
  End Function OutputProgramTime    

  Function ClockTime () RESULT (outputTime)
    ! -- Argument and result
    Real(kind=DoubleReal) :: outputTime 
    Call cpu_time(outputTime)
  End Function ClockTime    

End Module initialise