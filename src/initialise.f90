Module initialise

! Setup Modules
  Use mpi
  Use kinds
  Use constants
  Use general
  Use strings		!string functions
  Use maths

!Force declaration of all variables
  Implicit None
  
!Declare Global Variables  
  Real(kind=SingleReal) :: programStartTime
  Character(len=255) :: currentWorkingDirectory
  Character(len=255) :: outputFile
  Character(len=255) :: outputFileEnergies
  Character(len=255) :: outputFileForces
  Character(len=512), Dimension(1:100) :: fileCleanupList
  Character(len=255) :: outputDirectory
  Character(len=255) :: tempDirectory
  !Character(len=255) :: dataDirectory
!MPI Global Variables
  Integer(kind=StandardInteger) :: mpiProcessCount, mpiProcessID  

!Privacy of functions/subroutines/variables
  Private
  
!Variables
  Public :: programStartTime		!Variable
  Public :: outputFile      		!Variable
  Public :: outputFileEnergies    	!Variable
  Public :: outputFileForces    	!Variable
  Public :: currentWorkingDirectory	!Variable
  Public :: fileCleanupList
  Public :: outputDirectory
  Public :: tempDirectory
  !Public :: dataDirectory
!MPI Variables  
  Public :: mpiProcessCount
  Public :: mpiProcessID

!Subroutines
  Public :: runInitialise		    !Subroutine
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
!Internal subroutine variables
    Integer(kind=StandardInteger) :: error
!Set random seed		
	Call SetRandomSeedArray()      	
!get the working directory
	Call getcwd(currentWorkingDirectory)	
!Set output and temp/scratch directories
    outputDirectory = Trim(currentWorkingDirectory)//"/output"
    tempDirectory = Trim(currentWorkingDirectory)//"/temp"
!MPI variables (public)
	Call MPI_Init(error)
    Call MPI_Comm_size( MPI_COMM_WORLD ,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)	
!store start time
	Call cpu_time(programStartTime)
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
	Integer(kind=StandardInteger) :: i
!output directory
    Call makeDir(outputDirectory)    
	Call makeDir(tempDirectory)
	
  End Subroutine makeDirectories  
  

!Run all the input subroutines
  Subroutine initDataFiles()	
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
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
!close output file
	  close(999)
	End If
!save output file name
	outputFileForces = trim(outputDirectory)//"/"//"outputForces.dat"
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
!save output file name
	outputFileForces = trim(outputDirectory)//"/"//"rssLog.dat"
!Create output file
    If(mpiProcessID.eq.0)Then
	  open(unit=979,file=trim(outputFileForces))
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

function ProgramTime () RESULT (outputTime)
    ! -- Argument and result
	Real(kind=SingleReal) :: inputTime, outputTime 
	Call cpu_time(inputTime)
	outputTime = inputTime - programStartTime
End Function ProgramTime    

function OutputProgramTime () RESULT (outputTime)
    ! -- Argument and result
	Real(kind=SingleReal) :: outputTime 
	open(unit=9991,file=trim(outputFile),status="old",position="append",action="write")
	outputTime = ProgramTime()
	write(9991,"(A16,F8.4)") "Program time:   ",outputTime
	close(9991)
End Function OutputProgramTime    

function ClockTime () RESULT (outputTime)
    ! -- Argument and result
	Real(kind=DoubleReal) :: outputTime 
	Call cpu_time(outputTime)
End Function ClockTime    

End Module initialise