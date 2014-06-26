Module initialise

! Setup Modules
  Use mpi
  Use kinds
  Use constants
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
!Ini Variables
  Character(len=255) :: outputDirectory
  Character(len=255) :: tempDirectory
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
!Ini Variables
  Public :: outputDirectory
  Public :: tempDirectory
!MPI Variables  
  Public :: mpiProcessCount
  Public :: mpiProcessID

!Subroutines
  Public :: runInitialise		    !Subroutine
  
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
!get the working directory
	Call getcwd(currentWorkingDirectory)	
!MPI variables (public)
	Call MPI_Init(error)
    Call MPI_Comm_size( MPI_COMM_WORLD ,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)	
!store start time
	Call cpu_time(programStartTime)
!Call init subroutines
    Call readInitFile()
    Call initDataFiles()
	Call SetRandomSeedArray()      !Set random seed		
  End Subroutine runInitialise	
	
!Read INIT file
  Subroutine readInitFile()	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
    Integer(kind=StandardInteger) :: i, ios, rowCount  
	Character(len=128) :: fileRow
	Character(len=128) :: bufferA, bufferB, bufferC
!Open ini file and read	
	Open(unit=10,file="eampa.ini")
	rowCount = 0
	Do i=1,maxFileRows 
      rowCount = rowCount + 1
	  Read(1,"(A128)",IOSTAT=ios) fileRow
	  If (ios /= 0) Then
	    EXIT 
	  End If
	End Do  
	Close(10)
!Set defaults
    outputDirectory = Trim(currentWorkingDirectory)//"/output"
    tempDirectory = Trim(currentWorkingDirectory)//"/temp"
!Open ini file and reread	
	Open(unit=10,file="eampa.ini")
	Do i=1,rowCount 	
	  Read(1,"(A128)",IOSTAT=ios) fileRow
	  If(StrToUpper(fileRow(1:10)).eq."#OUTPUTDIR")Then
	    Read(fileRow,*) bufferA, bufferB
		If(StrToUpper(bufferB(1:7)).ne."DEFAULT")Then
		  outputDirectory = bufferB
		End If
	  End If
	End Do
	Close(10)
	

	
  End Subroutine readInitFile	


!Run all the input subroutines
  Subroutine initDataFiles()	
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger), Dimension(1:3) :: theTime, theDate
!Call date subroutines
	Call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    Call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
!save output file name
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
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
	  write(999,"(A1)") " "
!close output file
	  close(999)
	End If
!save output file name
	outputFileForces = trim(currentWorkingDirectory)//"/"//"outputForces.dat"
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