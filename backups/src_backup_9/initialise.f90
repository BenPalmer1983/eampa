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
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger), Dimension(1:3) :: theTime, theDate
	Character(len=255) :: outputFilePath
	Integer(kind=StandardInteger), Dimension(1:10) :: rnSeed
    Integer(kind=StandardInteger) :: error  
	
!MPI
	Call MPI_Init(error)
    Call MPI_Comm_size( MPI_COMM_WORLD ,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)
	
!store start time
	Call cpu_time(programStartTime)
	Call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    Call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
	
!get the working directory
	Call getcwd(currentWorkingDirectory)
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
	
!Set random seed
    Call SetRandomSeedArray()
	
  End Subroutine runInitialise
  
  
    
  
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