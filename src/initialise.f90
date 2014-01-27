Module initialise

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Real(kind=SingleReal) :: programStartTime
  Character(len=255) :: currentWorkingDirectory
  Character(len=255) :: outputFile

!Privacy of functions/subroutines/variables
  Private
  Public :: programStartTime		!Variable
  Public :: outputFile      		!Variable
  Public :: currentWorkingDirectory	!Variable
  Public :: runInitialise		    !Subroutine
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
	
!store start time
	Call cpu_time(programStartTime)
	Call idate(theDate)   ! theDate(1)=day, (2)=month, (3)=year
    Call itime(theTime)   ! theDate(1)=hour, (2)=minute, (3)=second
	
!get the working directory
	Call getcwd(currentWorkingDirectory)
	
	!Create output file
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
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
	
!random number seed from cpu
    do i=1,10
	  Call SYSTEM_CLOCK(rnSeed(i))
	enddo
    Call RANDOM_SEED(put=rnSeed)
	
!close output file
	close(999)

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
  

End Module initialise