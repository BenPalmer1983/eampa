Module clean

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
  
!declare global variables  

!Privacy of functions/subroutines/variables
  Private
!Variables  
!Subroutines  
  Public :: runClean   
!Functions
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

 
!------------------------------------------------------------------------!
! Run Clean Subroutines
!------------------------------------------------------------------------!  
  
  Subroutine runClean()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i
	Character(len=512) :: testLine
!Add to list
	Do i=1,100
	  testLine = fileCleanupList(i)
	  If(testLine(1:2).eq."  ")Then
		Exit
	  Else
!clear out files
        Call system("rm -f "//trim(adjustl(testLine)))	
	  End If
	End Do 
	
  
  End Subroutine runClean
  
  
  Subroutine cleanFiles()
!force declaration of all variables
	Implicit None	
!declare private variables
    
  
  
  End Subroutine cleanFiles
  
  
  
  
  

End Module clean