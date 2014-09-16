Module strings

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
  Private
  Public :: StrToUpper
  Public :: NumericOnly
  Public :: RemoveSpaces
  Public :: CorrectFilePath
  Public :: TrimSpaces
  
  character( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  

  function StrToUpper (input) RESULT (output)
    ! -- Argument and result
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: output
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, n

    ! -- Copy input string
    output = input
    ! -- Loop over string elements
    DO i = 1, LEN( output )
      ! -- Find location of letter in lower case constant string
      n = INDEX( LOWER_CASE, output( i:i ) )
      ! -- If current substring is a lower case letter, make it upper case
      IF ( n /= 0 ) output( i:i ) = UPPER_CASE( n:n )
    END DO
  END FUNCTION StrToUpper   
  
  function NumericOnly (input) RESULT (output)
    ! -- Argument and result
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: outputTemp
    CHARACTER(LEN(input)) :: output
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, j, n

    ! -- Copy input string
    outputTemp = input
	
	DO i = 1, LEN( outputTemp )
	  output( i:i ) = " "
	End Do
	n = 0	
	DO i = 1, LEN( outputTemp )
	  if(outputTemp( i:i ).eq.".".or.(iachar(outputTemp( i:i )).ge.48.and.iachar(outputTemp( i:i )).le.57))then
	    n = n + 1
	    output( n:n ) = outputTemp( i:i )
	  else
	    output( i:i ) = " "
	  endif
	End Do
	  
  END FUNCTION NumericOnly   
  
  
  function RemoveSpaces (input) RESULT (output)
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: outputTemp
    CHARACTER(LEN(input)) :: output
	
	! -- Local variables
    Integer(kind=StandardInteger) :: i, j, n

    ! -- Copy input string
    outputTemp = input
	
	!Blank output
	DO i = 1, LEN( outputTemp )
	  output( i:i ) = " "
	End Do
	
	!transfer outputtemp to output without spaces
	j = 0
	DO i = 1, LEN( outputTemp )
	  if(outputTemp( i:i ).ne." ")then
	    j = j + 1
		output( j:j ) = outputTemp( i:i )
	  endif
	End Do
	
  
  END FUNCTION RemoveSpaces   
  
  
  Function TrimSpaces (input) RESULT (output)
    CHARACTER(*), INTENT(IN) :: input
	CHARACTER(LEN(trim(adjustl(input)))) :: output
    output = trim(adjustl(input))
  End Function TrimSpaces   
  
  
  
  
  
  
  function CorrectFilePath (input) RESULT (output)
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: outputTemp
    CHARACTER(LEN(input)) :: output
  
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, j, n

    ! -- Copy input string
    outputTemp = input
	
	DO i = 1, LEN( outputTemp )
	  output( i:i ) = " "
	End Do
	n = 0	
	DO i = 1, LEN( outputTemp )
	  if((iachar(outputTemp( i:i )).ge.32.and.iachar(outputTemp( i:i )).le.126))then
	    if(outputTemp( i:i ).ne."?".and.outputTemp( i:i ).ne."*".and.&
		outputTemp( i:i ).ne."%".and.outputTemp( i:i ).ne."+")then
	      n = n + 1
	      output( n:n ) = outputTemp( i:i )
		endif
	  endif
	End Do
	output = trim (output)
  
  END FUNCTION CorrectFilePath   
  
  
End Module strings