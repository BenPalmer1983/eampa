Module units

! Setup Modules
  Use kinds
  Use strings

!force declaration of all variables
  Implicit None

!declare variables

  
!Privacy of functions/subroutines/variables
  Private
  Public :: NormaliseTime					!Function
  Public :: NormaliseLength					!Function
  Public :: NormaliseDensity				!Function
   
     
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
    
function NormaliseTime (inputTime, timeUnit) RESULT (outputTime)
!Convert times into seconds
    ! -- Argument and result
	Real(kind=DoubleReal) :: inputTime, factor
    CHARACTER(len=2) :: timeUnit
	Real(kind=DoubleReal):: outputTime
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, n
  
	timeUnit = StrToUpper(timeUnit)
	
	if(timeUnit(1:2).eq."HR")then
	  factor = 3600
	elseif(timeUnit(1:2).eq."MS")then
	  factor = 0.001
	elseif(timeUnit(1:1).eq."M")then
	  factor = 60
	elseif(timeUnit(1:1).eq."S")then
	  factor = 1.0
	elseif(timeUnit(1:1).eq."D")then
	  factor = 86400
	endif
  
    outputTime = factor * inputTime
  
End Function NormaliseTime  
  
  
function NormaliseLength (inputLength, lengthUnit) RESULT (outputLength)
!Convert lengths into metres
    ! -- Argument and result
	Real(kind=SingleReal) :: inputLength, factor
    CHARACTER(len=2) :: lengthUnit
	Real(kind=DoubleReal) :: outputLength
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, n
  
	lengthUnit = StrToUpper(lengthUnit)
	
	if(lengthUnit(1:2).eq."MM")then
	  factor = 1.0E-03
	elseif(lengthUnit(1:2).eq."CM")then
	  factor = 1.0E-02
	elseif(lengthUnit(1:1).eq."M")then
	  factor = 1.0
	elseif(lengthUnit(1:1).eq."A")then
	  factor = 1.0E-10
	endif
  
    outputLength = factor * inputLength
  
End Function NormaliseLength    
  
  
  
function NormaliseDensity (inputDensity, densityUnit) RESULT (outputDensity)
    ! -- Argument and result
	Real(kind=SingleReal) :: inputDensity, factor
    CHARACTER(len=5) :: densityUnit
	Real(kind=DoubleReal) :: outputDensity
    ! -- Local variables
    Integer(kind=StandardInteger) :: i, n
  
	densityUnit = StrToUpper(densityUnit)
	
	if(densityUnit(1:4).eq."KGM3")then
	  factor = 1.0
	elseif(densityUnit(1:4).eq."GCM3")then
	  factor = 1000
	endif
  
    outputDensity = factor * inputDensity
  
End Function NormaliseDensity  
  
  
  
  
  
End Module units