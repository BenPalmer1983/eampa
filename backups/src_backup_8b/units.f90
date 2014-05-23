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
  Public :: UnitConvert    		            !Function
   
     
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
  
  
  
  
  
  Function UnitConvert(inputValue, inputUnitArg, outputUnitArg) RESULT (outputValue)
  !force declaration of all variables
	Implicit None
!declare private variables
    Character(*) :: inputUnitArg
    Character(*) :: outputUnitArg
    Character(len=8) :: inputUnit
    Character(len=8) :: outputUnit
	Real(kind=DoubleReal) :: inputValue
	Real(kind=DoubleReal) :: outputValue
	Real(kind=DoubleReal) :: factorInput, factorOutput
!Set variables
    inputUnit = "        "
    inputUnit = StrToUpper(inputUnitArg)
    outputUnit = "        "
	outputUnit = StrToUpper(outputUnitArg)
	
!set factors

!-------------------------
! Energy
!-------------------------
  
!Energies - input factor - convert to joules
    If(inputUnit(1:4).eq."RY  ")Then
	  factorInput = 2.179872172697D-18
	End If
    If(inputUnit(1:4).eq."EV  ")Then
	  factorInput = 1.602176568D-19
	End If
    If(inputUnit(1:4).eq."J   ")Then
	  factorInput = 1
	End If
!Energies - output factor - convert from joules	
    If(outputUnit(1:4).eq."RY  ")Then
	  factorOutput = 1.0D0/(2.179872172697D-18)
	End If 
    If(outputUnit(1:4).eq."EV  ")Then
	  factorOutput = 1.0D0/(1.602176568D-19)
	End If 
    If(outputUnit(1:4).eq."J   ")Then
	  factorOutput = 1.0D0
	End If 
 
 
!-------------------------
! Length
!-------------------------
 
!Lengths - input factor - convert to meters
    If(inputUnit(1:4).eq."M   ")Then
	  factorInput = 1.0D0
	End If
	If(inputUnit(1:4).eq."A   ".or.inputUnit(1:4).eq."ANGS")Then
	  factorInput = 1.0D-10
	End If
	If(inputUnit(1:4).eq."B   ".or.inputUnit(1:4).eq."BOHR")Then
	  factorInput = 5.291772109217D-11
	End If
!Lengths - output factor - convert from meters
    If(outputUnit(1:4).eq."M   ")Then
	  factorOutput = 1.0D0
	End If
	If(outputUnit(1:4).eq."A   ".or.outputUnit(1:4).eq."ANGS")Then
	  factorOutput = 1.0D10
	End If 
	If(outputUnit(1:4).eq."B   ".or.outputUnit(1:4).eq."BOHR")Then
	  factorOutput = 1.0D0/5.291772109217D-11
	End If
 
 
 
!-------------------------
! Stresses
!-------------------------

!Stresses - input factor - convert to Pa
    If(inputUnit(1:6).eq."RYBOH3")Then
	  factorInput = 1.471050658D13
	End If
    If(inputUnit(1:4).eq."PA  ")Then
	  factorInput = 1.0D0
	End If
    If(inputUnit(1:4).eq."GPA ")Then
	  factorInput = 1.0D9
	End If
    If(inputUnit(1:5).eq."EVAN3")Then
	  factorInput = 1.602176568D11
	End If
 !Stresses - output factor - convert from Pa
    If(outputUnit(1:6).eq."RYBOH3")Then
	  factorOutput = 1.0D0/1.471050658D13
	End If
    If(outputUnit(1:4).eq."PA  ")Then
	  factorOutput = 1.0D0
	End If
    If(outputUnit(1:4).eq."GPA ")Then
	  factorOutput = 1.0D-9
	End If
    If(outputUnit(1:5).eq."EVAN3")Then
	  factorOutput = 1.0D0/1.602176568D-11
	End If
 
 
 
  
!-------------------------
! Forces
!-------------------------

!Forces - input factor - convert to N
    If(inputUnit(1:6).eq."RYBOHR")Then
	  factorInput = 4.11936139295036D-8
	End If
    If(inputUnit(1:4).eq."N   ")Then
	  factorInput = 1.602176568D-9
	End If
    If(inputUnit(1:5).eq."EVANG")Then
	  factorInput = 1.602176568D-9
	End If
 !Stresses - output factor - convert from Pa
    If(outputUnit(1:6).eq."RYBOHR")Then
	  factorOutput = 1.0D0/4.11936139295036D-8
	End If
    If(outputUnit(1:4).eq."N   ")Then
	  factorOutput = 1.0D0
	End If
    If(outputUnit(1:5).eq."EVANG")Then
	  factorOutput = 1.0D0/1.602176568D-9
	End If
 
 
 
!-------------------------
! Output
!------------------------- 
 
    outputValue = inputValue * factorInput * factorOutput
  
  End Function UnitConvert  
  
  
  
  
End Module units