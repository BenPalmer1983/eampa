Module units

! Setup Modules
  Use kinds
  Use general

!force declaration of all variables
  Implicit None

!declare variables

  
!Privacy of functions/subroutines/variables
  Private
  Public :: NormaliseTime          !Function
  Public :: NormaliseLength          !Function
  Public :: NormaliseDensity        !Function
  Public :: UnitConvert                    !Function
   
     
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
    
  Function NormaliseTime (inputTime, timeUnit) RESULT (outputTime)
! Convert times into seconds
! Declare variables
    Real(kind=DoubleReal) :: inputTime, factor
    Character(len=2) :: timeUnit
    Real(kind=DoubleReal):: outputTime
! Initialise variables   
    factor = 1.0D0    
! Convert
    timeUnit = StrToUpper(timeUnit)  
    if(timeUnit(1:2).eq."HR")then
      factor = 3600D0
    elseif(timeUnit(1:2).eq."MS")then
      factor = 0.001D0
    elseif(timeUnit(1:1).eq."M")then
      factor = 60D0
    elseif(timeUnit(1:1).eq."S")then
      factor = 1.0D0
    elseif(timeUnit(1:1).eq."D")then
      factor = 86400D0
    endif
    outputTime = factor * inputTime  
 End Function NormaliseTime  
  
  
  Function NormaliseLength (inputLength, lengthUnit) RESULT (outputLength)
! Convert lengths into metres
! Declare variables
    Real(kind=DoubleReal) :: inputLength, factor
    CHARACTER(len=2) :: lengthUnit
    Real(kind=DoubleReal) :: outputLength
! Initialise variables   
    factor = 1.0D0  
! Convert    
    lengthUnit = StrToUpper(lengthUnit)  
    if(lengthUnit(1:2).eq."MM")then
      factor = 1.0D-03
    elseif(lengthUnit(1:2).eq."CM")then
      factor = 1.0D-02
    elseif(lengthUnit(1:1).eq."M")then
      factor = 1.0D0
    elseif(lengthUnit(1:1).eq."A")then
      factor = 1.0D-10
    endif  
    outputLength = factor * inputLength  
  End Function NormaliseLength    
  
  
  
  Function NormaliseDensity (inputDensity, densityUnit) RESULT (outputDensity)
! Normalise density input
    Real(kind=DoubleReal) :: inputDensity, factor
    CHARACTER(len=5) :: densityUnit
    Real(kind=DoubleReal) :: outputDensity
! Initialise variables   
    factor = 1.0D0      
! Normalise
    densityUnit = StrToUpper(densityUnit)  
    if(densityUnit(1:4).eq."KGM3")then
      factor = 1.0D0
    elseif(densityUnit(1:4).eq."GCM3")then
      factor = 1000D0
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
    outputValue = 0.0D0
    factorInput = 0.0D0
    factorOutput = 0.0D0
    inputUnit = "        "
    inputUnit = StrToUpper(inputUnitArg)
    outputUnit = "        "
    outputUnit = StrToUpper(outputUnitArg)  
!set factors
!-------------------------
! Default - if not found, assume default factors
!-------------------------
    factorInput = 0.0D0
    factorOutput = 0.0D0
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
! Volume
!-------------------------
!Volumes - input factor - convert to meters cubed
  If(inputUnit(1:4).eq."M3  ")Then
    factorInput = 1.0D0
  End If
    If(inputUnit(1:4).eq."ANG3")Then
    factorInput = 1.0D-30
  End If
!Volumes - output factor - convert from meters cubed
    If(outputUnit(1:4).eq."ANG3")Then
    factorOutput = 1.0D30
  End If 
  If(outputUnit(1:4).eq."M3  ")Then
    factorOutput = 1.0D0
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