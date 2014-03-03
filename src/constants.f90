Module constants

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None

!declare variables
  Real(kind=SingleReal), Parameter :: avogadrosConstant = 6.02214129E23 
  Real(kind=SingleReal), Parameter :: elementaryCharge = 1.60217656E-19 
  Real(kind=SingleReal), Parameter :: pi = 3.1415926535898 
  Real(kind=DoubleReal), Parameter :: lnTwo = 0.693147180559945
  
!Privacy of functions/subroutines/variables
  Private
  Public :: avogadrosConstant				!Variable
  Public :: elementaryCharge				!Variable
  Public :: pi				                !Variable
  Public :: lnTwo				            !Variable


End Module constants