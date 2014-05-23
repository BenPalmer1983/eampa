Module constants

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None

!declare variables
  Real(kind=SingleReal), Parameter :: avogadrosConstant = 6.02214129E23 
  Real(kind=SingleReal), Parameter :: elementaryCharge = 1.60217656E-19 
  Real(kind=SingleReal), Parameter :: pi = 3.1415926535898 

  
!Privacy of functions/subroutines/variables
  Private
  Public :: avogadrosConstant				!Variable
  Public :: elementaryCharge				!Variable


End Module constants