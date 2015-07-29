Module types
! --------------------------------------------------------------!
! Data types
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Define data types
! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------
! Setup Modules  
  Use kinds
! Force declaration of all variables
  Implicit None

! 
  
  Type :: chart
    Character(len=32) ::     title
    Character(len=32) ::     xAxis
    Character(len=32) ::     yAxis
    Real(kind=DoubleReal) :: xMin=1.1D99
    Real(kind=DoubleReal) :: xMax=-1.1D99
    Real(kind=DoubleReal) :: yMin=1.1D99
    Real(kind=DoubleReal) :: yMax=-1.1D99
    Logical ::               cleanPyFile=.true.
  End Type  

  Type :: bulkProperties
    Real(kind=DoubleReal) :: aLat=-2.1D20
    Real(kind=DoubleReal) :: v0=-2.1D20
    Real(kind=DoubleReal) :: e0=-2.1D20
    Real(kind=DoubleReal) :: b0=-2.1D20
    Real(kind=DoubleReal) :: bp0=-2.1D20
    Real(kind=DoubleReal) :: c11=-2.1D20
    Real(kind=DoubleReal) :: c12=-2.1D20
    Real(kind=DoubleReal) :: c44=-2.1D20
  End Type  

  Type :: rssConfig
    Real(kind=DoubleReal) :: total=0.0D0
    Real(kind=DoubleReal) :: energy=0.0D0
    Real(kind=DoubleReal) :: force=0.0D0
    Real(kind=DoubleReal) :: stress=0.0D0
  End Type 
  
  Type :: rssBP
    Real(kind=DoubleReal) :: total=0.0D0
    Real(kind=DoubleReal) :: aLat=0.0D0
    Real(kind=DoubleReal) :: v0=0.0D0
    Real(kind=DoubleReal) :: e0=0.0D0
    Real(kind=DoubleReal) :: b0=0.0D0
    Real(kind=DoubleReal) :: bp0=0.0D0
    Real(kind=DoubleReal) :: c11=0.0D0
    Real(kind=DoubleReal) :: c12=0.0D0
    Real(kind=DoubleReal) :: c44=0.0D0
  End Type  
  
  Type :: saConfig
    Real(kind=DoubleReal) ::         temp=10.0D0
    Integer(kind=StandardInteger) :: tempLoops=5
    Integer(kind=StandardInteger) :: varLoops=100
    Real(kind=DoubleReal) ::         maxVar=0.005D0
    Integer(kind=StandardInteger) :: refinementLoops=3
  End Type 
  
  
End Module types
