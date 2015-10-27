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
    Character(len=128) ::    tempDirectory 
    Character(len=128) ::    outputDirectory 
    Character(len=64) ::     outputName 
    Character(len=32) ::     title
    Character(len=32) ::     xAxis
    Character(len=32) ::     yAxis
    Real(kind=DoubleReal) :: xMin=1.1D99
    Real(kind=DoubleReal) :: xMax=-1.1D99
    Real(kind=DoubleReal) :: yMin=1.1D99
    Real(kind=DoubleReal) :: yMax=-1.1D99
    Logical ::               cleanPyFile=.true.
    Integer(kind=StandardInteger), Dimension(1:10) :: rowStart = -1
    Integer(kind=StandardInteger), Dimension(1:10) :: rowEnd = -1
    Integer(kind=StandardInteger), Dimension(1:10) :: colX = 1
    Integer(kind=StandardInteger), Dimension(1:10) :: colY = 2
  End Type  
  
  Type :: bpIn
    Character(len=3) :: structure="   "
    Integer(kind=StandardInteger) :: size=2          ! no unit cells
    Integer(kind=StandardInteger) :: atomsPerUnit=2  ! atoms per unit cell
    Character(len=2) :: element="  "
    Real(kind=DoubleReal) :: aLat=-2.1D20
    Real(kind=DoubleReal) :: v0=-2.1D20
    Real(kind=DoubleReal) :: e0=-2.1D20
    Real(kind=DoubleReal) :: b0=-2.1D20
    Real(kind=DoubleReal) :: bp0=-2.1D20
    Real(kind=DoubleReal) :: c11=-2.1D20
    Real(kind=DoubleReal) :: c12=-2.1D20
    Real(kind=DoubleReal) :: c44=-2.1D20
    Real(kind=DoubleReal) :: shearConstant=-2.1D20
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
    Real(kind=DoubleReal) :: shearConstant=-2.1D20
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
    Real(kind=DoubleReal) :: eos=0.0D0
    Real(kind=DoubleReal) :: c11=0.0D0
    Real(kind=DoubleReal) :: c12=0.0D0
    Real(kind=DoubleReal) :: c44=0.0D0
  End Type 
  
  Type :: eos
    Real(kind=DoubleReal) :: v0=0.0D0
    Real(kind=DoubleReal) :: e0=0.0D0
    Real(kind=DoubleReal) :: b0=0.0D0
    Real(kind=DoubleReal) :: bp0=0.0D0
  End Type   
  
  Type :: saConfig
    Real(kind=DoubleReal) ::         temp=10.0D0
    Real(kind=DoubleReal) ::         tempEnd=0.1D0
    Integer(kind=StandardInteger) :: tempLoops=5
    Integer(kind=StandardInteger) :: varLoops=100
    Real(kind=DoubleReal) ::         maxVar=0.05D0
    Real(kind=DoubleReal) ::         minVar=0.0005D0
    Integer(kind=StandardInteger) :: refinementLoops=3
  End Type 
  
  
End Module types
