Module clean

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Read user input file 

!----------------------------------------
! Updated: 12th Aug 2014
!----------------------------------------

! Setup Modules
  Use kinds
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use initialise
  Use loadData  
  Use globals
  Use output
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: runClean
  
Contains
  Subroutine runClean()
    Implicit None   ! Force declaration of all variables
! Private variables    
! Synchronise all processes
    Call M_synchProcesses()
! Clean up files
    Call cleanFiles()
    
    
! Testing subroutine for testing new code developments
    Call testing()
  End Subroutine runClean
  
  
  
  
  Subroutine cleanFiles()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i
    Character(len=512) :: testLine
! Clean files on root process only
    If(mpiProcessID.eq.0)Then
      Call outputCleanupList()
      Do i=1,100
        testLine = fileCleanupList(i)
        If(testLine(1:1).eq." ")Then       
          Exit
        Else  
          Call rmFile(testLine)
        End If
      End Do     
    End If
  End Subroutine cleanFiles
  
  Subroutine testing()
    Implicit None   ! Force declaration of all variables
! Private variables    
    !Integer(kind=StandardInteger) :: i  
    Real(kind=DoubleReal), Dimension(1:4) :: pointA
    Real(kind=DoubleReal), Dimension(1:4) :: pointB
    Real(kind=DoubleReal), Dimension(1:6) :: coeffs
    
    pointA(1) = 0.2860000000D00 
    pointA(2) = 0.1000989695D03
    pointA(3) = -0.6178915350D03
    pointA(4) = 0.6903121141D04
    pointB(1) = 0.5590000000D00
    pointB(2) = 0.2880548468D02 
    pointB(3) = -0.9791722262D02  
    pointB(4) = 0.4997629802D03
    
    If(mpiProcessID.eq.0)Then
      coeffs = SplineAB(pointA, pointB)
    End If
    !print *,coeffs(1)
    !print *,coeffs(2)
    !print *,coeffs(3)
    !print *,coeffs(4)
    !print *,coeffs(5)
    !print *,coeffs(6)
      
  
  End Subroutine testing
  
  
  
End Module clean    
  