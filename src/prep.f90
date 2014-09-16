Module prep

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
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: runPrep
  
Contains
  Subroutine runPrep()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Call setProcessMap()

  End Subroutine runPrep
  
  
  
  
  Subroutine setProcessMap()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i
    Do i=1,configCount
      processMap(i,1) = mod(i-1,mpiProcessCount)
    End Do
    
  End Subroutine setProcessMap
  
  
  
End Module prep    
  