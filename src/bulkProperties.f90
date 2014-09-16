Module bulkProperties

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
  Use globals
  Use initialise
  Use loadData  
  Use readEAM
  Use calcEAM  
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: calcEquilibrium
  
Contains
  
  
  
  Subroutine calcEquilibrium()
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger) :: configID  
    
    Do configID=1,configCount
      
      If(mpiProcessID.eq.0)Then
        print *,configID,configRefEV(configID)
      End If
    
    End Do
    
  
    If(calcEqVol(1:3).eq."ALL")Then
    
    End If
  
  
  
  End Subroutine calcEquilibrium   
  
  

  
  
End Module bulkProperties  