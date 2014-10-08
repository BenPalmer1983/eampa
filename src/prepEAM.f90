Module prepEAM

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Prepare EAM file

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
  Use readEAM
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: runPrepEAM
  
Contains
  Subroutine runPrepEAM()
    Implicit None   ! Force declaration of all variables
! Private variables    

! Alloy EAM - make from input single element EAM potential
    If(eamMakeAlloy(1).ne."  ")Then
      print *,eamMakeAlloy(1),eamMakeAlloy(2),eamMakeAlloy(3)
      Call makeAlloyEAM()
    End If


  End Subroutine runPrepEAM
  
  
  Subroutine makeAlloyEAM()
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: i, j, k, m, n, functionCount
    Integer(kind=StandardInteger) :: startPoint, endPoint, dataLength
    Integer(kind=StandardInteger) :: pairComplete, densityComplete
    Integer(kind=StandardInteger) :: embeddingComplete, alloyElements
    Character(Len=64) :: eamAlloyFile
! Init    
    pairComplete = 0
    densityComplete = 0
    embeddingComplete = 0
    alloyElements = 0
    functionCount = 0
    eamAlloyFile = "eamAlloy.pot"
! Alloy elements    
    Do i=1,size(eamMakeAlloy,1)
      If(eamMakeAlloy(i).eq."  ")Then
        Exit
      End If
      alloyElements = alloyElements + 1
      Call AddUniqueElement(eamMakeAlloy(i))
    End Do    
! Make functions
    n = 0
    Do i=1,size(eamKey,1)
      If(eamKey(i,1).eq.-1)Then
        Exit ! break out
      End If
! Pair Functions
      If(eamKey(i,3).eq.1.and.pairComplete.eq.0)Then
        pairComplete = 1        
        Do j=1,alloyElements
          Do k=j,alloyElements
! Loop through data points
            functionCount = functionCount + 1
            startPoint = n + 1
            Do m=eamKey(i,4),eamKey(i,6)
              n = n + 1
              eamDataOpt(n,1) = eamData(m,1)
              eamDataOpt(n,2) = eamData(m,2)
              eamDataOpt(n,3) = eamData(m,3)
              eamDataOpt(n,4) = eamData(m,4)
            End Do
            endPoint = n
            dataLength = endPoint - startPoint + 1
            ! End looping through points
! Store key data
            eamKeyOpt(functionCount,1) = QueryUniqueElement(eamMakeAlloy(j))
            eamKeyOpt(functionCount,2) = QueryUniqueElement(eamMakeAlloy(k))
            eamKeyOpt(functionCount,3) = 1
            eamKeyOpt(functionCount,4) = startPoint
            eamKeyOpt(functionCount,5) = dataLength
            eamKeyOpt(functionCount,6) = endPoint
            If(mpiProcessID.eq.0)Then
              print *,functionCount,eamMakeAlloy(j),eamMakeAlloy(k)
            End If  
          End Do  
        End Do
      End If
! Density Functions
      If(eamKey(i,3).eq.2.and.densityComplete.eq.0)Then
        densityComplete = 1        
        Do j=1,alloyElements
! Loop through data points
          functionCount = functionCount + 1
          startPoint = n + 1
          Do m=eamKey(i,4),eamKey(i,6)
            n = n + 1
            eamDataOpt(n,1) = eamData(m,1)
            eamDataOpt(n,2) = eamData(m,2)
            eamDataOpt(n,3) = eamData(m,3)
            eamDataOpt(n,4) = eamData(m,4)
          End Do
          endPoint = n
          dataLength = endPoint - startPoint + 1
          ! End looping through points
! Store key data
          eamKeyOpt(functionCount,1) = QueryUniqueElement(eamMakeAlloy(j))
          eamKeyOpt(functionCount,2) = 0
          eamKeyOpt(functionCount,3) = 2
          eamKeyOpt(functionCount,4) = startPoint
          eamKeyOpt(functionCount,5) = dataLength
          eamKeyOpt(functionCount,6) = endPoint
            If(mpiProcessID.eq.0)Then
              print *,functionCount,eamMakeAlloy(j),eamMakeAlloy(k)
            End If  
        End Do
      End If     
! Embedding Functions
      If(eamKey(i,3).eq.3.and.embeddingComplete.eq.0)Then
        embeddingComplete = 1        
        Do j=1,alloyElements
! Loop through data points
          functionCount = functionCount + 1
          startPoint = n + 1
          Do m=eamKey(i,4),eamKey(i,6)
            n = n + 1
            eamDataOpt(n,1) = eamData(m,1)
            eamDataOpt(n,2) = eamData(m,2)
            eamDataOpt(n,3) = eamData(m,3)
            eamDataOpt(n,4) = eamData(m,4)
          End Do
          endPoint = n
          dataLength = endPoint - startPoint + 1
          ! End looping through points
! Store key data
          eamKeyOpt(functionCount,1) = QueryUniqueElement(eamMakeAlloy(j))
          eamKeyOpt(functionCount,2) = 0
          eamKeyOpt(functionCount,3) = 3
          eamKeyOpt(functionCount,4) = startPoint
          eamKeyOpt(functionCount,5) = dataLength
          eamKeyOpt(functionCount,6) = endPoint
            If(mpiProcessID.eq.0)Then
              print *,functionCount,eamMakeAlloy(j),eamMakeAlloy(k)
            End If  
        End Do
      End If   
    End Do
    eamFunctionCount = functionCount
! Store into original arrays and clear opt array
    eamKey = eamKeyOpt
    eamData = eamDataOpt
    eamKeyOpt = -1
    eamDataOpt = 0.0D0
    Call saveEamFile(eamAlloyFile)
  
  End Subroutine makeAlloyEAM
  
  
  
  
End Module prepEAM    
  