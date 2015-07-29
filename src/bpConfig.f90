Module bpConfig

! --------------------------------------------------------------!
! Bulk Property Configurations
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------

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
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: prepareBPConfig

  Contains
  Subroutine prepareBPConfig()
    Implicit None   ! Force declaration of all variables
! Private variables
! Print out
    If(TerminalPrint())Then
      print *,"Build Configurations for Bulk Property Testing"
    End If
    Call resetBPConfig() 
! Make
    Call makeConfigs()
! Output
    Call outputConfigBPPoints()
! Synch MPI processes
    Call M_synchProcesses()
  End Subroutine prepareBPConfig

  
  Subroutine resetBPConfig()  
    Implicit None   ! Force declaration of all variables
! Private variables
! Reset arrays
    configCountBP = 0
    configurationCoordsKeyBP = 0
    configurationCoordsIBP = 0
    configurationCoordsRBP = 0.0D0 
  End Subroutine resetBPConfig
  
  
  Subroutine makeConfigs()
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: coordBP, elementID, i, configIDBP
    Integer(kind=StandardInteger) :: xLoop, yLoop, zLoop
    Integer(kind=StandardInteger) :: coordStartBP, coordEndBP, coordLengthBP
    Real(kind=DoubleReal), Dimension(1:4,1:3) :: fccUnit 
        
! FCC atom 1
    fccUnit(1,1) = 0.0D0
    fccUnit(1,2) = 0.0D0
    fccUnit(1,3) = 0.0D0    
! FCC atom 1
    fccUnit(2,1) = 0.5D0
    fccUnit(2,2) = 0.5D0
    fccUnit(2,3) = 0.0D0    
! FCC atom 1
    fccUnit(3,1) = 0.5D0
    fccUnit(3,2) = 0.0D0
    fccUnit(3,3) = 0.5D0    
! FCC atom 1
    fccUnit(4,1) = 0.0D0
    fccUnit(4,2) = 0.5D0
    fccUnit(4,3) = 0.5D0    
    
    
! Counter    
    configIDBP = 0
! Build an FCC structure for each element
    configsAtomTotalBP = 0
    coordBP = 0
    coordStartBP = 1
    Do elementID=1,elementsCount!   
      configIDBP = configIDBP + 1    
      aLatBP(configIDBP) = fccAlatBP
      configCountBP = configCountBP + 1
      bpUnitCellCount(configIDBP) = size(fccUnit,1)
! BP Copies      
      Do xLoop=1,bpCopies
        Do yLoop=1,bpCopies
          Do zLoop=1,bpCopies
            Do i=1,size(fccUnit,1)              
              coordBP = coordBP + 1
              configurationCoordsIBP(coordBP,1) = elementID
! Fractional coords
              configurationCoordsRBP(coordBP,4) = &
              (xLoop + fccUnit(i,1) - 1.0D0)/(1.0D0*bpCopies)
              configurationCoordsRBP(coordBP,5) = &
              (yLoop + fccUnit(i,2) - 1.0D0)/(1.0D0*bpCopies)
              configurationCoordsRBP(coordBP,6) = &
              (zLoop + fccUnit(i,3) - 1.0D0)/(1.0D0*bpCopies)
! Real co-ordinates
              configurationCoordsRBP(coordBP,1) = &
              configurationCoordsRBP(coordBP,4) * (1.0D0*bpCopies) * aLatBP(configIDBP)
              configurationCoordsRBP(coordBP,2) = &
              configurationCoordsRBP(coordBP,5) * (1.0D0*bpCopies) * aLatBP(configIDBP)
              configurationCoordsRBP(coordBP,3) = &
              configurationCoordsRBP(coordBP,6) * (1.0D0*bpCopies) * aLatBP(configIDBP)
            End Do
          End Do
        End Do          
      End Do        
! Calculate keys
      coordEndBP = coordBP
      coordLengthBP = coordEndBP - coordStartBP + 1
      configsAtomTotalBP = configsAtomTotalBP + coordLengthBP
! Store keys  
      configurationCoordsKeyBP(elementID,1) = coordStartBP
      configurationCoordsKeyBP(elementID,2) = coordLengthBP
      configurationCoordsKeyBP(elementID,3) = coordEndBP
! reset start coord count
      coordStartBP = coordEndBP + 1      
    End Do
    
  
  End Subroutine makeConfigs
  

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!
  Function QueryUniqueElement (element) RESULT (output)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: output
    Integer(kind=StandardInteger) :: i, k, found
! convert to uppercase
    element = adjustl(StrToUpper(element))
! loop through elements array
    k = 0
    found = 0
    Do i=1,size(elements,1)
      k = k + 1
      If(elements(i).eq.element)Then
        found = 1
        exit
      End If
    End Do
! save element if not found
    If(found.eq.1)Then
      output = k
    Else
      output = 0
    End If
  End Function QueryUniqueElement

  Function QueryFunctionType (functionType) RESULT (output)
    Character(len=4) :: functionType
    Integer(kind=StandardInteger) :: output
    Integer(kind=StandardInteger) :: i
! convert to uppercase
    FunctionType = adjustl(StrToUpper(functionType))
    Do i=1,size(eamFunctionTypes,1)
      If(functionType.eq.eamFunctionTypes(i))Then
        output = i
        Exit
      End If
    End Do
  End Function QueryFunctionType

End Module bpConfig
