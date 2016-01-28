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
    Integer(kind=StandardInteger) :: coordBP, configIDBP
    Integer(kind=StandardInteger) :: coordStartBP, coordEndBP, coordLengthBP
    Real(kind=DoubleReal), Dimension(1:2,1:3) :: bccUnit 
    Real(kind=DoubleReal), Dimension(1:4,1:3) :: fccUnit         
! -------------------
! Structure atom arrangements       
! -------------------       
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
! BCC atom 1
    bccUnit(1,1) = 0.0D0
    bccUnit(1,2) = 0.0D0
    bccUnit(1,3) = 0.0D0    
! BCC atom 1
    bccUnit(2,1) = 0.5D0
    bccUnit(2,2) = 0.5D0
    bccUnit(2,3) = 0.5D0     
! Init vars    
    configsAtomTotalBP = 0
    coordBP = 0
    coordStartBP = 1
! loop through bp configs
    Do configIDBP=1,maxConfigsBP
! break if no more configs stored
      If(bpInArr(configIDBP)%structure.eq."   ")Then
        Exit  
      End If
! store id into configCountBP
      configCountBP = configIDBP
! BP Copies           
      If(bpInArr(configIDBP)%structure.eq."FCC")Then    ! FCC
        Call makeConfigCoords(configIDBP, coordStartBP, coordEndBP, fccUnit)
        bpInArr(configIDBP)%atomsPerUnit = 4
      End If          
      If(bpInArr(configIDBP)%structure.eq."BCC")Then    ! BCC
        Call makeConfigCoords(configIDBP, coordStartBP, coordEndBP, bccUnit)
        bpInArr(configIDBP)%atomsPerUnit = 2
      End If
! store key data
      coordLengthBP = coordEndBP - coordStartBP + 1
      configurationCoordsKeyBP(configIDBP,1) = coordStartBP
      configurationCoordsKeyBP(configIDBP,2) = coordLengthBP
      configurationCoordsKeyBP(configIDBP,3) = coordEndBP
! increment starting coord id
      coordStartBP = coordEndBP + 1
    End Do
  End Subroutine makeConfigs
    
  Subroutine makeConfigCoords(configIDBP, coordStartBP, coordEndBP, unitCell)
    Implicit None   ! Force declaration of all variables
! Private variables  
    Integer(kind=StandardInteger) :: configIDBP
    Integer(kind=StandardInteger) :: coordStartBP, coordEndBP
    Integer(kind=StandardInteger) :: coordID, i
    Integer(kind=StandardInteger) :: xLoop, yLoop, zLoop
    Real(kind=DoubleReal), Dimension(:,:) :: unitCell
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: elementID
    Real(kind=DoubleReal) :: alat
    Integer(kind=StandardInteger) :: unitCopies
! Init
    coordID = coordStartBP
    coordEndBP = coordStartBP
    element = bpInArr(configIDBP)%element
    elementID = QueryUniqueElement(element)
    alat = bpInArr(configIDBP)%alat
    unitCopies = bpInArr(configIDBP)%size
! Store unit cell size    
    bpUnitCellCount(configIDBP) = unitCopies
! Only make if element is in EAM potential
    If(elementID.gt.0)Then
! Loop
      Do xLoop=1,unitCopies
        Do yLoop=1,unitCopies
          Do zLoop=1,unitCopies 
            Do i=1,size(unitCell,1)
! Element ID           
              configurationCoordsIBP(coordID,1) = elementID
! Fractional coords
              configurationCoordsRBP(coordID,4) = &
              (xLoop + unitCell(i,1) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRBP(coordID,5) = &
              (yLoop + unitCell(i,2) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRBP(coordID,6) = &
              (zLoop + unitCell(i,3) - 1.0D0)/(1.0D0*unitCopies)
! Real co-ordinates
              configurationCoordsRBP(coordID,1) = &
              configurationCoordsRBP(coordID,4) * (1.0D0*unitCopies) * alat
              configurationCoordsRBP(coordID,2) = &
              configurationCoordsRBP(coordID,5) * (1.0D0*unitCopies) * alat
              configurationCoordsRBP(coordID,3) = &
              configurationCoordsRBP(coordID,6) * (1.0D0*unitCopies) * alat
! Increment
              coordID = coordID + 1 
            End Do  
          End Do
        End Do   
      End Do
      coordEndBP = coordID - 1    
    End If
    
  End Subroutine makeConfigCoords
  

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!
  Function QueryUniqueElement (element) RESULT (output)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: output
    Integer(kind=StandardInteger) :: i, k, found
! Init    
    output = -1
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
