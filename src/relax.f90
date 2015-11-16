Module relax
! --------------------------------------------------------------!
! Relax
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Subroutines and functions to calculate energy, forces and stresses
! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------
! Setup Modules
  Use kinds
  Use types
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use plotTypes
  Use plot
  Use globals
  Use initialise
  Use loadData
  Use output
  Use eamGen
  Use neighbourListRelax
  Use relaxCalcEAM
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: relaxAtoms
  Public :: makeRelaxConfigs
  Contains
! ---------------------------------------------------------------------------------------------------  
  Subroutine resetRelaxConfig()  
    Implicit None   ! Force declaration of all variables
! Private variables
! Reset arrays
    configCountRelax = 0
    configurationCoordsKeyRelax = 0
    configurationCoordsIRelax = 0
    configurationCoordsRRelax = 0.0D0 
  End Subroutine resetRelaxConfig
  
  
! ---------------------------------------------------------------------------------------------------  

  Subroutine makeRelaxConfigs()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, configID    
    Integer(kind=StandardInteger) :: configStart, configLength, configEnd
    Integer(kind=StandardInteger) :: coordID, unitCopies, elementID
    Real(kind=DoubleReal) :: alat
    Integer(kind=StandardInteger) :: typeCalc

    Real(kind=DoubleReal), Dimension(1:2,1:3) :: bccUnit 
    Real(kind=DoubleReal), Dimension(1:4,1:3) :: fccUnit         
! ----------
! Init
! ----------
    elementID = 1
    unitCopies = 4
    alat = 4.04D0
    configurationCoordsKeyRelax = 0
    
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
    
    ! 1 = normal, 2 = missing 1 atom
    typeCalc = 2
    
    configID = 1
    Call makeRelaxCoords(configID, elementID, fccUnit, unitCopies, alat, configLength, typeCalc)
    configurationCoordsKeyRelax(1,1) = 1
    configurationCoordsKeyRelax(1,2) = configLength
    configurationCoordsKeyRelax(1,3) = configLength

    
    configCountRelax = 1    
    configsAtomTotalRelax = configLength
    
    
  End Subroutine makeRelaxConfigs
  
  Subroutine makeRelaxCoords(configID, elementID, unitCellIn, unitCopies, alat, configLength, typeCalc)
    Implicit None   ! Force declaration of all variables
! In/Out
    Integer(kind=StandardInteger) :: configID, elementID
    Real(kind=DoubleReal), Dimension(:,:) :: unitCellIn
    Integer(kind=StandardInteger) :: unitCopies
    Real(kind=DoubleReal) :: alat
    Integer(kind=StandardInteger) :: configLength
    Integer(kind=StandardInteger) :: typeCalc
! Private variables  
    Real(kind=DoubleReal), Dimension(1:size(unitCellIn,1),1:size(unitCellIn,2)) :: unitCell
    Real(kind=DoubleReal) :: randNumber, denominator
    Integer(kind=StandardInteger) :: coordID, i, j, atomCount
    Integer(kind=StandardInteger) :: xLoop, yLoop, zLoop
! Init
    coordID = 1
    unitCell = unitCellIn
! typeCalc
    If(typeCalc.eq.1)Then
      aLatRelax(1) = 1.0D0*unitCopies*alat
      ! Loop
      Do xLoop=1,unitCopies
        Do yLoop=1,unitCopies
          Do zLoop=1,unitCopies 
            Do i=1,size(unitCell,1)   
! Element ID           
              configurationCoordsIRelax(coordID,1) = elementID
! Fractional coords
              configurationCoordsRRelax(coordID,4) = &
              (xLoop + unitCell(i,1) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRRelax(coordID,5) = &
              (yLoop + unitCell(i,2) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRRelax(coordID,6) = &
              (zLoop + unitCell(i,3) - 1.0D0)/(1.0D0*unitCopies)
! Real co-ordinates
              configurationCoordsRRelax(coordID,1) = &
              configurationCoordsRRelax(coordID,4) * (1.0D0*unitCopies) * alat
              configurationCoordsRRelax(coordID,2) = &
              configurationCoordsRRelax(coordID,5) * (1.0D0*unitCopies) * alat
              configurationCoordsRRelax(coordID,3) = &
              configurationCoordsRRelax(coordID,6) * (1.0D0*unitCopies) * alat
! Increment
              coordID = coordID + 1 
            End Do  
          End Do
        End Do   
      End Do
      configLength = coordID - 1    
    End If
    
    
    
    If(typeCalc.eq.2)Then
      aLatRelax(1) = 1.0D0*unitCopies*alat
      ! Loop
      Do xLoop=1,unitCopies
        Do yLoop=1,unitCopies
          Do zLoop=1,unitCopies 
            Do i=1,size(unitCell,1)     
              If(coordID.eq.1)Then
                denominator = 5.0D0
              Else
                denominator = 100.0D0
              End If               
              Do j=1,3
                randNumber = RandomLCG()            
                unitCell(i,j) = (0.5D0-randNumber)/denominator + unitCellIn(i,j)
              End Do  
! Element ID           
              configurationCoordsIRelax(coordID,1) = elementID
! Fractional coords
              configurationCoordsRRelax(coordID,4) = &
              (xLoop + unitCell(i,1) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRRelax(coordID,5) = &
              (yLoop + unitCell(i,2) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRRelax(coordID,6) = &
              (zLoop + unitCell(i,3) - 1.0D0)/(1.0D0*unitCopies)
! Real co-ordinates
              configurationCoordsRRelax(coordID,1) = &
              configurationCoordsRRelax(coordID,4) * (1.0D0*unitCopies) * alat
              configurationCoordsRRelax(coordID,2) = &
              configurationCoordsRRelax(coordID,5) * (1.0D0*unitCopies) * alat
              configurationCoordsRRelax(coordID,3) = &
              configurationCoordsRRelax(coordID,6) * (1.0D0*unitCopies) * alat
! Increment
              coordID = coordID + 1 
            End Do  
          End Do
        End Do   
      End Do
      configLength = coordID - 1    
    End If
    
    
    
    If(typeCalc.eq.3)Then
      atomCount = size(unitCell,1)*unitCopies**3
      alat = alat*((atomCount-1.0D0)/(1.0D0*atomCount))**(1.0D0/3.0D0)
      aLatRelax(1) = 1.0D0*unitCopies*alat
! Loop
      Do xLoop=1,unitCopies
        Do yLoop=1,unitCopies
          Do zLoop=1,unitCopies 
            Do i=1,size(unitCell,1)   
              If(coordID.lt.256)Then
! Element ID           
              configurationCoordsIRelax(coordID,1) = elementID
! Fractional coords
              configurationCoordsRRelax(coordID,4) = &
              (xLoop + unitCell(i,1) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRRelax(coordID,5) = &
              (yLoop + unitCell(i,2) - 1.0D0)/(1.0D0*unitCopies)
              configurationCoordsRRelax(coordID,6) = &
              (zLoop + unitCell(i,3) - 1.0D0)/(1.0D0*unitCopies)
! Real co-ordinates
              configurationCoordsRRelax(coordID,1) = &
              configurationCoordsRRelax(coordID,4) * (1.0D0*unitCopies) * alat
              configurationCoordsRRelax(coordID,2) = &
              configurationCoordsRRelax(coordID,5) * (1.0D0*unitCopies) * alat
              configurationCoordsRRelax(coordID,3) = &
              configurationCoordsRRelax(coordID,6) * (1.0D0*unitCopies) * alat
! Increment
              coordID = coordID + 1 
              End If
            End Do  
          End Do
        End Do   
      End Do
      configLength = coordID - 1    
    End If
  End Subroutine makeRelaxCoords
  
  Subroutine relaxAtoms()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, k, n
    Real(kind=DoubleReal), Dimension(1:8192,1:6) :: configurationCoordsRRelax_T
    
! init
    mdCoordChange = 0.0D0
! make  
    Call makeRelaxConfigs()  
! store    
    Do i=1,configsAtomTotalRelax
      Do j=1,6
        configurationCoordsRRelax_T(i,j) = configurationCoordsRRelax(i,j)
      End Do
    End Do
    
    
    
    Do n = 1,10
      Call makeNeighbourListRelax()
      Call relaxCalcEnergies()      
      print *,n,configCalcEnergiesRelax(1)   
      Do i=1,configsAtomTotalRelax
        Do j=1,3
          configurationCoordsRRelax(i,j) = &
          configurationCoordsRRelax(i,j)+(0.5D0-RandomLCG())*0.1D0*configCalcForcesRelax(i,j)
        End Do
      End Do
    End Do
        
    
! load    
    Do i=1,configsAtomTotalRelax
      Do j=1,6
        configurationCoordsRRelax(i,j) = configurationCoordsRRelax_T(i,j)
      End Do
    End Do

    
    
    !Call updateNeighbourListRelax()
  
  End Subroutine relaxAtoms
  
  
  
  
  

  Subroutine relaxAtomsOld()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, configID, timeStep  
    Integer(kind=StandardInteger) :: coordStart, coordLength, coordEnd
    Real(kind=DoubleReal) :: totalForce, atomicMass
    Real(kind=DoubleReal) :: time, timeInc
    !Real(kind=DoubleReal), Dimension(1: :: 
! MD run
    print *,"Relax"
! init variables
    atomicMass = 106.0D0
    timeInc = 1.0D-5  ! in nanoseconds    
    mdVelocity = 0.0D0
    mdCoordChange = 0.0D0
! Make neighbour list
    Call makeNeighbourListRelax()
    
    Do timeStep = 1,50
      Call relaxCalcEnergies()
      print *,timeStep,configCalcEnergiesRelax(1)   
      Do i=1,configsAtomTotalRelax
! atom acceleration      
        mdAcceleration(i,1) = 9.6485D9*(configCalcForcesRelax(i,1)/atomicMass)
        mdAcceleration(i,2) = 9.6485D9*(configCalcForcesRelax(i,2)/atomicMass)
        mdAcceleration(i,3) = 9.6485D9*(configCalcForcesRelax(i,3)/atomicMass)
! atom velocity   
        mdVelocity(i,1) = mdVelocity(i,1) + mdAcceleration(i,1)*timeInc
        mdVelocity(i,2) = mdVelocity(i,2) + mdAcceleration(i,2)*timeInc
        mdVelocity(i,3) = mdVelocity(i,3) + mdAcceleration(i,3)*timeInc 
! Store coord-change        
        mdCoordChange(i,1) = 0.5*mdAcceleration(i,1)*&
        timeInc**2+mdVelocity(i,1)*timeInc
        mdCoordChange(i,2) = 0.5*mdAcceleration(i,2)*&
        timeInc**2+mdVelocity(i,2)*timeInc
        mdCoordChange(i,3) = 0.5*mdAcceleration(i,3)*&
        timeInc**2+mdVelocity(i,3)*timeInc
!         
        configurationCoordsRRelax(i,1) = configurationCoordsRRelax(i,1) + mdCoordChange(i,1)
        configurationCoordsRRelax(i,2) = configurationCoordsRRelax(i,2) + mdCoordChange(i,2)
        configurationCoordsRRelax(i,3) = configurationCoordsRRelax(i,3) + mdCoordChange(i,3)
        If(i.eq.1)Then
          print *,i,configCalcForcesRelax(i,1),configCalcForcesRelax(i,2),configCalcForcesRelax(i,3)
          print *,i,configurationCoordsRRelax(i,1),configurationCoordsRRelax(i,2),configurationCoordsRRelax(i,3)
        End If
      End Do
! Update neighbour list      
      Call updateNeighbourListRelax()
    End Do          
  End Subroutine relaxAtomsOld

End Module relax