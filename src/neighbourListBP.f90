Module neighbourListBP

! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: neighbourListBP
! Updated: 19th June 2015
! --------------------------------------------------------------!
! Description:
! Loop through all configurations
! Build neighbour list for each configuration
! --------------------------------------------------------------!

! Setup Modules
  Use mpi
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
  Use geomTypes
  Use geom
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: makeNeighbourListBP
  Public :: applyDistortion
  Public :: bpStoreNL
  Public :: bpLoadNL

  Contains
  Subroutine makeNeighbourListBP()
! Make neighbour list of all atom pairs for all configs where separation is .le. rcutoff
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, coordStart, coordLength, coordEnd
    Integer(kind=StandardInteger) :: atomA, atomB, nlKey, asKey
    Integer(kind=StandardInteger) :: configStart, configEnd, configLength
    Integer(kind=StandardInteger) :: i, l, m, n
    Integer(kind=StandardInteger) :: unitCopies
    Real(kind=DoubleReal) :: rCutoff, rCutoffSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
    Real(kind=DoubleReal) :: rMin, rMax
    Real(kind=DoubleReal) :: startTime, endTime
    Type(nlType) :: nl
    Call initNL(nl)
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCountBP = 0
! Prepare NL key array - clear all above configIDStart
    Do configID=1,maxConfigsBP
      neighbourListKeyBP(configID,1) = 0
      neighbourListKeyBP(configID,2) = 0
      neighbourListKeyBP(configID,3) = 0
    End Do    
! Loop through configurations
    nlKey = 0
    configStart = 1
    Do configID=1,configCountBP
! Init config specific variables
      rMin = 2.0D21
      rMax = -2.0D21
! Start time
      startTime = MPI_Wtime()
! Check config is there
      If(configurationCoordsKeyBP(configID,1).gt.0)Then
! Config init values      
        coordStart = configurationCoordsKeyBP(configID,1)
        coordLength = configurationCoordsKeyBP(configID,2)
        coordEnd = configurationCoordsKeyBP(configID,3)
        aLat = (bpInArr(configID)%alat)*(bpInArr(configID)%size)  
        rCutoff = bpCutoffNL
        rCutoffSq = rCutoff**2
! make nl        
        Call makeNL(nl, configurationCoordsIBP, configurationCoordsRBP, coordStart, coordEnd, 6.5D0, aLat) 
        Do i=1,nl%length
          nlKey = nlKey + 1
! key/type
          neighbourListIBP(nlKey,1) = nl%i(i,3)  !Atom A type
          neighbourListIBP(nlKey,2) = nl%i(i,4)  !Atom B type
          neighbourListIBP(nlKey,3) = nl%i(i,1)  !Atom A id
          neighbourListIBP(nlKey,4) = nl%i(i,2)  !Atom B id
          neighbourListIBP(nlKey,5) = 0  !Atom A-B Key
          neighbourListIBP(nlKey,6) = nl%i(i,5)
! Store atom separation
          neighbourListRBP(nlKey) = nl%r(i,1)
! Atom coordinate data
          neighbourListCoordsBP(nlKey,1) = nl%r(i,5)  ! xA
          neighbourListCoordsBP(nlKey,2) = nl%r(i,6)
          neighbourListCoordsBP(nlKey,3) = nl%r(i,7)
          neighbourListCoordsBP(nlKey,4) = nl%r(i,8)
          neighbourListCoordsBP(nlKey,5) = nl%r(i,9)
          neighbourListCoordsBP(nlKey,6) = nl%r(i,10) ! zB
          neighbourListCoordsBP(nlKey,7) = 0.0D0
          neighbourListCoordsBP(nlKey,8) = 0.0D0
          neighbourListCoordsBP(nlKey,9) = 0.0D0
          neighbourListCoordsBP(nlKey,10) = nl%r(i,2)
          neighbourListCoordsBP(nlKey,11) = nl%r(i,3)
          neighbourListCoordsBP(nlKey,12) = nl%r(i,4)
        End Do        
! Store nl key
        neighbourListKeyBP(configID,1) = configStart
        neighbourListKeyBP(configID,2) = nl%length
        neighbourListKeyBP(configID,3) = configStart+nl%length-1
! Store other data
        neighbourListKeyRBP(configID,1) = nl%rVerlet
        neighbourListKeyRBP(configID,2) = nl%rMin
        neighbourListKeyRBP(configID,3) = nl%rMax
! Increment configStart for next loop/config
        configStart = configStart + nl%length
      End If  
! End time
      endTime = MPI_Wtime()    
      Call timeAcc(nlTime,startTime,endTime)
    End Do  ! End loop configs
! Output to file
    If(mpiProcessID.eq.0.and.saveNLToFile)Then
! save to file
      Open(UNIT=1,FILE=Trim(outputDirectory)//"/"//"nlFileBP.dat",&
      status="old",position="append",action="write")
      Do configID=1,configCountBP
        write(1,"(A15,I8)") "Configuration: ",configID
        configStart = neighbourListKeyBP(configID,1)
        configEnd = neighbourListKeyBP(configID,3)
        Do i=configStart,configEnd
          write(1,"(I8,I8,I8,F14.7)") i,neighbourListIBP(i,3),neighbourListIBP(i,4),&
          neighbourListRBP(i)
        End Do
      End Do
    End If
! Synch MPI processes
    Call M_synchProcesses()
  End Subroutine makeNeighbourListBP
! ---------------------------------------------------------------------------------------------------
  Subroutine applyDistortion(configID, dMatrix)
! Distorts crystal, changes position of atoms for enregy calculation only
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, configStart, configEnd
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: dMatrix
    Real(kind=DoubleReal), Dimension(1:3) :: xVect, xVectP
    Real(kind=DoubleReal), Dimension(1:3) :: yVect, yVectP
! Prepare vars
    configStart = neighbourListKeyBP(configID,1)
    configEnd = neighbourListKeyBP(configID,3)    
! Needs isometric check      
! Loop through and adjust values
    Do i=configStart,configEnd
      xVect(1) = neighbourListCoordsBP(i,1)
      xVect(2) = neighbourListCoordsBP(i,2)
      xVect(3) = neighbourListCoordsBP(i,3)
      yVect(1) = neighbourListCoordsBP(i,4)
      yVect(2) = neighbourListCoordsBP(i,5)
      yVect(3) = neighbourListCoordsBP(i,6)
! Transform coords
      xVectP = TransformCoords(xVect, dMatrix)
      yVectP = TransformCoords(yVect, dMatrix)
! Calculate rd      
      neighbourListRBP(i) = RdCoords (xVectP, yVectP)      
    End Do
  End Subroutine applyDistortion
    
  Subroutine bpStoreNL(configID)
! Stores coords and distance only (calculating energy, not forces)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, configStart, configEnd
    Integer(kind=StandardInteger) :: i, j
! Prepare vars
    i = 0
    j = 0
    configStart = neighbourListKeyBP(configID,1)
    configEnd = neighbourListKeyBP(configID,3)
! Store rij
    Do i=configStart,configEnd
      neighbourListRBP_T(i) = neighbourListRBP(i)
    End Do  
! Store coords (currently switched off
    !Do i=configStart,configEnd
      !Do j=1,6
      !  neighbourListCoordsBP_T(i,j) = neighbourListCoordsBP(i,j)
      !End Do  
    !End Do  
  End Subroutine bpStoreNL
  
  Subroutine bpLoadNL(configID)
! Stores coords and distance only (calculating energy, not forces)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, nlStart, nlEnd
    Integer(kind=StandardInteger) :: i, j
! Prepare vars
    j = 0
    nlStart = neighbourListKeyBP(configID,1)
    nlEnd = neighbourListKeyBP(configID,3)
! Store rij
    Do i=nlStart,nlEnd
      neighbourListRBP(i) = neighbourListRBP_T(i)
    End Do  
! Store coords (currently switched off
    !Do i=nlStart,nlEnd
      !Do j=1,6
      !  neighbourListCoordsBP(i,j) = neighbourListCoordsBP_T(i,j)
      !End Do  
    !End Do  
  End Subroutine bpLoadNL
  
End Module neighbourListBP
