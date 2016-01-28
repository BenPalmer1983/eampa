Module neighbourList

! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! Module: neighbourList
! Updated: 18th May 2015
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
  Public :: makeNeighbourList

  Contains
  Subroutine makeNeighbourList()
! Make neighbour list of all atom pairs for all configs where separation is .le. rcutoff
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, coordStart, coordLength, coordEnd
    Integer(kind=StandardInteger) :: atomA, atomB, nlKey, neighbourListCount, asKey
    Integer(kind=StandardInteger) :: configStart, configEnd, configLength
    Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
    Integer(kind=StandardInteger) :: i, l, m, n
    Real(kind=DoubleReal) :: rCutoff, rCutoffSq
    Real(kind=DoubleReal) :: rVerlet, rVerletSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
    Real(kind=DoubleReal) :: rMin, rMax
    Real(kind=DoubleReal) :: startTime, endTime
    Type(nlType) :: nl
    Call initNL(nl)
! Init variables
    neighbourListCount = 0
    configStart = 1
! Prepare NL key array - clear all above configIDStart
    Do configID=1,maxConfigs
      neighbourListKey(configID,1) = 0
      neighbourListKey(configID,2) = 0
      neighbourListKey(configID,3) = 0
    End Do
! Loop through configurations
    nlKey = 0
    configStart = 1
    Do configID=1,configCount
! Init config specific variables
      rMin = 2.0D21
      rMax = -2.0D21
! Start time
      startTime = MPI_Wtime()
! Check config is there
      If(configurationCoordsKeyG(configID,1).gt.0)Then
! Config init values      
        coordStart = configurationCoordsKeyG(configID,1)
        coordLength = configurationCoordsKeyG(configID,2)
        coordEnd = configurationCoordsKeyG(configID,3)
        xCopy = configurationsI(configID,1)
        yCopy = configurationsI(configID,2)
        zCopy = configurationsI(configID,3)
        aLat = configurationsR(configID,1)        
        rCutoff = configurationsR(configID,11)
        rCutoffSq = rCutoff**2    
        rVerlet = configurationsR(configID,12)
        rVerletSq = rVerlet**2
! make nl        
        Call makeNL(nl, configurationCoordsIG, configurationCoordsRG, coordStart, coordEnd, rVerlet, aLat * xCopy) 
        Do i=1,nl%length
          nlKey = nlKey + 1
! key/type
          neighbourListI(nlKey,1) = nl%i(i,3)  !Atom A type
          neighbourListI(nlKey,2) = nl%i(i,4)  !Atom B type
          neighbourListI(nlKey,3) = nl%i(i,1)  !Atom A id
          neighbourListI(nlKey,4) = nl%i(i,2)  !Atom B id
          neighbourListI(nlKey,5) = 0  !Atom A-B Key
          neighbourListI(nlKey,6) = nl%i(i,5)
! Store atom separation
          neighbourListR(nlKey) = nl%r(i,1)
! Atom coordinate data
          neighbourListCoords(nlKey,1) = nl%r(i,5)  ! xA
          neighbourListCoords(nlKey,2) = nl%r(i,6)
          neighbourListCoords(nlKey,3) = nl%r(i,7)
          neighbourListCoords(nlKey,4) = nl%r(i,8)
          neighbourListCoords(nlKey,5) = nl%r(i,9)
          neighbourListCoords(nlKey,6) = nl%r(i,10) ! zB
          neighbourListCoords(nlKey,7) = 0.0D0
          neighbourListCoords(nlKey,8) = 0.0D0
          neighbourListCoords(nlKey,9) = 0.0D0
          neighbourListCoords(nlKey,10) = nl%r(i,2)
          neighbourListCoords(nlKey,11) = nl%r(i,3)
          neighbourListCoords(nlKey,12) = nl%r(i,4)
        End Do        
! Store nl key
        neighbourListKey(configID,1) = configStart
        neighbourListKey(configID,2) = nl%length
        neighbourListKey(configID,3) = configStart+nl%length-1
! Store other data
        neighbourListKeyR(configID,1) = rCutoff
        neighbourListKeyR(configID,2) = nl%rMin
        neighbourListKeyR(configID,3) = nl%rMax
        neighbourListKeyR(configID,6) = rVerlet
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
      Open(UNIT=1,FILE=Trim(outputDirectory)//"/"//"nlFile.dat",&
      status="old",position="append",action="write")
      Do configID=1,configCount
        write(1,"(A15,I8)") "Configuration: ",configID
        configStart = neighbourListKey(configID,1)
        configEnd = neighbourListKey(configID,3)
        Do i=configStart,configEnd
          write(1,"(I8,I8,I8,F14.7)") i,neighbourListI(i,3),neighbourListI(i,4),&
          neighbourListR(i)
        End Do
      End Do
    End If
! Synch MPI processes
    Call M_synchProcesses()
  End Subroutine makeNeighbourList

End Module neighbourList
