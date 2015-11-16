Module neighbourListRelax

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
  Public :: makeNeighbourListRelax
  Public :: updateNeighbourListRelax

  Contains
  Subroutine makeNeighbourListRelax()
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
    neighbourListCountRelax = 0
    configStart = 1
! Prepare NL key array - clear all above configIDStart
    Do configID=1,maxConfigsRelax
      neighbourListKeyRelax(configID,1) = 0
      neighbourListKeyRelax(configID,2) = 0
      neighbourListKeyRelax(configID,3) = 0
    End Do
! Loop through configurations
    nlKey = 0
    configStart = 1
    Do configID=1,configCountRelax
! Init config specific variables
      rMin = 2.0D21
      rMax = -2.0D21
! Start time
      startTime = MPI_Wtime()
! Check config is there
      If(configurationCoordsKeyRelax(configID,1).gt.0)Then
! Config init values      
        coordStart = configurationCoordsKeyRelax(configID,1)
        coordLength = configurationCoordsKeyRelax(configID,2)
        coordEnd = configurationCoordsKeyRelax(configID,3)
        aLat = aLatRelax(configID)
        rCutoff = bpCutoffNL
        rCutoffSq = rCutoff**2
! make nl        
        Call makeNL(nl, configurationCoordsIRelax, configurationCoordsRRelax, coordStart, coordEnd, 6.5D0, aLat) 
        Do i=1,nl%length
          nlKey = nlKey + 1
! key/type
          neighbourListIRelax(nlKey,1) = nl%i(i,3)  !Atom A type
          neighbourListIRelax(nlKey,2) = nl%i(i,4)  !Atom B type
          neighbourListIRelax(nlKey,3) = nl%i(i,1)  !Atom A id
          neighbourListIRelax(nlKey,4) = nl%i(i,2)  !Atom B id
          neighbourListIRelax(nlKey,5) = 0  !Atom A-B Key
          neighbourListIRelax(nlKey,6) = nl%i(i,5)
! Store atom separation
          neighbourListRRelax(nlKey) = nl%r(i,1)
! Atom coordinate data
          neighbourListCoordsRelax(nlKey,1) = nl%r(i,5)  ! xA
          neighbourListCoordsRelax(nlKey,2) = nl%r(i,6)
          neighbourListCoordsRelax(nlKey,3) = nl%r(i,7)
          neighbourListCoordsRelax(nlKey,4) = nl%r(i,8)
          neighbourListCoordsRelax(nlKey,5) = nl%r(i,9)
          neighbourListCoordsRelax(nlKey,6) = nl%r(i,10) ! zB
          neighbourListCoordsRelax(nlKey,7) = 0.0D0
          neighbourListCoordsRelax(nlKey,8) = 0.0D0
          neighbourListCoordsRelax(nlKey,9) = 0.0D0
          neighbourListCoordsRelax(nlKey,10) = nl%r(i,2)
          neighbourListCoordsRelax(nlKey,11) = nl%r(i,3)
          neighbourListCoordsRelax(nlKey,12) = nl%r(i,4)
        End Do        
! Store nl key
        neighbourListKeyRelax(configID,1) = configStart
        neighbourListKeyRelax(configID,2) = nl%length
        neighbourListKeyRelax(configID,3) = configStart+nl%length-1
! Store other data
        neighbourListKeyRRelax(configID,1) = nl%rVerlet
        neighbourListKeyRRelax(configID,2) = nl%rMin
        neighbourListKeyRRelax(configID,3) = nl%rMax
! Increment configStart for next loop/config
        configStart = configStart + nl%length
      End If  
! End time
      endTime = MPI_Wtime()    
      Call timeAcc(nlTime,startTime,endTime)
    End Do  ! End loop configs
! Synch MPI processes
    Call M_synchProcesses()
  End Subroutine makeNeighbourListRelax
!
!
  Subroutine updateNeighbourListRelax()
! Make neighbour list of all atom pairs for all configs where separation is .le. rcutoff
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, nlKey
    Integer(kind=StandardInteger) :: nlKeyStart, nlKeyEnd
    Integer(kind=StandardInteger) :: atomAKey, atomBKey
    Real(kind=DoubleReal) :: rd, xd, yd, zd

            

    Do configID=1,configCountRelax

      nlKeyStart = neighbourListKeyRelax(configID,1)
      nlKeyEnd = neighbourListKeyRelax(configID,3)
      
      Do nlKey=nlKeyStart,nlKeyEnd
        atomAKey = neighbourListIRelax(nlKey,7)
        atomBKey = neighbourListIRelax(nlKey,8)        
        
! Update atom positions
        neighbourListCoordsRelax(nlKey,1) = neighbourListCoordsRelax(nlKey,1) + mdCoordChange(atomAKey,1) ! xA
        neighbourListCoordsRelax(nlKey,2) = neighbourListCoordsRelax(nlKey,2) + mdCoordChange(atomAKey,2) ! yA
        neighbourListCoordsRelax(nlKey,3) = neighbourListCoordsRelax(nlKey,3) + mdCoordChange(atomAKey,3) ! zA
        neighbourListCoordsRelax(nlKey,4) = neighbourListCoordsRelax(nlKey,4) + mdCoordChange(atomBKey,1) ! xB
        neighbourListCoordsRelax(nlKey,5) = neighbourListCoordsRelax(nlKey,5) + mdCoordChange(atomBKey,2) ! yB
        neighbourListCoordsRelax(nlKey,6) = neighbourListCoordsRelax(nlKey,6) + mdCoordChange(atomBKey,3) ! zB
! Calculate displacements
        xd = neighbourListCoordsRelax(nlKey,1) - neighbourListCoordsRelax(nlKey,4)
        yd = neighbourListCoordsRelax(nlKey,2) - neighbourListCoordsRelax(nlKey,5)
        zd = neighbourListCoordsRelax(nlKey,3) - neighbourListCoordsRelax(nlKey,6)
        rd = (xd**2+yd**2+zd**2)**0.5D0              
! Store
        neighbourListCoordsRelax(neighbourListCountRelax,7) = xd        
        neighbourListCoordsRelax(neighbourListCountRelax,8) = yd
        neighbourListCoordsRelax(neighbourListCountRelax,9) = zd
        If(atomAKey.eq.1)Then 
        !  print *,neighbourListRRelax(neighbourListCountRelax),rd
        End If 
        neighbourListRRelax(neighbourListCountRelax) = rd
        !mdCoordChange
        !print *,nlKey,atomAKey,atomBKey,neighbourListCoordsRelax(nlKey,1),neighbourListCoordsRelax(nlKey,4)
        If(nlKey.eq.nlKeyStart)Then
        !   print *,rd
        !  print *,nlKey,atomAKey,atomBKey,neighbourListCoordsRelax(nlKey,1),neighbourListCoordsRelax(nlKey,4),rd
        End If
        
        
      End Do  
      
    End Do


    

  End Subroutine updateNeighbourListRelax

    
  Subroutine relaxStoreNL(configID)
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
      neighbourListRRelax_T(i) = neighbourListRRelax(i)
    End Do  
! Store coords (currently switched off
    !Do i=configStart,configEnd
      !Do j=1,6
      !  neighbourListCoordsBP_T(i,j) = neighbourListCoordsBP(i,j)
      !End Do  
    !End Do  
  End Subroutine relaxStoreNL
  
  Subroutine relaxLoadNL(configID)
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
      neighbourListRRelax(i) = neighbourListRRelax_T(i)
    End Do  
! Store coords (currently switched off
    !Do i=nlStart,nlEnd
      !Do j=1,6
      !  neighbourListCoordsBP(i,j) = neighbourListCoordsBP_T(i,j)
      !End Do  
    !End Do  
  End Subroutine relaxLoadNL

End Module neighbourListRelax
