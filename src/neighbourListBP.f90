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
    Real(kind=DoubleReal) :: rCutoffSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
    Real(kind=DoubleReal) :: rMin, rMax
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCountBP = 0
    configStart = 1
! Prepare NL key array - clear all above configIDStart
    Do i=1,maxConfigsBP
      neighbourListKeyBP(i,1) = 0
      neighbourListKeyBP(i,2) = 0
      neighbourListKeyBP(i,3) = 0
    End Do
! print *,"NL",configStart
! Loop through configurations
    Do configID=1,configCountBP
! Init config specific variables
      rMin = 2.0D21
      rMax = -2.0D21
! Check config is there
      If(configurationCoordsKeyBP(configID,1).gt.0)Then
! Config data      
        aLat = bpInArr(configID)%alat
        unitCopies = bpInArr(configID)%size
! Init looping variables
        atomA = 0
        atomB = 0
        nlKey = 0
        configLength = 0
        coordStart = configurationCoordsKeyBP(configID,1)
        coordLength = configurationCoordsKeyBP(configID,2)
        coordEnd = configurationCoordsKeyBP(configID,3)
! separation cutoff
        rCutoffSq = bpCutoffNL**2
        !print *,aLat,coordStart,coordLength,coordEnd,unitCopies
! loop through Atom B 3x3x3
        Do l=-1,1
          Do m=-1,1
            Do n=-1,1
! Set co-ordinate shift
              xShift = aLat * unitCopies * l
              yShift = aLat * unitCopies * m
              zShift = aLat * unitCopies * n
! Reset unique key list
              nlUniqueKeysBP = 0
! Loop through atom pairs
              Do atomA=1,coordLength
                Do atomB=1,coordLength
                  If(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)Then  ! Don't self count atom
                  Else
! calculate the key of the atom A atom B combination
! half length list A=i,B=j == A=j,B=i
                    If(atomA.lt.atomB)Then
                      nlKey = (atomB-1)*(atomB-2)/2+atomA
                    Else
                      nlKey = (atomA-1)*(atomA-2)/2+atomB
                    End If
!
                    If(nlUniqueKeysBP(nlKey).eq.0)Then
                      nlUniqueKeysBP(nlKey) = 1
                      xA = 1.0D0*configurationCoordsRBP(coordStart+atomA-1,1)
                      xB = 1.0D0*(xshift + configurationCoordsRBP(coordStart+atomB-1,1))
                      yA = 1.0D0*configurationCoordsRBP(coordStart+atomA-1,2)
                      yB = 1.0D0*(yshift + configurationCoordsRBP(coordStart+atomB-1,2))
                      zA = 1.0D0*configurationCoordsRBP(coordStart+atomA-1,3)
                      zB = 1.0D0*(zshift + configurationCoordsRBP(coordStart+atomB-1,3))
                      xdSq = (xA-xB)**2
                      If(xdSq.le.rCutoffSq)Then
                        ydSq = (yA-yB)**2
                        If(ydSq.le.rCutoffSq)Then
                          zdSq = (zA-zB)**2
                          If(zdSq.le.rCutoffSq)Then
                            rdSq = xdSq + ydSq + zdSq
                            If(rdSq.le.rCutoffSq)Then
                              neighbourListCountBP = neighbourListCountBP + 1
                              configLength = configLength + 1
! Store atom type/id data
                              neighbourListIBP(neighbourListCountBP,1) = &
                              configurationCoordsIBP(coordStart+atomA-1,1)  !Atom A type
                              neighbourListIBP(neighbourListCountBP,2) = &
                              configurationCoordsIBP(coordStart+atomB-1,1)  !Atom B type
                              neighbourListIBP(neighbourListCountBP,3) = atomA  !Atom A id
                              neighbourListIBP(neighbourListCountBP,4) = atomB  !Atom B id
                              neighbourListIBP(neighbourListCountBP,5) = nlKey  !Atom A-B Key
                              If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
                                neighbourListIBP(neighbourListCountBP,6) = 1
                              Else
                                neighbourListIBP(neighbourListCountBP,6) = 0
                              End If
! Store atom separation
                              neighbourListRBP(neighbourListCountBP) = rdSq**0.5
! Atom coordinate data
                              neighbourListCoordsBP(neighbourListCountBP,1) = xA
                              neighbourListCoordsBP(neighbourListCountBP,2) = yA
                              neighbourListCoordsBP(neighbourListCountBP,3) = zA
                              neighbourListCoordsBP(neighbourListCountBP,4) = xB
                              neighbourListCoordsBP(neighbourListCountBP,5) = yB
                              neighbourListCoordsBP(neighbourListCountBP,6) = zB
                              neighbourListCoordsBP(neighbourListCountBP,7) = xA-xB
                              neighbourListCoordsBP(neighbourListCountBP,8) = yA-yB
                              neighbourListCoordsBP(neighbourListCountBP,9) = zA-zB
                              neighbourListCoordsBP(neighbourListCountBP,10) = &
                              (xA-xB)/neighbourListRBP(neighbourListCountBP)
                              neighbourListCoordsBP(neighbourListCountBP,11) = &
                              (yA-yB)/neighbourListRBP(neighbourListCountBP)
                              neighbourListCoordsBP(neighbourListCountBP,12) = &
                              (zA-zB)/neighbourListRBP(neighbourListCountBP)
! Tally atom separation
                              asKey = Ceiling(neighbourListRBP(neighbourListCountBP)*100)
                              If(asKey.lt.1)Then
                                asKey = 1
                              End If
                              !atomSeparationSpread(asKey) = atomSeparationSpread(asKey) + 1
                              If(neighbourListRBP(neighbourListCountBP).lt.rMin)Then
                                rMin = neighbourListRBP(neighbourListCountBP)
                              End If
                              If(neighbourListRBP(neighbourListCountBP).gt.rMax)Then
                                rMax = neighbourListRBP(neighbourListCountBP)
                              End If
                            End If
                          End If
                        End If
                      End If
                    End If
                  End If
                End Do
              End Do
            End Do
          End Do
        End Do
! Store nl key
        neighbourListKeyBP(configID,1) = configStart
        neighbourListKeyBP(configID,2) = configLength
        neighbourListKeyBP(configID,3) = configStart+configLength-1
! Store other data
        neighbourListKeyRBP(configID,1) = rCutoffSq**0.5D0
        neighbourListKeyRBP(configID,2) = rMin
        neighbourListKeyRBP(configID,3) = rMax
! Increment configStart for next loop/config
        configStart = configStart + configLength
      End If
      Call cpu_time(timeEnd)
      Call timeAcc(nlTimeBP,timeStart,timeEnd)
    End Do  ! End loop configs
! Output
    If(mpiProcessID.eq.0)Then
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
