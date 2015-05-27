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
    Real(kind=DoubleReal) :: rCutoffSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
    Real(kind=DoubleReal) :: rMin, rMax
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: crystalUnitCellTemp
    Real(kind=DoubleReal), Dimension(1:3) :: aVect, bVect
! Output
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      print *,"Make neighbour list"
    End If
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCount = 0
    configStart = 1
    rMin = 2.0D21
    rMax = -2.0D21
! Prepare NL key array - clear all above configIDStart
    Do i=1,maxConfigs
      neighbourListKey(i,1) = 0
      neighbourListKey(i,2) = 0
      neighbourListKey(i,3) = 0
    End Do
! print *,"NL",configStart
! Loop through configurations
    Do configID=1,configCount
! Check config is there
      If(configurationCoordsKeyG(configID,1).gt.0)Then
! Load transformation matrix      
        crystalUnitCellTemp(1,1) = crystalUnitCell(configID,1)
        crystalUnitCellTemp(1,2) = crystalUnitCell(configID,2)
        crystalUnitCellTemp(1,3) = crystalUnitCell(configID,3)
        crystalUnitCellTemp(2,1) = crystalUnitCell(configID,4)
        crystalUnitCellTemp(2,2) = crystalUnitCell(configID,5)
        crystalUnitCellTemp(2,3) = crystalUnitCell(configID,6)
        crystalUnitCellTemp(3,1) = crystalUnitCell(configID,7)
        crystalUnitCellTemp(3,2) = crystalUnitCell(configID,8)
        crystalUnitCellTemp(3,3) = crystalUnitCell(configID,9)
! Init looping variables
        atomA = 0
        atomB = 0
        nlKey = 0
        configLength = 0
        coordStart = configurationCoordsKeyG(configID,1)
        coordLength = configurationCoordsKeyG(configID,2)
        coordEnd = configurationCoordsKeyG(configID,3)
        !If(nlCutoff.lt.0.0D0)Then
          rCutoffSq = configurationsR(configID,11)**2
        !Else
        !  rCutoffSq = nlCutoff**2
        !End If
        xCopy = configurationsI(configID,1)
        yCopy = configurationsI(configID,2)
        zCopy = configurationsI(configID,3)
        aLat = configurationsR(configID,1)
! loop through Atom B 3x3x3
        Do l=-1,1
          Do m=-1,1
            Do n=-1,1      
! Set co-ordinate shift
              xShift = aLat * xCopy * l
              yShift = aLat * yCopy * m
              zShift = aLat * zCopy * n
! Reset unique key list
              nlUniqueKeys = 0
! Loop through atom pairs
              Do atomA=1,coordLength
                Do atomB=1,coordLength
                  If(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)Then  ! Don't self count atom
                  Else
                    If(atomA.lt.atomB)Then
                      nlKey = (atomB-1)*(atomB-2)/2+atomA
                    Else
                      nlKey = (atomA-1)*(atomA-2)/2+atomB
                    End If
                    If(nlUniqueKeys(nlKey).eq.0)Then
                      nlUniqueKeys(nlKey) = 1
! Atom A vector
                      aVect(1) = 1.0D0*configurationCoordsRG(coordStart+atomA-1,1)
                      aVect(2) = 1.0D0*configurationCoordsRG(coordStart+atomA-1,2)
                      aVect(3) = 1.0D0*configurationCoordsRG(coordStart+atomA-1,3)
                      aVect = TransformCoords (aVect, crystalUnitCellTemp)
! Atom A vector
                      bVect(1) = 1.0D0*(xshift + configurationCoordsRG(coordStart+atomB-1,1))
                      bVect(2) = 1.0D0*(yshift + configurationCoordsRG(coordStart+atomB-1,2))
                      bVect(3) = 1.0D0*(zshift + configurationCoordsRG(coordStart+atomB-1,3))
                      bVect = TransformCoords (bVect, crystalUnitCellTemp)
! Assign to variables
                      xA = aVect(1)
                      yA = aVect(2)
                      zA = aVect(3)
                      xB = bVect(1)
                      yB = bVect(2)
                      zB = bVect(3)
                      
                      xdSq = (xA-xB)**2
                      If(xdSq.le.rCutoffSq)Then
                        ydSq = (yA-yB)**2
                        If(ydSq.le.rCutoffSq)Then
                          zdSq = (zA-zB)**2
                          If(zdSq.le.rCutoffSq)Then
                            rdSq = xdSq + ydSq + zdSq
                            If(rdSq.le.rCutoffSq)Then
                              neighbourListCount = neighbourListCount + 1
                              configLength = configLength + 1
! Store atom type/id data
                              neighbourListI(neighbourListCount,1) = &
                              configurationCoordsIG(coordStart+atomA-1,1)  !Atom A type
                              neighbourListI(neighbourListCount,2) = &
                              configurationCoordsIG(coordStart+atomB-1,1)  !Atom B type
                              neighbourListI(neighbourListCount,3) = atomA  !Atom A id
                              neighbourListI(neighbourListCount,4) = atomB  !Atom B id
                              neighbourListI(neighbourListCount,5) = nlKey  !Atom A-B Key
                              If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
                                neighbourListI(neighbourListCount,6) = 1
                              Else
                                neighbourListI(neighbourListCount,6) = 0
                              End If
! Store atom separation
                              neighbourListR(neighbourListCount) = rdSq**0.5
! Atom coordinate data
                              neighbourListCoords(neighbourListCount,1) = xA
                              neighbourListCoords(neighbourListCount,2) = yA
                              neighbourListCoords(neighbourListCount,3) = zA
                              neighbourListCoords(neighbourListCount,4) = xB
                              neighbourListCoords(neighbourListCount,5) = yB
                              neighbourListCoords(neighbourListCount,6) = zB
                              neighbourListCoords(neighbourListCount,7) = xB-xA
                              neighbourListCoords(neighbourListCount,8) = yB-yA
                              neighbourListCoords(neighbourListCount,9) = zB-zA
                              neighbourListCoords(neighbourListCount,10) = &
                              (xB-xA)/neighbourListR(neighbourListCount)
                              neighbourListCoords(neighbourListCount,11) = &
                              (yB-yA)/neighbourListR(neighbourListCount)
                              neighbourListCoords(neighbourListCount,12) = &
                              (zB-zA)/neighbourListR(neighbourListCount)
! Tally atom separation
                              asKey = Ceiling(neighbourListR(neighbourListCount)*100)
                              If(asKey.lt.1)Then
                                asKey = 1
                              End If
                              atomSeparationSpread(asKey) = atomSeparationSpread(asKey) + 1
                              If(neighbourListR(neighbourListCount).lt.rMin)Then
                                rMin = neighbourListR(neighbourListCount)
                              End If
                              If(neighbourListR(neighbourListCount).gt.rMax)Then
                                rMax = neighbourListR(neighbourListCount)
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
        neighbourListKey(configID,1) = configStart
        neighbourListKey(configID,2) = configLength
        neighbourListKey(configID,3) = configStart+configLength-1
        neighbourListKeyR(configID,1) = rCutoffSq**0.5D0
        If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
          print "(A3,I4,I8,I8,A3,I8,A1,A1,F10.5,A1)",&
          "NL ",configID,configStart,(configStart+configLength-1),"  (",&
          configLength,")","[",(rCutoffSq**0.5D0),"]"
        End If
        configStart = configStart + configLength
      End If
      Call cpu_time(timeEnd)
      Call timeAcc(nlTime,timeStart,timeEnd)
    End Do  ! End loop configs
    If(mpiProcessID.eq.0)Then
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
