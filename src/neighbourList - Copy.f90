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
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCount = 0
    configStart = 1
! Prepare NL key array - clear all above configIDStart
    Do i=1,maxConfigs
      neighbourListKey(i,1) = 0
      neighbourListKey(i,2) = 0
      neighbourListKey(i,3) = 0
    End Do
! print *,"NL",configStart
! Loop through configurations
    Do configID=1,configCount
! Init config specific variables
      rMin = 2.0D21
      rMax = -2.0D21
! Check config is there
      If(configurationCoordsKeyG(configID,1).gt.0)Then
! Init looping variables
        atomA = 0
        atomB = 0
        nlKey = 0
        configLength = 0
        coordStart = configurationCoordsKeyG(configID,1)
        coordLength = configurationCoordsKeyG(configID,2)
        coordEnd = configurationCoordsKeyG(configID,3)
        Do i=coordStart,coordEnd
          print *,i,configurationCoordsRG(i,1),configurationCoordsRG(i,2),configurationCoordsRG(i,3)
        End Do
! If(nlCutoff.lt.0.0D0)Then
        rCutoffSq = configurationsR(configID,11)**2
! Else
!  rCutoffSq = nlCutoff**2
! End If
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
! calculate the key of the atom A atom B combination
! half length list A=i,B=j == A=j,B=i
                    If(atomA.lt.atomB)Then
                      nlKey = (atomB-1)*(atomB)/2+atomA
                    Else
                      nlKey = (atomA-1)*(atomA)/2+atomB
                    End If
!
                    If(nlUniqueKeys(nlKey).eq.0)Then
                      nlUniqueKeys(nlKey) = 1
                      xA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,1)
                      xB = 1.0D0*(xshift + configurationCoordsRG(coordStart+atomB-1,1))
                      yA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,2)
                      yB = 1.0D0*(yshift + configurationCoordsRG(coordStart+atomB-1,2))
                      zA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,3)
                      zB = 1.0D0*(zshift + configurationCoordsRG(coordStart+atomB-1,3))
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
                              neighbourListCoords(neighbourListCount,7) = xA-xB
                              neighbourListCoords(neighbourListCount,8) = yA-yB
                              neighbourListCoords(neighbourListCount,9) = zA-zB
                              neighbourListCoords(neighbourListCount,10) = &
                              (xA-xB)/neighbourListR(neighbourListCount)
                              neighbourListCoords(neighbourListCount,11) = &
                              (yA-yB)/neighbourListR(neighbourListCount)
                              neighbourListCoords(neighbourListCount,12) = &
                              (zA-zB)/neighbourListR(neighbourListCount)
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
        print *,"NL Length: ",neighbourListCount," ",aLat," ",xCopy," ",rCutoffSq 
! Store nl key
        neighbourListKey(configID,1) = configStart
        neighbourListKey(configID,2) = configLength
        neighbourListKey(configID,3) = configStart+configLength-1
! Store other data
        neighbourListKeyR(configID,1) = rCutoffSq**0.5D0
        neighbourListKeyR(configID,2) = rMin
        neighbourListKeyR(configID,3) = rMax
! Increment configStart for next loop/config
        configStart = configStart + configLength
      End If
      Call cpu_time(timeEnd)
      Call timeAcc(nlTime,timeStart,timeEnd)
    End Do  ! End loop configs
! Output
! If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
!  print *,"----------------------------------------------------------------------"
! End If
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
