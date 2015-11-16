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
    Real(kind=DoubleReal) :: rCutoffSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
    Real(kind=DoubleReal) :: rMin, rMax
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCountRelax = 0
    configStart = 1
! Prepare NL key array - clear all above configIDStart
    Do i=1,maxConfigsRelax
      neighbourListKeyRelax(i,1) = 0
      neighbourListKeyRelax(i,2) = 0
      neighbourListKeyRelax(i,3) = 0
    End Do
! print *,"NL",configStart
! Loop through configurations
    Do configID=1,configCountRelax
! Init config specific variables
      rMin = 2.0D21
      rMax = -2.0D21
! Check config is there
      If(configurationCoordsKeyRelax(configID,1).gt.0)Then
! Init looping variables
        atomA = 0
        atomB = 0
        nlKey = 0
        configLength = 0
        coordStart = configurationCoordsKeyRelax(configID,1)
        coordLength = configurationCoordsKeyRelax(configID,2)
        coordEnd = configurationCoordsKeyRelax(configID,3)
! separation cutoff
        rCutoffSq = bpCutoffNL**2
        !print *,configID,coordStart,coordLength,coordEnd,aLatRelax(configID)
! loop through Atom B 3x3x3
        Do l=-1,1
          Do m=-1,1
            Do n=-1,1
! Set co-ordinate shift
              xShift = aLatRelax(configID) * l
              yShift = aLatRelax(configID) * m
              zShift = aLatRelax(configID) * n
! Reset unique key list
              nlUniqueKeysRelax = 0
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
                    If(nlUniqueKeysRelax(nlKey).eq.0)Then
                      nlUniqueKeysRelax(nlKey) = 1
                      xA = 1.0D0*configurationCoordsRRelax(coordStart+atomA-1,1)
                      xB = 1.0D0*(xshift + configurationCoordsRRelax(coordStart+atomB-1,1))
                      yA = 1.0D0*configurationCoordsRRelax(coordStart+atomA-1,2)
                      yB = 1.0D0*(yshift + configurationCoordsRRelax(coordStart+atomB-1,2))
                      zA = 1.0D0*configurationCoordsRRelax(coordStart+atomA-1,3)
                      zB = 1.0D0*(zshift + configurationCoordsRRelax(coordStart+atomB-1,3))
                      xdSq = (xA-xB)**2
                      If(xdSq.le.rCutoffSq)Then
                        ydSq = (yA-yB)**2
                        If(ydSq.le.rCutoffSq)Then
                          zdSq = (zA-zB)**2
                          If(zdSq.le.rCutoffSq)Then
                            rdSq = xdSq + ydSq + zdSq
                            If(rdSq.le.rCutoffSq)Then
                              neighbourListCountRelax = neighbourListCountRelax + 1
                              configLength = configLength + 1
! Store atom type/id data
                              neighbourListIRelax(neighbourListCountRelax,1) = &
                              configurationCoordsIRelax(coordStart+atomA-1,1)  !Atom A type
                              neighbourListIRelax(neighbourListCountRelax,2) = &
                              configurationCoordsIRelax(coordStart+atomB-1,1)  !Atom B type
                              neighbourListIRelax(neighbourListCountRelax,3) = atomA  !Atom A id
                              neighbourListIRelax(neighbourListCountRelax,4) = atomB  !Atom B id
                              neighbourListIRelax(neighbourListCountRelax,5) = nlKey  !Atom A-B Key
                              If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
                                neighbourListIRelax(neighbourListCountRelax,6) = 1
                              Else
                                neighbourListIRelax(neighbourListCountRelax,6) = 0
                              End If
                              neighbourListIRelax(neighbourListCountRelax,7) = coordStart+atomA-1
                              neighbourListIRelax(neighbourListCountRelax,8) = coordStart+atomB-1
! large cube coords for atom B
                              neighbourListIRelax(neighbourListCountRelax,9) = l
                              neighbourListIRelax(neighbourListCountRelax,10) = m
                              neighbourListIRelax(neighbourListCountRelax,11) = n
! Store atom separation
                              neighbourListRRelax(neighbourListCountRelax) = rdSq**0.5
! Atom coordinate data
                              neighbourListCoordsRelax(neighbourListCountRelax,1) = xA
                              neighbourListCoordsRelax(neighbourListCountRelax,2) = yA
                              neighbourListCoordsRelax(neighbourListCountRelax,3) = zA
                              neighbourListCoordsRelax(neighbourListCountRelax,4) = xB
                              neighbourListCoordsRelax(neighbourListCountRelax,5) = yB
                              neighbourListCoordsRelax(neighbourListCountRelax,6) = zB
                              neighbourListCoordsRelax(neighbourListCountRelax,7) = xA-xB
                              neighbourListCoordsRelax(neighbourListCountRelax,8) = yA-yB
                              neighbourListCoordsRelax(neighbourListCountRelax,9) = zA-zB
                              neighbourListCoordsRelax(neighbourListCountRelax,10) = &
                              (xA-xB)/neighbourListRRelax(neighbourListCountRelax)
                              neighbourListCoordsRelax(neighbourListCountRelax,11) = &
                              (yA-yB)/neighbourListRRelax(neighbourListCountRelax)
                              neighbourListCoordsRelax(neighbourListCountRelax,12) = &
                              (zA-zB)/neighbourListRRelax(neighbourListCountRelax)
! Tally atom separation
                              asKey = Ceiling(neighbourListRRelax(neighbourListCountRelax)*100)
                              If(asKey.lt.1)Then
                                asKey = 1
                              End If
                              !atomSeparationSpread(asKey) = atomSeparationSpread(asKey) + 1
                              If(neighbourListRRelax(neighbourListCountRelax).lt.rMin)Then
                                rMin = neighbourListRRelax(neighbourListCountRelax)
                              End If
                              If(neighbourListRRelax(neighbourListCountRelax).gt.rMax)Then
                                rMax = neighbourListRRelax(neighbourListCountRelax)
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
        neighbourListKeyRelax(configID,1) = configStart
        neighbourListKeyRelax(configID,2) = configLength
        neighbourListKeyRelax(configID,3) = configStart+configLength-1
        !print *,configID,configStart,configLength 
! Store other data
        neighbourListKeyRRelax(configID,1) = rCutoffSq**0.5D0
        neighbourListKeyRRelax(configID,2) = rMin
        neighbourListKeyRRelax(configID,3) = rMax
        
        neighbourListKeyRRelax(configID,6) = aLatRelax(configID)
! Increment configStart for next loop/config
        configStart = configStart + configLength
      End If
      Call cpu_time(timeEnd)
      !Call timeAcc(nlTimeBP,timeStart,timeEnd)
    End Do  ! End loop configs
! Output
    !If(mpiProcessID.eq.0)Then
! save to file
    !  Open(UNIT=1,FILE=Trim(outputDirectory)//"/"//"nlFileRelax.dat",&
    !  status="old",position="append",action="write")
      Do configID=1,configCountRelax
      
    !    write(1,"(A15,I8)") "Configuration: ",configID
    !    configStart = neighbourListKeyRelax(configID,1)
    !    configEnd = neighbourListKeyRelax(configID,3)
        
        Do i=configStart,configEnd
    !      write(1,"(I8,I8,I8,F14.7)") i,neighbourListIRelax(i,3),neighbourListIRelax(i,4),&
    !      neighbourListRRelax(i)
          !Print *,i,neighbourListIRelax(i,3),neighbourListIRelax(i,4),neighbourListRRelax(i)
        End Do
      End Do
    !End If
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
