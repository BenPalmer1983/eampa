Module neighbourList

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Read user input file 

!----------------------------------------
! Updated: 12th Aug 2014
!----------------------------------------

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
    Implicit None   ! Force declaration of all variables
! Private variables    
    Integer(kind=StandardInteger) :: configID, coordStart, coordLength, coordEnd
    Integer(kind=StandardInteger) :: atomA, atomB, nlKey, neighbourListCount, asKey
    Integer(kind=StandardInteger) :: configStart, configLength
    Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
    Integer(kind=StandardInteger) :: l, m, n
    Real(kind=DoubleReal) :: rCutoffSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCount = 0
    configStart = 1
! Loop through configurations     
    Do configID=1,configCount
! Init looping variables   
      atomA = 0
      atomB = 0
      nlKey = 0
      configLength = 0
      coordStart = configurationCoordsKeyG(configID,1)
      coordLength = configurationCoordsKeyG(configID,2)
      coordEnd = configurationCoordsKeyG(configID,3)
      rCutoffSq = configurationsR(configID,11)**2
      xCopy = configurationsI(configID,1) 
      yCopy = configurationsI(configID,2) 
      zCopy = configurationsI(configID,3) 
      aLat = configurationsR(configID,1)   
      nlUniqueKeys = 0    
!loop through Atom B 3x3x3
      Do l=-1,1
        Do m=-1,1
          Do n=-1,1   
!Set co-ordinate shift
            xShift = aLat * xCopy * l
            yShift = aLat * yCopy * m
            zShift = aLat * zCopy * n
!Reset unique key list
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
!check range one co-ord at a time then distance squared
                    xA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,1)
                    xB = 1.0D0*(xshift + configurationCoordsRG(coordStart+atomB-1,1))
                    yA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,2)
                    yB = 1.0D0*(yshift + configurationCoordsRG(coordStart+atomB-1,2))
                    zA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,3)
                    zB = 1.0D0*(zshift + configurationCoordsRG(coordStart+atomB-1,3))
                    xdSq = (xA-xB)**2
                    If(xdSq.le.rCutoffSq)then
                      ydSq = (yA-yB)**2
                      If(ydSq.le.rCutoffSq)then
                        zdSq = (zA-zB)**2
                        If(zdSq.le.rCutoffSq)then
                          rdSq = xdSq + ydSq + zdSq          
                          If(rdSq.le.rCutoffSq)then
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
                            atomSeparationSpread(asKey) = atomSeparationSpread(asKey) + 1
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
      configStart = configStart + configLength
    End Do  ! End loop configs
! Save tally to file    
    Open(UNIT=1,FILE=Trim(outputDirectory)//"/"//"nlSeparation.dat",&
    status="old",position="append",action="write") 
    Do n=1,size(atomSeparationSpread,1)
      Write(1,"(I8,A1,F6.3,A1,I8)") n," ",1.0D0*(n/100.0D0)," ",atomSeparationSpread(n)
    End Do
    Close(1)
! End time
    Call cpu_time(timeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("Neighbour List",timeEnd-timeStart)
  End Subroutine makeNeighbourList 

  
  
  
End Module neighbourList  