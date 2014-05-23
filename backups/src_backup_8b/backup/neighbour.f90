Module neighbour

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths
  Use initialise
  Use input


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKey
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListI
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtoms
  Real(kind=SingleReal), Dimension( : ), Allocatable :: neighbourListR
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: neighbourListCoords

!Privacy of functions/subroutines/variables
  Private
  Public :: runNeighbour	
  Public :: neighbourListKey	
  Public :: neighbourListI	
  Public :: neighbourListR		
  Public :: configAtoms			
  Public :: neighbourListCoords			
			    
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runNeighbour()
!force declaration of all variables
	Implicit None	
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k, l, m, n, x, y, z
	Integer(kind=StandardInteger) :: atoms
	Integer(kind=StandardInteger) :: atomA, atomB
	Integer(kind=StandardInteger) :: neighbourListLength, tempNeighbourListLength
	Integer(kind=StandardInteger) :: neighbourListCount, configStart, configLength
	Integer(kind=StandardInteger) :: coordsStart, coordsLength
	Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
	Real(kind=SingleReal) :: rCutoff, rCutoffNeg, rCutoffSq
	Real(kind=SingleReal) :: xcoord, ycoord, zcoord
	Real(kind=SingleReal) :: xshift, yshift, zshift
	Real(kind=SingleReal) :: xlat, ylat, zlat, alat
	Real(kind=SingleReal) :: xA, xB, yA, yB, zA, zB
	Real(kind=SingleReal) :: rdSq, xdSq, ydSq, zdSq
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKeyTemp
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListITemp
    Real(kind=SingleReal), Dimension( : ), Allocatable :: neighbourListRTemp
    Real(kind=SingleReal), Dimension( : , : ), Allocatable :: neighbourListCoordsTemp
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: coordsITemp
	Real(kind=SingleReal), Dimension( : , : ), Allocatable :: coordsRTemp
	
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	
!open coords file	
    if(saveFileCoords.eq."Y")then
	  outputFile = trim(currentWorkingDirectory)//"/"//"coords.dat"
	  open(unit=10,file=trim(outputFile))	
	endif
	
!open neighbour list file	
	!outputFile = trim(currentWorkingDirectory)//"/"//"neighboutlist.dat"
	!open(unit=11,file=trim(outputFile))
	
	
!write to output file
    write(999,"(A30,F8.4)") "Run Neighbour List Subroutine ",ProgramTime()
	
	!configHeaderR
	!1 LP, 2 RC
	!configHeaderI
	!1 X1, 2 X2, 3 X3, 4 Y1, 5 Y2, 6 Y3, 7 Z1, 8 Z2, 9 Z3, 10 CC1, 11 CC2, 12 CC3, 13 Start, 14 Length    
	!configCoordsI
	!configCoordsR

!Allocate temp neighbour list array - assume 100 neighbours per atom
	tempNeighbourListLength = 0
	do i=1,configCount
      tempNeighbourListLength = tempNeighbourListLength + &
	  configHeaderI(i,10) * configHeaderI(i,11)* configHeaderI(i,12) * &
	  configHeaderI(i,headerWidth) * 100
	enddo
	Allocate(neighbourListKeyTemp(1:size(configHeaderI)/headerWidth,1:2)) 
	Allocate(neighbourListITemp(1:tempNeighbourListLength,1:4))
	Allocate(neighbourListRTemp(1:tempNeighbourListLength)) 
	Allocate(neighbourListCoordsTemp(1:tempNeighbourListLength,1:6)) 
	Allocate(configAtoms(1:configCount)) 
	

!loop through configurations - estimate size of
	configStart = 1
	neighbourListCount = 0
	do i=1,configCount
!set variables      
      alat = configHeaderR(i,1)
	  xCopy = configHeaderI(i,10)
	  yCopy = configHeaderI(i,11)
	  zCopy = configHeaderI(i,12)
	
!expand unit cell	    
	  atoms = xCopy * yCopy * zCopy * configHeaderI(i,headerWidth)
	  configAtoms(i) = atoms
!Allocate temp array
	  if(Allocated(coordsITemp))then
	    Deallocate(coordsITemp)
	  endif
	  if(Allocated(coordsRTemp))then
	    Deallocate(coordsRTemp)
	  endif
	  Allocate(coordsITemp(1:atoms))
	  Allocate(coordsRTemp(1:atoms,1:3))
	  m = 0
!output to file	  
      write(999,"(A27,I4,F8.4)") "Make full cell for config: ",i,ProgramTime()	  
!loop over copies to make full crystal	  
	  coordsStart = configHeaderI(i,headerWidth-1)
	  coordsLength = configHeaderI(i,headerWidth)	  
	  do x=1,xCopy
	    do y=1,yCopy
	      do z=1,zCopy
		    do n=coordsStart,(coordsStart+coordsLength-1)  
			  m = m + 1
			  coordsITemp(m) = configCoordsI(n)
			  coordsRTemp(m,1) = alat * (x + configCoordsR(n,1) - 1)
			  coordsRTemp(m,2) = alat * (y + configCoordsR(n,2) - 1)
			  coordsRTemp(m,3) = alat * (z + configCoordsR(n,3) - 1)
              if(saveFileCoords.eq."Y")then
	            write(10,"(I8,I8,I8,F8.4,F8.4,F8.4)") i,m,coordsITemp(m),&
			    coordsRTemp(m,1),coordsRTemp(m,2),coordsRTemp(m,3)
			  endif
			enddo
		  enddo
		enddo
      enddo
!output to file	  
      write(999,"(A33,I4,F8.4)") "Build neighbour list for config: ",i,ProgramTime()	 
!build neighbour list
      rCutoff = configHeaderR(i,2)
      rCutoffNeg = -1 * configHeaderR(i,2)
      rCutoffSq = configHeaderR(i,2)**2
	  xlat = alat * xCopy
	  ylat = alat * yCopy
	  zlat = alat * zCopy
!loop through atom pairs    Atom A in inner cell, Atom B in 3x3x3 supercell      
	  configLength = 0
      do atomA=1,atoms
!loop through Atom B 3x3x3
	    do l=-1,1
	      do m=-1,1
	        do n=-1,1	              			
			  xshift = xlat * l
			  yshift = ylat * m
			  zshift = zlat * n
	          do atomB=1,atoms 			  
!don't self count atom
			    if(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)then
			    else
!check range one co-ord at a time then distance squared
                  xA = coordsRTemp(atomA,1)
                  xB = (xshift + coordsRTemp(atomB,1))
                  yA = coordsRTemp(atomA,2)
                  yB = (yshift + coordsRTemp(atomB,2))
                  zA = coordsRTemp(atomA,3)
                  zB = (zshift + coordsRTemp(atomB,3))
                  xdSq = (xA-xB)**2
				  if(xdSq.le.rCutoffSq)then
				    ydSq = (yA-yB)**2
				    if(ydSq.le.rCutoffSq)then
				      zdSq = (zA-zB)**2
				      if(zdSq.le.rCutoffSq)then
					    rdSq = xdSq + ydSq + zdSq
					    if(rdSq.le.rCutoffSq)then
					      neighbourListCount = neighbourListCount + 1
						  configLength = configLength + 1
!atom type data
						  neighbourListITemp(neighbourListCount,1) = coordsITemp(atomA) !Atom Type
						  neighbourListITemp(neighbourListCount,2) = coordsITemp(atomB)	!Atom Type
!atom id
						  neighbourListITemp(neighbourListCount,3) = atomA 	!Atom ID
						  neighbourListITemp(neighbourListCount,4) = atomB 	!Atom ID
!atom A and B co-ords (to recalculate rd if strain matrix applied)
						  neighbourListCoordsTemp(neighbourListCount,1) = xA	
						  neighbourListCoordsTemp(neighbourListCount,2) = yA
						  neighbourListCoordsTemp(neighbourListCount,3) = zA
						  neighbourListCoordsTemp(neighbourListCount,1) = xB
						  neighbourListCoordsTemp(neighbourListCount,2) = yB
						  neighbourListCoordsTemp(neighbourListCount,3) = zB
						  neighbourListRTemp(neighbourListCount) = rdSq**0.5
!Save neighbour list file
					    endif
					  endif
					endif
			      endif
				endif
			  enddo
			enddo
          enddo			
!end loop through Atom B 3x3x3	  
		enddo
	  enddo
!output to file	  
      write(999,"(A36,I4,F8.4)") "Neighbour list complete for config: ",i,ProgramTime()	 
      write(999,"(A28,I8,I8,F8.4)") "NL Start/Length for config: ",configStart,configLength,ProgramTime()	
!store configuration key information and update start point
      neighbourListKeyTemp(i,1) = configStart	  
      neighbourListKeyTemp(i,2) = configLength	  
	  configStart = configStart + configLength
	enddo
	
!output to file
    write(999,"(A35,F8.4)") "Move NL from temp array to global: ",ProgramTime()	 
!Move from temp arrays - Allocate arrays
	Allocate(neighbourListKey(1:size(configHeaderI)/headerWidth,1:2))
	Allocate(neighbourListI(1:tempNeighbourListLength,1:4))
	Allocate(neighbourListR(1:tempNeighbourListLength)) 
	Allocate(neighbourListCoords(1:tempNeighbourListLength,1:6)) 
!Transfer key data
    do i=1,size(configHeaderI)/headerWidth
	  neighbourListKey(i,1) = neighbourListKeyTemp(i,1)
	  neighbourListKey(i,2) = neighbourListKeyTemp(i,2)
	  !print *,i,neighbourListKey(i,1),neighbourListKey(i,2)
	enddo
!Transfer neighbour list data	
    do i=1,neighbourListCount
	  neighbourListI(i,1) = neighbourListITemp(i,1)
	  neighbourListI(i,2) = neighbourListITemp(i,2)
	  neighbourListI(i,3) = neighbourListITemp(i,3)
	  neighbourListI(i,4) = neighbourListITemp(i,4)
	  neighbourListR(i) = neighbourListRTemp(i)
	  neighbourListCoords(i,1) = neighbourListCoordsTemp(i,1)
	  neighbourListCoords(i,2) = neighbourListCoordsTemp(i,2)
	  neighbourListCoords(i,3) = neighbourListCoordsTemp(i,3)
	  neighbourListCoords(i,4) = neighbourListCoordsTemp(i,4)
	  neighbourListCoords(i,5) = neighbourListCoordsTemp(i,5)
	  neighbourListCoords(i,6) = neighbourListCoordsTemp(i,6)
	  !print *,i,neighbourListI(i,1),neighbourListI(i,2),&
	  !neighbourListI(i,3),neighbourListI(i,4),neighbourListR(i)	  
	enddo
!Deallocate arrays
    Deallocate(neighbourListITemp)
    Deallocate(neighbourListRTemp)
    Deallocate(neighbourListCoordsTemp)
	
!write to output file
    write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	write(999,"(A30,I8)") "Total pairs for all configs:  ",neighbourListCount
	
!close neighbour list file
    !close(11)	 
!close coords file
    if(saveFileCoords.eq."Y")then
      close(10)	 	 
    endif	  
!close output file
    close(999)	 

  End Subroutine runNeighbour
  
  
  
  Subroutine makeNeighbourList()
	
	!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
	
	
	
	
	
	
	
	
	

  End Subroutine makeNeighbourList
  
  

End Module neighbour