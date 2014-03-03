Module prep

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise
  Use input


!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!declare global variables  
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKey		!key for atom neighbour list
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListI          !neighbour list integer data
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtoms					!count of total atoms for each config
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtomsUnitCell			!
  Integer(kind=StandardInteger) :: configAtomsTotal
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListR
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: neighbourListCoords
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: configLatticeParameters
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationVolume
  
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationEnergy
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationBM
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptVolume
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptEnergy
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configurationForces
  

!Privacy of functions/subroutines/variables
  Private
  
!Variables - Define majority of global variables here
  Public :: runPrep	
  Public :: neighbourListKey	
  Public :: neighbourListI	
  Public :: neighbourListR		
  Public :: configAtoms	
  Public :: configAtomsUnitCell  
  Public :: neighbourListCoords		
  Public :: configLatticeParameters	
  Public :: configurationVolume
  Public :: configurationEnergy	
  Public :: configurationBM	
  Public :: configurationOptVolume	
  Public :: configurationOptEnergy	
  Public :: configAtomsTotal
  
!Subroutines  
  Public :: applyUnitVector
  
!Functions  
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runPrep()
!force declaration of all variables
	Implicit None	
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
!Make neighbour list
	Call makeNeighbourList()
	Call applyUnitVector(0)
	If(mpiProcessID.eq.0)Then
	  Call summariseConfigurations()
    End If	
	Call orderEamPotentials()
	Call calcEAMDerivative()
	If(mpiProcessID.eq.0)Then
	  Call storeEAM()
    End If
	Call synchMpiProcesses()
	Call prepDataStore()
	
	Call synchMpiProcesses()
  End Subroutine runPrep
  
  
!------------------------------------------------------------------------!
!                                                                        
! Neighbour List Subroutines                                              
!                                                                        
!------------------------------------------------------------------------!  
  
  Subroutine makeNeighbourList()
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
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListRTemp
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
!write to output file
    If(mpiProcessID.eq.0)Then
      write(999,"(F8.4,A2,A37)") ProgramTime(),"  ",&
	  "Start building config neighbour lists"
	End If
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
	Allocate(configAtomsUnitCell(1:configCount)) 
	Allocate(configLatticeParameters(1:configCount,1:3))

	
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
	  configAtomsUnitCell(i) = configHeaderI(i,headerWidth)
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
!build neighbour list
      rCutoff = configHeaderR(i,2)
      rCutoffNeg = -1 * configHeaderR(i,2)
      rCutoffSq = configHeaderR(i,2)**2
	  xlat = alat * xCopy
	  ylat = alat * yCopy
	  zlat = alat * zCopy
	  configLatticeParameters(i,1) = xlat
	  configLatticeParameters(i,2) = ylat
	  configLatticeParameters(i,3) = zlat
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
                  xA = 1.0D0*coordsRTemp(atomA,1)
                  xB = 1.0D0*(xshift + coordsRTemp(atomB,1))
                  yA = 1.0D0*coordsRTemp(atomA,2)
                  yB = 1.0D0*(yshift + coordsRTemp(atomB,2))
                  zA = 1.0D0*coordsRTemp(atomA,3)
                  zB = 1.0D0*(zshift + coordsRTemp(atomB,3))
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
						  neighbourListITemp(neighbourListCount,1) = coordsITemp(atomA) !Atom A Type
						  neighbourListITemp(neighbourListCount,2) = coordsITemp(atomB)	!Atom B Type
!atom id
						  neighbourListITemp(neighbourListCount,3) = atomA 	!Atom ID
						  neighbourListITemp(neighbourListCount,4) = atomB 	!Atom ID
!atom A and B co-ords (to recalculate rd if strain matrix applied)
						  neighbourListCoordsTemp(neighbourListCount,1) = xA	
						  neighbourListCoordsTemp(neighbourListCount,2) = yA
						  neighbourListCoordsTemp(neighbourListCount,3) = zA
						  neighbourListCoordsTemp(neighbourListCount,4) = xB
						  neighbourListCoordsTemp(neighbourListCount,5) = yB
						  neighbourListCoordsTemp(neighbourListCount,6) = zB
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
!store configuration key information and update start point
      neighbourListKeyTemp(i,1) = configStart	  
      neighbourListKeyTemp(i,2) = configLength	  
	  configStart = configStart + configLength
	enddo
!end loop through configurations


    
	
!Move from temp arrays - Allocate arrays
	Allocate(neighbourListKey(1:size(configHeaderI)/headerWidth,1:2))
	Allocate(neighbourListI(1:tempNeighbourListLength,1:4))
	Allocate(neighbourListR(1:tempNeighbourListLength)) 
	Allocate(neighbourListCoords(1:tempNeighbourListLength,1:9)) !Ax Ay Az, Bx By Bz, R_ABi R_ABj R_ABk 
!Transfer key data 
    !do i=1,size(configHeaderI)/headerWidth
    do i=1,configCount
	  neighbourListKey(i,1) = neighbourListKeyTemp(i,1)
	  neighbourListKey(i,2) = neighbourListKeyTemp(i,2)
	enddo
!Transfer neighbour list data	
    do i=1,neighbourListCount
	  neighbourListI(i,1) = neighbourListITemp(i,1)
	  neighbourListI(i,2) = neighbourListITemp(i,2)
	  neighbourListI(i,3) = neighbourListITemp(i,3)
	  neighbourListI(i,4) = neighbourListITemp(i,4)
	  neighbourListR(i) = neighbourListRTemp(i)
	  neighbourListCoords(i,1) = neighbourListCoordsTemp(i,1) !xA
	  neighbourListCoords(i,2) = neighbourListCoordsTemp(i,2) !yA
	  neighbourListCoords(i,3) = neighbourListCoordsTemp(i,3) !zA
	  neighbourListCoords(i,4) = neighbourListCoordsTemp(i,4) !xB
	  neighbourListCoords(i,5) = neighbourListCoordsTemp(i,5) !yB
	  neighbourListCoords(i,6) = neighbourListCoordsTemp(i,6) !zB 
	  neighbourListCoords(i,7) = neighbourListCoordsTemp(i,4)-neighbourListCoordsTemp(i,1)
	  neighbourListCoords(i,8) = neighbourListCoordsTemp(i,5)-neighbourListCoordsTemp(i,2)
	  neighbourListCoords(i,9) = neighbourListCoordsTemp(i,6)-neighbourListCoordsTemp(i,3)
	enddo
!Deallocate arrays
    Deallocate(neighbourListITemp)
    Deallocate(neighbourListRTemp)
    Deallocate(neighbourListCoordsTemp)

!close coords file
    if(saveFileCoords.eq."Y")then
      close(10)	 	 
    endif	
!output to file	
    If(mpiProcessID.eq.0)Then
      write(999,"(A6,A31,I8)") "      ","Total configurations prepared: "
	  write(999,"(F8.4,A2,A53,I8)") ProgramTime(),"  ",&
	  "Finish building config neighbour lists, total pairs: ",neighbourListCount
	End If
!close output file
    close(999)	 	

  End Subroutine makeNeighbourList
  
  
!------------------------------------------------------------------------!
!                                                                        
! Neighbour List and Configurations Prep                                              
!                                                                        
!------------------------------------------------------------------------!    
  
!------------------------------------------------------------------------!
! applyUnitVector 
!------------------------------------------------------------------------!
  
  Subroutine applyUnitVector(configurationID)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k				
!loop through configurations
    if(configurationID.eq.0)then
	  do i=1,configCount
!write to output file
        Call applyUnitVectorAction(i)
	  enddo
	else
	  Call applyUnitVectorAction(configurationID)
	endif	
  End Subroutine applyUnitVector
  
  
!------------------------------------------------------------------------!
! applyUnitVectorAction 
!------------------------------------------------------------------------!
      
  Subroutine applyUnitVectorAction(configurationID)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger) :: configStart, configLength
	Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, rD		
	configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)	
!Loop through neighbour list for the configuration
	do i=configStart,configStart+configLength-1 
!calculate atom positions with unitVector applied
	  xA = neighbourListCoords(i,1) * unitVector(1,1) + &
	  neighbourListCoords(i,2) * unitVector(2,1) + &
	  neighbourListCoords(i,3) * unitVector(3,1)
	  yA = neighbourListCoords(i,1) * unitVector(1,2) + &
	  neighbourListCoords(i,2) * unitVector(2,2) + &
	  neighbourListCoords(i,3) * unitVector(3,2)
	  zA = neighbourListCoords(i,1) * unitVector(1,3) + &
	  neighbourListCoords(i,2) * unitVector(2,3) + &
	  neighbourListCoords(i,3) * unitVector(3,3)
	  xB = neighbourListCoords(i,4) * unitVector(1,1) + &
	  neighbourListCoords(i,5) * unitVector(2,1) + &
	  neighbourListCoords(i,6) * unitVector(3,1)
	  yB = neighbourListCoords(i,4) * unitVector(1,2) + &
	  neighbourListCoords(i,5) * unitVector(2,2) + &
	  neighbourListCoords(i,6) * unitVector(3,2)
	  zB = neighbourListCoords(i,4) * unitVector(1,3) + &
	  neighbourListCoords(i,5) * unitVector(2,3) + &
	  neighbourListCoords(i,6) * unitVector(3,3)	  	  
!store new displacement vector between atom A and B
	  neighbourListCoords(i,7) = xB-xA
	  neighbourListCoords(i,8) = yB-yA
	  neighbourListCoords(i,9) = zB-zA
!store new distance between atoms		  
	  neighbourListR(i) = ((xA-xB)**2+(yA-yB)**2+(zA-zB)**2)**0.5
	enddo	
!recalculate configuration volume
    Call calcConfigurationVolume(configurationID)	
  End Subroutine applyUnitVectorAction
	
!------------------------------------------------------------------------!
! ***calcConfigurationVolume*** 
!------------------------------------------------------------------------!
   
  Subroutine calcConfigurationVolume(configurationID)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k
    Real(kind=SingleReal) :: xA, yA, zA
    Real(kind=SingleReal) :: xB, yB, zB
    Real(kind=SingleReal) :: xC, yC, zC
	Real(kind=SingleReal) :: crossProductI,crossProductJ,crossProductK
	Real(kind=SingleReal) :: volume
!if array not allocated, allocate it
    If(Allocated(configurationVolume))Then
	Else
	  Allocate(configurationVolume(1:configCount))
	End If  
!loop through configurations if 0 and calc volume for all configurations
    If(configurationID.eq.0)Then
	  Do i=1,configCount
	    xA = configLatticeParameters(i,1) * unitVector(1,1)
        xB = configLatticeParameters(i,1) * unitVector(1,2)
        xC = configLatticeParameters(i,1) * unitVector(1,3)
        yA = configLatticeParameters(i,2) * unitVector(2,1)
        yB = configLatticeParameters(i,2) * unitVector(2,2)
        yC = configLatticeParameters(i,2) * unitVector(2,3)
        zA = configLatticeParameters(i,3) * unitVector(3,1)
        zB = configLatticeParameters(i,3) * unitVector(3,2)
        zC = configLatticeParameters(i,3) * unitVector(3,3)
	    crossProductI = (xB*yC-xC*yB)
	    crossProductJ = (xC*yA-xA*yC)
	    crossProductK = (xA*yB-xB*yA)
	    volume = zA*crossProductI+zB*crossProductJ+zC*crossProductK	
        configurationVolume(i) = 1.0D0 * volume
	  End Do
	Else
      xA = configLatticeParameters(configurationID,1) * unitVector(1,1)
      xB = configLatticeParameters(configurationID,1) * unitVector(1,2)
      xC = configLatticeParameters(configurationID,1) * unitVector(1,3)
      yA = configLatticeParameters(configurationID,2) * unitVector(2,1)
      yB = configLatticeParameters(configurationID,2) * unitVector(2,2)
      yC = configLatticeParameters(configurationID,2) * unitVector(2,3)
      zA = configLatticeParameters(configurationID,3) * unitVector(3,1)
      zB = configLatticeParameters(configurationID,3) * unitVector(3,2)
      zC = configLatticeParameters(configurationID,3) * unitVector(3,3)
	  crossProductI = (xB*yC-xC*yB)
	  crossProductJ = (xC*yA-xA*yC)
	  crossProductK = (xA*yB-xB*yA)
	  volume = zA*crossProductI+zB*crossProductJ+zC*crossProductK	
      configurationVolume(configurationID) = 1.0D0 * volume
	End If
  
  End Subroutine calcConfigurationVolume	
  
  
  Subroutine summariseConfigurations()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i
!If master process
	If(mpiProcessID.eq.0)Then
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
!save summary of configurations to output file
      write(999,"(A6,A40)") "      ","----------------------------------------"
      write(999,"(A6,A25)") "      ","Summary of Configurations"
      write(999,"(A6,A40)") "      ","----------------------------------------"
      Do i=1,configCount
        write(999,"(A6,A15,I8,A12,F14.6,A11,I8,A15,I8)") "      ",&
	    "Configuration: ",i,"    Volume: ",configurationVolume(i),&
	    "    Atoms: ",configAtoms(i),"    NL Length: ",neighbourListKey(i,2)
	  End Do	
!close the output file
      close(999) 	
	End If
  End Subroutine summariseConfigurations
  
  
  
  
  
  
!Order eam potentials
  Subroutine orderEamPotentials()  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, sortLoop, sorted, reordered
    Integer(kind=StandardInteger) :: potStart, potLength
	Real(kind=SingleReal) :: tempXA, tempYA, tempXB, tempYB	
!write if master process
    If(mpiProcessID.eq.0)Then	
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
!write to output file
	  write(999,"(F8.4,A2,A20)") ProgramTime(),"  ","Order EAM Potentials"	
    End If
	reordered = 0
	do i=1,size(eamKey,1)
	  potStart = eamKey(i,4)
	  potLength = eamKey(i,5)
	  If(mpiProcessID.eq.0)Then
	    write(999,"(I8,I8)") potStart, potLength	
	  End If
	  sorted = 1
	  do while(sorted.eq.1)
	    sorted = 0
	    do j=potStart,potStart+potLength-2	    
	      if(eamData(j,1).gt.eamData(j+1,1))then
	        sorted = 1
			reordered = reordered + 1
!temporarily store variables
			tempXA = eamData(j,1)
			tempYA = eamData(j,2)
			tempXB = eamData(j+1,1)
			tempYB = eamData(j+1,2)
!move variables
			eamData(j,1) = tempXB
			eamData(j,2) = tempYB
			eamData(j+1,1) = tempXA
			eamData(j+1,2) = tempYA
	      endif
	    enddo
	  enddo		
	enddo	
	If(mpiProcessID.eq.0)Then
!close output file
      close(999)	
    End If	  
  End Subroutine orderEamPotentials

  
  
  
    
  Subroutine calcEAMDerivative()
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k  
	Integer(kind=StandardInteger) :: eamCount, eamDataCount
	Integer(kind=StandardInteger) :: eamDataStart, eamDataEnd
	Integer(kind=StandardInteger) :: interpPointCount, posOffset
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: interpPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: interpResults
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArraySend
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArrayRecv
	Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processID,processCount
	Integer(kind=StandardInteger) :: processFrom, processTo
	Real(kind=DoubleReal) :: testSend, testRecv
!calculate derivatives for pair, density and embedding functions

!PAIR 1, DENS/SDEN 2, EMBE/SEMB 3, DDEN 4, DEMB 5

!call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)

!Set variables
    interpPointCount = eamInterpType + 1 !set interpolation count
    eamCount = size(eamKey,1)
	eamDataCount = size(eamData,1)
!Allocate mpi buffer array
	Allocate(bufferArraySend(1:eamDataCount))
	Allocate(bufferArrayRecv(1:eamDataCount))
!Initialise arrays
    Do i=1,eamDataCount
	  bufferArraySend(i) = 0.0D0
	  bufferArrayRecv(i) = 0.0D0
	End Do
!Assign MPI processes 
    Do i=1,eamCount
	  selectProcess = mod(i-1,mpiProcessCount)
	  If(selectProcess.eq.mpiProcessID)Then
	    !print *,i,selectProcess,interpPointCount
!data position start and end points
        eamDataStart = eamKey(i,4)
		eamDataEnd = eamKey(i,4)+eamKey(i,5)-1
!Calculate derivative for function - loop through data points
        !print *,i,eamDataStart,eamDataEnd
		Do j=eamDataStart,eamDataEnd
!Allocate array if not allocated
          If(Allocated(interpPoints))Then
		  Else 
            Allocate(interpPoints(1:interpPointCount,1:2))
		  End If
!prepare set of data points used in interp
		  posOffset = -1*floor(interpPointCount/2.0D0)
		  If((j+posOffset).lt.eamDataStart)Then
		    posOffset = 0  
		  End If
		  If((j-posOffset).gt.eamDataEnd)Then
		    posOffset = -1*interpPointCount  
		  End If
!set x value		  
		  x = eamData(j,1)
!make interp array
		  Do k=1,interpPointCount
		    interpPoints(k,1) = eamData(j+posOffset+k,1)
		    interpPoints(k,2) = eamData(j+posOffset+k,2)
		  End Do
!interpolate and get derivative at the point
          interpResults = PointInterpolationFull(interpPoints,x)
		  If(mpiProcessID.eq.0)Then
		    eamData(j,3) = interpResults(2)  !store Master process directly
		  Else 
		    bufferArraySend(j) = interpResults(2) !store Worker process into buffer to send
		  End If	
		  !print *,selectProcess,j,x,eamData(j,2),interpResults(1),eamData(j,3),bufferArray(j)
		End Do		
	  End If
	End Do  


!Gather data from MPI worker processes
    If(mpiProcessID.gt.0) Then
	  testSend = 1.3D0
!send buffers from all worker processes
      processTo = 0
      Do i=1,(mpiProcessCount-1)
        If(i.eq.mpiProcessID)Then
          tag = 1100 + i
		  !print *,"send from ",mpiProcessID,eamDataCount,tag	
          Call MPI_send(bufferArraySend,eamDataCount,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
		  !Call MPI_send(testSend,1,&
		  !MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End If
      End Do
    End If
	If(mpiProcessID.eq.0) Then
!read buffers from worker processes
      Do i=1,(mpiProcessCount-1)
		processFrom = i		!set processFrom in MPI as a separate var for mpi call only
        tag = 1100 + i
	    !print *,"recv by ",mpiProcessID,eamDataCount,tag
        Call MPI_recv(bufferArrayRecv,eamDataCount,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
		Do j=1,eamCount
		  selectProcess = mod(j-1,mpiProcessCount)
		  !print *,i,selectProcess
	      If(selectProcess.eq.i)Then
		    eamDataStart = eamKey(j,4)
		    eamDataEnd = eamKey(j,4)+eamKey(j,5)-1
			!print *,"...",j,eamDataStart,eamDataEnd
		    Do k=eamDataStart,eamDataEnd
		      eamData(k,3) = bufferArrayRecv(k)
			End Do  
	      End If
	    End Do 
      End Do
	End If  
!Send out data to all processes
    
!store column in buffer - from main process only
	If(mpiProcessID.eq.0) Then
	  Do i=1,eamDataCount
	    bufferArraySend(i) = eamData(i,3)
	  End Do
	  Do i=1,(mpiProcessCount-1)
        tag = 1200 + i	
		processFrom = i
        Call MPI_send(bufferArraySend,eamDataCount,MPI_double_precision,processFrom,tag,&
		MPI_comm_world,error)
      End Do
	Else
!retrieve array
      tag = 1200 + processID
	  processTo = 0
	  Call MPI_recv(bufferArrayRecv,eamDataCount,MPI_double_precision,processTo,tag,&
      MPI_comm_world,status,error)  
	  Do i=1,eamDataCount
	    eamData(i,3) = bufferArrayRecv(i)
	  End Do
    End If
  
  End Subroutine calcEAMDerivative
  
  
  
  
  
  
  
  
  Subroutine storeEAM()
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k  
!store eam data and close file
    if(saveFilePot.eq."Y".and.mpiProcessID.eq.0)then
	  do i=1,size(eamKey,1)
	    if(eamKey(i,3).eq.1)then
	      write(21,"(A5,A2,A1,A2,A4,I4,I4)") "PAIR ",elements(eamKey(i,1)),&
		  " ",elements(eamKey(i,2)),"    ",eamKey(i,1),eamKey(i,2)
		endif  
		if(eamKey(i,3).eq.2)then
		  if(eamType.eq.1)then
	        write(21,"(A5,A2,A4,I4)") "DENS ",elements(eamKey(i,1)),"    ",&
			eamKey(i,1)
		  endif
		  if(eamType.eq.2)then
	        write(21,"(A5,A2,A1,A2,A4,I4,I4)") "DENS ",elements(eamKey(i,1))," ",&
			elements(eamKey(i,2)),"    ",eamKey(i,1),eamKey(i,2)
		  endif
		endif 
		if(eamKey(i,3).eq.3)then
	      write(21,"(A5,A2,A4,I4)") "EMBE ",elements(eamKey(i,1)),"    ",&
		  eamKey(i,1)
		endif 
		if(eamKey(i,4).eq.4)then
	      write(21,"(A5,A2,A4,I4)") "DEND ",elements(eamKey(i,1)),"    ",&
		  eamKey(i,1)
		endif 
		if(eamKey(i,4).eq.4)then
	      write(21,"(A5,A2,A4,I4)") "EMBD ",elements(eamKey(i,1)),"    ",&
		  eamKey(i,1)
		endif 		
		k = 1
		do j=eamKey(i,4),(eamKey(i,4)+eamKey(i,5)-1)
		  write(21,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") eamData(j,1),"  ",&
		  eamData(j,2),"  ",eamData(j,3),"    ",i,k,j
		  k = k + 1
		enddo
      enddo
      close(21)
    endif	
  
  End Subroutine storeEAM
  
  
  
  
  
  
  
  Subroutine prepDataStore()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i, totalAtoms   
!Configuration Energies	
!Allocate array
	If(Allocated(configurationEnergy))Then
	  Deallocate(configurationEnergy)
	End If	  
    Allocate(configurationEnergy(1:configCount))
!zero out energies
    Do i=1,configCount
      configurationEnergy(i) = 0.0D0	
	End Do
!count total atoms	
	configAtomsTotal = 0
	Do i=1,configCount
	  configAtomsTotal = configAtomsTotal + configAtoms(i)
	End Do
!configuration atom forces 
    If(calcForcesOnOff.eq.1)Then
	  If(Allocated(configurationForces))Then
	    Deallocate(configurationForces)
	  End If
      Allocate(configurationForces(1:configCount,1:configAtomsTotal))
    End If

    
  End Subroutine prepDataStore 
   
  
  Subroutine synchMpiProcesses()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i,send,receive,tag 
    Integer(kind=StandardInteger) :: processID,processCount,error,status 
!call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!Send out data to all processes
    send = 0
    If(processID.eq.0) Then
      Do i=1,(processCount-1)
        tag = 1000 + i	
        Call MPI_send(send,1,MPI_integer,i,tag,&
		MPI_comm_world,error)
      End Do
    Else
      tag = 1000 + processID
	  Call MPI_recv(receive,1,MPI_integer,0,tag,&
      MPI_comm_world,status,error)  
    End If  
  End Subroutine synchMpiProcesses 
  
  
  
  
  
  
  
  
  
  
  
  
  
  

End Module prep