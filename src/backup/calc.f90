Module calc

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use neighbour


!force declaration of all variables
  Implicit None
  
!declare global variables 
  Real(kind=SingleReal), Dimension( : ), Allocatable :: configurationEnergy

!Privacy of functions/subroutines/variables
  Private
  Public :: runCalc		
  Public :: configurationEnergy				    
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runCalc()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
	
!write to output file
    write(999,"(A36,F8.4)") "Start calculation:                  ",ProgramTime()
	
	
    if(calcRunType(1:6).eq."ENERGY")then
!write to output file
      write(999,"(A36,F8.4)") "Calculating configuration energies: ",ProgramTime()
	  Call calcConfigEnergies()
	endif
	
	if(calcRunType(1:6).eq."BULKMODULUS")then
!write to output file
      write(999,"(A26,F8.4)") "Calculating bulk modulus: ",ProgramTime()
	  Call calcConfigEnergies()
	endif
	
  
  
!write to output file
    write(999,"(A36,F8.4)") "End calculation:                  ",ProgramTime()
  
    close(999)
	
  End Subroutine runCalc
  
  
  
  Subroutine calcConfigEnergies()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
	
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	
!write to output file
	write(999,"(A32,F8.4)") "Calculating configuration energy",ProgramTime()	
	
!Allocate array
    Allocate(configurationEnergy(1:configCount))
	
!zero out energies
    do i=1,configCount
      configurationEnergy(i) = 0.0	
	enddo	
	
	do i=1,configCount
	  configurationEnergy(i) = CalcEnergy(i)
	enddo
	
	
	write(999,"(A44,I4)") "--------------------------------------------"
	write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	write(999,"(A22,I4)") "Configuration Details "
	write(999,"(A44,I4)") "--------------------------------------------"
	do i=1,configCount
	  write(999,"(I4,A1,F12.6,A4,I8,A6)") i," ",(1.0*configurationEnergy(i))," eV ",&
	  configAtoms(i)," atoms"
	  write(999,"(I4,F12.6,A12)") i,(configurationEnergy(i)/configAtoms(i))," eV per atom"
	  write(999,"(A1)") " "
	enddo
    
	close(999)
	
  End Subroutine calcConfigEnergies
  
  
  
  
  
  
  
  
  

  Subroutine calcEnergyOld()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l  
	Integer(kind=StandardInteger) :: configStart, configLength
	Integer(kind=StandardInteger) :: atoms, elements
	Real(kind=SingleReal) :: pairEnergy, embeddingEnergy
	Real(kind=SingleReal), Dimension( : ), Allocatable :: density
	
	
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
	
!write to output file
    write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	write(999,"(A32)") "Calculating configuration energy"	
	
	
!Allocate array
    Allocate(configurationEnergy(1:size(neighbourListKey)/2))	
	
!zero out energies
    do i=1,configCount
      configurationEnergy(i) = 0.0	
	enddo	
		  
    do i=1,configCount
!loop through each configuration separately
      configStart = neighbourListKey(i,1)
      configLength = neighbourListKey(i,2)
      atoms = configAtoms(i)
	  
!Allocate density array
      if(Allocated(density))then
	    Deallocate(density)
	  endif
      Allocate(density(1:(atoms*elementCount)))
!zero density
      do j=1,(atoms*elementCount)
	    density(j) = 0.0
	  enddo
	  
!pair energy contribution and density calculation
	  pairEnergy = 0.0
	  embeddingEnergy = 0.0
!sum pair energy and density	
	  do j=configStart,configStart+configLength-1 
		pairEnergy = pairEnergy + &
		SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),1,neighbourListR(j))
		k = neighbourListI(j,3) + (atoms * (neighbourListI(j,2) - 1))
		density(k) = density(k) + &
		SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),2,neighbourListR(j))		
	  enddo	
	  pairEnergy = 0.5 * pairEnergy
!sum embedding energy
      do j=1,atoms
	    do l=1,elementCount
	      k = atoms + (atoms*(l-1))
		  embeddingEnergy = embeddingEnergy + &
		  SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),3,density(k))
	    enddo
	  enddo
!total energy of configuration
      configurationEnergy(i) = pairEnergy + embeddingEnergy
	enddo
	
	
	write(999,"(A44,I4)") "--------------------------------------------"
	write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	write(999,"(A22,I4)") "Configuration Details "
	write(999,"(A44,I4)") "--------------------------------------------"
	do i=1,configCount
	  write(999,"(I4,A1,F12.6,A4,I8,A6)") i," ",(1.0*configurationEnergy(i))," eV ",&
	  configAtoms(i)," atoms"
	  write(999,"(I4,F12.6,A12)") i,(configurationEnergy(i)/configAtoms(i))," eV per atom"
	write(999,"(A1)") " "
      
	enddo
	
  	
!close output file
    close(999)	
  
  End Subroutine calcEnergyOld 
  
  
  
  
  
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!   
  
  Function CalcEnergy (configurationID) RESULT (energy)
  
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, configurationID  
	Integer(kind=StandardInteger) :: configStart, configLength
	Integer(kind=StandardInteger) :: atoms
	Real(kind=SingleReal) :: pairEnergy, embeddingEnergy, energy
	Real(kind=SingleReal), Dimension( : ), Allocatable :: density  
	
    configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)
    atoms = configAtoms(configurationID)
    
!Allocate density array
    if(Allocated(density))then
	  Deallocate(density)
	endif
    Allocate(density(1:(atoms*elementCount)))
!zero density
    do j=1,(atoms*elementCount)
	  density(j) = 0.0
	enddo
	  
!pair energy contribution and density calculation
	pairEnergy = 0.0
	embeddingEnergy = 0.0
!sum pair energy and density	
	do j=configStart,configStart+configLength-1 
	  pairEnergy = pairEnergy + &
	  SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),1,neighbourListR(j))
	  k = neighbourListI(j,3) + (atoms * (neighbourListI(j,2) - 1))
	  density(k) = density(k) + &
	  SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),2,neighbourListR(j))		
	enddo	
	pairEnergy = 0.5 * pairEnergy
!sum embedding energy
    do j=1,atoms
	  do l=1,elementCount
	    k = atoms + (atoms*(l-1))
	    embeddingEnergy = embeddingEnergy + &
	    SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),3,density(k))
	  enddo
	enddo
!total energy of configuration
    energy = pairEnergy + embeddingEnergy   
  
  End function CalcEnergy    
  
  
  Function SearchPotentialPoint (atomA, atomB, potType, x) RESULT (y)
    Integer(kind=StandardInteger) :: output 
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger) :: atomA, atomB, potType, potKey
	Integer(kind=StandardInteger) :: atomMin, atomMax
	Integer(kind=StandardInteger) :: potStart, potLength, dataPos
	Real(kind=SingleReal) :: x, y, xLower, xUpper, yLower, yUpper, xPos, yPos
	
	!potType = 1	Pair
	!potType = 2	Dens or S-Band Dens
	!potType = 3	Embe or S-Band Embe
	!potType = 4	D-Band Dens (dend)
	!potType = 5	D-Band Embe (embd)

	!eamType = 1	Standard EAM
	!eamType = 2	2BMEAM
	
	!pairCount, densCount, dendCount, embeCount, embdCount, elementCount
	
!Get the potential key	
	if(potType.eq.1)then
	  atomMax = max(atomA,atomB)-1
      atomMin = min(atomA,atomB)-1
	  potKey = 1+atomMin+(atomMax*(atomMax+1))/2
	endif
    if(potType.eq.2)then
	  if(eamType.eq.1)then
	    potKey = pairCount + atomA
	  endif
	  if(eamType.eq.2)then
	    if(elementCount.eq.1)then
	      potKey = pairCount + 1
		else
		  atomMax = max(atomA,atomB)-1
          atomMin = min(atomA,atomB)-1
          potKey = pairCount+1+atomMin+(atomMax*(atomMax-1))/2
		endif
	  endif
	endif
    if(potType.eq.3)then
	  if(eamType.eq.1)then
	    potKey = pairCount + densCount + atomA
	  endif
	  if(eamType.eq.2)then
	    potKey = pairCount + densCount + dendCount + atomA
	  endif
	endif
    if(potType.eq.4)then
	  if(eamType.eq.2)then
	    potKey = pairCount + densCount + atomA
	  endif
	endif
    if(potType.eq.5)then
	  if(eamType.eq.2)then
	    potKey = pairCount + densCount + dendCount + embeCount + atomA
	  endif
	endif
	
!get the start point and length of potential data
    potStart = eamKey(potKey,4)
    potLength = eamKey(potKey,5)
!estimate location
    xLower = eamData(potStart,1)
    xUpper = eamData(potStart+potLength-1,1)
	dataPos = (potStart-1)+floor(potLength*(x/(xUpper-xLower)))
	!print *,potType,potStart,potLength,x,xLower,xUpper,dataPos
	!if(potType.ne.3)then
	xPos = eamData(dataPos,1)
	yPos = eamData(dataPos,2)
	
!narrow down location
    if(x.eq.xPos)then
!done
      y = eamData(dataPos,2)
	else
!interpolate between data points
      if(x.gt.xPos)then
	    do j=1,100
	      if(x.ge.eamData(dataPos,1).and.x.le.eamData(dataPos+1,1))then
		    exit  !break out of loop
		  else
		    dataPos = dataPos + 1
		  endif
	    enddo
		xPos = eamData(dataPos,1)
		xLower = eamData(dataPos,1)
		xUpper = eamData(dataPos+1,1)
		yLower = eamData(dataPos,2)
		yUpper = eamData(dataPos+1,2)
		y = yLower + ((x-xLower)/(xUpper-xLower))*(yUpper-yLower)
	  endif
	  if(x.lt.xPos)then
	    do j=1,100
	      if(x.le.eamData(dataPos,1).and.x.ge.eamData(dataPos-1,1))then
		    exit  !break out of loop
		  else
		    dataPos = dataPos + 1
		  endif
	    enddo
		xPos = eamData(dataPos,1)
		xLower = eamData(dataPos-1,1)
		xUpper = eamData(dataPos,1)
		yLower = eamData(dataPos-1,2)
		yUpper = eamData(dataPos,2)
		y = yLower + ((x-xLower)/(xUpper-xLower))*(yUpper-yLower)
	  endif
	endif
	
	
	!endif
	
  End function SearchPotentialPoint    
  

End Module calc