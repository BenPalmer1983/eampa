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
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationEnergy
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationVolume
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationBM
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptVolume
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptEnergy
  Integer(kind=StandardInteger) :: functionOutOfRangeCounter

!Privacy of functions/subroutines/variables
  Private
  Public :: runCalc		
  Public :: configurationEnergy	
  Public :: configurationVolume	
  Public :: configurationBM	
  Public :: configurationOptVolume	
  Public :: configurationOptEnergy	
  Public :: functionOutOfRangeCounter  
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!------------------------------------------------------------------------!
! ***runCalc*** Run all the input subroutines
!------------------------------------------------------------------------!

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
	
!initialise any arrays etc
	
	
    if(calcRunType(1:6).eq."ENERGY")then
!write to output file
      write(999,"(A36,F8.4)") "Calculating configuration energies: ",ProgramTime()
	  Call applyUnitVector(0)
	  Call calcConfigEnergies()
	endif
	
	if(calcRunType(1:11).eq."BULKMODULUS")then
!write to output file
      write(999,"(A26,F8.4)") "Calculating bulk modulus: ",ProgramTime()
	  Call calcBulkModulus()
	endif
	
  
  
!write to output file
    write(999,"(A36,F8.4)") "End calculation:                  ",ProgramTime()  
    close(999)	
  End Subroutine runCalc
  
  
  
!------------------------------------------------------------------------!
! ***runCalcInitialise*** 
!------------------------------------------------------------------------!
  
  Subroutine runCalcInitialise()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")		
!write to output file
    write(999,"(A36,F8.4)") "Start runCalcInitialise:            ",ProgramTime()
!allocate configurationVolume and fill
    Allocate(configurationVolume(1:configCount))	
	do i=1,configCount
	  configurationVolume(i) = 0.0D0
	enddo
!Apply unit vector to neighbour list 
    write(999,"(A36,F8.4)") "Apply unit vector to neighbour list ",ProgramTime()
    Call applyUnitVector(0)
  
  
!write to output file
    write(999,"(A36,F8.4)") "End runCalcInitialise:              ",ProgramTime()  
    close(999)	 
  End Subroutine runCalcInitialise 
  
  
  
  
!------------------------------------------------------------------------!
! ***calcConfigEnergies*** 
!------------------------------------------------------------------------!

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
      configurationEnergy(i) = 0.0D0	
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
  
    
  
!------------------------------------------------------------------------!
! ***calcBulkModulus*** 
!------------------------------------------------------------------------!

Subroutine calcBulkModulus()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,polyOrder,bmPoints,bmStrain	
	Real(kind=SingleReal), Dimension( : , : ), Allocatable :: unitVectorOriginal
	Real(kind=SingleReal) :: systemEnergy, volume
	Real(kind=DoubleReal) :: minimumVolume, bulkModulus
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bmFit
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: volumeCoeffs
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPoints	
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A17,F8.4)") "Calc Bulk Modulus",ProgramTime()
!Allocate array
    Allocate(unitVectorOriginal(1:3,1:3))
    Allocate(configurationOptEnergy(1:configCount))
	Allocate(configurationBM(1:configCount))
	Allocate(configurationOptVolume(1:configCount))
!zero out energies
    do i=1,configCount
      !configurationEnergy(i) = 0.0	
	enddo
!store original diagonal components of unit vectors
    do i=1,3
	  do j=1,3
        unitVectorOriginal(i,j) = unitVector(i,j)
	  enddo
	enddo
!allocate BM fit array	and datapoints array
    polyOrder = 3
	bmPoints = 20 !(2*bmPoints+1)
	bmStrain = 8 !max strain percent
	Allocate(bmFit(0:polyOrder))
	Allocate(volumeCoeffs(0:(polyOrder-1)))
	Allocate(dataPoints(1:(2*bmPoints+1),1:2))
!loop through configurations and calculate bulk modulus
	do i=1,configCount	  
	  do j=(-1*bmPoints),bmPoints
!homogemeous strain 
        unitVector(1,1) = unitVectorOriginal(1,1) + &
		    unitVectorOriginal(1,1)*((0.1*j)/100)
        unitVector(2,2) = unitVectorOriginal(2,2) + &
		    unitVectorOriginal(2,2)*((0.1*j)/100)
        unitVector(3,3) = unitVectorOriginal(3,3) + &
		    unitVectorOriginal(3,3)*((0.1*j)/100)
		Call applyUnitVector(i)
		volume = configurationVolume(i)
		systemEnergy = CalcEnergy(i)
		dataPoints((j+bmPoints+1),1) = volume 
		dataPoints((j+bmPoints+1),2) = systemEnergy 
	    		
      enddo
	  do j=1,(2*bmPoints+1)
	    write(999,"(A26,I8,I8,F17.10,F17.10)") "Config/strain/vol/energy: ",&
		i,j,dataPoints(j,1),dataPoints(j,2)
	  enddo
	  
	  bmFit = polynomialFit(dataPoints,polyOrder)
!find minimum volume
      do j=1,polyOrder
	    volumeCoeffs(j-1) = j*bmFit(j)
	  enddo
      minimumVolume=SolvePolynomial(volumeCoeffs,1.0D0*&
	  dataPoints(1,1),1.0D0*dataPoints((2*bmPoints+1),1))	  
!calc bulk modulus
	  bulkModulus = minimumVolume * calcPolynomial (bmFit, minimumVolume, 2)
	  bulkModulus = bulkModulus * elementaryCharge * 1.0D30 * 1.0D-9
	  write(999,"(A7,I8)") "Config ",i
	  write(999,"(A32,F12.6)") "Minimum volume/Angstom:         ",minimumVolume
	  write(999,"(A32,F12.6)") "Bulk modulus/GPa:               ",bulkModulus	  
!store data in global arrays
	  configurationOptEnergy(i) = calcPolynomial (bmFit, minimumVolume, 0)
	  configurationBM(i) = bulkModulus
	  configurationOptVolume(i) = minimumVolume
	enddo
!reload original unit vector
	do i=1,3
	  do j=1,3
        unitVector(i,j) = unitVectorOriginal(i,j)
	  enddo
	enddo
!output data
	write(999,"(A44,I4)") "--------------------------------------------"
	write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	write(999,"(A22,I4)") "Configuration Details "
	write(999,"(A44,I4)") "--------------------------------------------"
	do i=1,configCount
	  !write(999,"(I4,A1,F12.6,A4,I8,A6)") i," ",(1.0*configurationEnergy(i))," eV ",&
	  !configAtoms(i)," atoms"
	  !write(999,"(I4,F12.6,A12)") i,(configurationEnergy(i)/configAtoms(i))," eV per atom"
	  !write(999,"(A1)") " "
	enddo    
	close(999)	
  End Subroutine calcBulkModulus
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
! ***applyUnitVector*** 
!------------------------------------------------------------------------!
  
  Subroutine applyUnitVector(configurationID)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k			
!open output file	
	!outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	!open(unit=9991,file=trim(outputFile),status="old",position="append",action="write")	
!loop through configurations
    if(configurationID.eq.0)then
	  do i=1,configCount
!write to output file
        write(999,"(A40,I8,F8.4)") "Apply unit vector to neighbour list rd: ",i,ProgramTime()
        Call applyUnitVectorAction(i)
	  enddo
	else
      write(999,"(A40,I8,F8.4)") "Apply unit vector to neighbour list rd: ",configurationID,ProgramTime()
	  Call applyUnitVectorAction(configurationID)
	endif
!close output file
    !close(9991)	 	
  End Subroutine applyUnitVector
  
  
!------------------------------------------------------------------------!
! ***applyUnitVectorAction*** 
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
!recalculate displacement between atoms
	do i=configStart,configStart+configLength-1 
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
    if(Allocated(configurationVolume))then
	else
	  Allocate(configurationVolume(1:configCount))
	endif
  
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
    configurationVolume(configurationID) = 1.0 * volume
    
    !unitVector(1,1)
  
  End Subroutine calcConfigurationVolume
    
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
	Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, energy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: density  
	
    configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)
    atoms = configAtoms(configurationID)
	functionOutOfRangeCounter = 0
    
!Allocate density array
    if(Allocated(density))then
	  Deallocate(density)
	endif
    Allocate(density(1:(atoms*elementCount)))
!zero density
    do j=1,(atoms*elementCount)
	  density(j) = 0.0D0
	enddo
	  
!pair energy contribution and density calculation
	pairEnergy = 0.0D0
	embeddingEnergy = 0.0D0
!sum pair energy and density	
	do j=configStart,configStart+configLength-1 
	  pairEnergy = pairEnergy + &
	  SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),1,neighbourListR(j))
	  k = neighbourListI(j,3) + (atoms * (neighbourListI(j,2) - 1))
	  density(k) = density(k) + &
	  SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),2,neighbourListR(j))		
	enddo	
	pairEnergy = 0.5D0 * pairEnergy
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
	Integer(kind=StandardInteger) :: dataPosMin, dataPosMax
	Real(kind=DoubleReal) :: x, y, xLower, xUpper, yLower, yUpper, xPos, yPos
	Real(kind=DoubleReal) :: xDouble, yDouble
	Real(kind=DoubleReal) :: xA,xB,xC,xD,yA,yB,yC,yD,aA,aB,aC,A,B,C
	
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
	dataPosMin = potStart
	dataPosMax = potStart+potLength-1
!check if out of range
    if(dataPos.lt.dataPosMin.or.(dataPos+1).gt.dataPosMax)then
!out of range, return zero
      functionOutOfRangeCounter = functionOutOfRangeCounter + 1
	  y = 0
	else
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
!nearest data point
          if(eamInterpType.eq.0)then
		    xLower = eamData(dataPos,1)
		    xUpper = eamData(dataPos+1,1)
			if((x-xLower).gt.(xUpper-x))then
			  y = eamData(dataPos+1,2) !y upper
			else
			  y = eamData(dataPos,2) !y lower
			endif
		  endif
!linear interpolation
          if(eamInterpType.eq.1)then
		    xPos = eamData(dataPos,1)
		    xLower = eamData(dataPos,1)
		    xUpper = eamData(dataPos+1,1)
		    yLower = eamData(dataPos,2)
		    yUpper = eamData(dataPos+1,2)
		    y = yLower + ((x-xLower)/(xUpper-xLower))*(yUpper-yLower)
		  endif
!threepoint interpolation
          if(eamInterpType.eq.2)then
!dataPos-1, dataPos, dataPos+1
!choose data points to interpolate between
		    if((dataPos-1).lt.dataPosMin)then
			  xA = eamData(dataPos,1)
			  xB = eamData(dataPos+1,1)
			  xC = eamData(dataPos+2,1)
			  yA = eamData(dataPos,2)
			  yB = eamData(dataPos+1,2)
			  yC = eamData(dataPos+2,2)
			elseif((dataPos+1).gt.dataPosMax)then
			  xA = eamData(dataPos-2,1)
			  xB = eamData(dataPos-1,1)
			  xC = eamData(dataPos,1)
			  yA = eamData(dataPos-2,2)
			  yB = eamData(dataPos-1,2)
			  yC = eamData(dataPos,2)
			else
			  xA = eamData(dataPos-1,1)
			  xB = eamData(dataPos,1)
			  xC = eamData(dataPos+1,1)
			  yA = eamData(dataPos-1,2)
			  yB = eamData(dataPos,2)
			  yC = eamData(dataPos+1,2)
		    endif			
!interpolate
			y = 1.0E0*QuadraticInterpolationCalc&
			   (1.0D0*xA,1.0D0*yA,1.0D0*xB,1.0D0*yB,1.0D0*xC,1.0D0*yC,1.0D0*x)
		  endif		  
!fourpoint interpolation
          if(eamInterpType.eq.3)then
!dataPos-1, dataPos, dataPos+1
!choose data points to interpolate between
            if((dataPos-2).lt.dataPosMin)then
			  dataPos = dataPos + 2
			endif
			if((dataPos-1).lt.dataPosMin)then
			  dataPos = dataPos + 1
			endif			
			if((dataPos+2).gt.dataPosMax)then
			  dataPos = dataPos - 2
			endif		
			if((dataPos+1).gt.dataPosMax)then
			  dataPos = dataPos - 1
			endif			
			xA = eamData(dataPos,1)
			xB = eamData(dataPos+1,1)
			xC = eamData(dataPos+2,1)
			xD = eamData(dataPos+3,1)
			yA = eamData(dataPos,2)
			yB = eamData(dataPos+1,2)
			yC = eamData(dataPos+2,2)	
			yD = eamData(dataPos+3,2)			
!interpolate
			y = 1.0E0*CubicInterpolationCalc&
			   (1.0D0*xA,1.0D0*yA,1.0D0*xB,1.0D0*yB,1.0D0*&
			   xC,1.0D0*yC,1.0D0*xD,1.0D0*yD,1.0D0*x)
		  endif
		  
		  
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
	endif
	
	
	!endif
	
  End function SearchPotentialPoint    
  

End Module calc