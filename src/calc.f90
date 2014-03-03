Module calc

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use prep


!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'   
!declare global variables 
  Integer(kind=StandardInteger) :: functionOutOfRangeCounter
!Privacy of functions/subroutines/variables
  Private
!variables
  Public :: functionOutOfRangeCounter 
!subroutines
  Public :: runCalc		
  Public :: calcConfigEnergies
!functions
  

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
	!outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	!open(unit=999,file=trim(outputFile),status="old",position="append",action="write")		
!write to output file
    !write(999,"(A36,F8.4)") "Start calculation:                  ",ProgramTime()	
!initialise any arrays etc	
	!Call applyUnitVector(0)
	!Call calcConfigEnergies()	
    !if(calcRunType(1:6).eq."ENERGY")then
!write to output file
      !write(999,"(A36,F8.4)") "Calculating configuration energies: ",ProgramTime()
	  !Call applyUnitVector(0)
	  !Call calcConfigEnergies()
	!endif	
	!if(calcRunType(1:11).eq."BULKMODULUS")then
!write to output file
      !write(999,"(A26,F8.4)") "Calculating bulk modulus: ",ProgramTime()
	  !Call calcBulkModulus()
	!endif
!write to output file
    !write(999,"(A36,F8.4)") "End calculation:                  ",ProgramTime()  
    !close(999)	
  End Subroutine runCalc
  

  
  
!------------------------------------------------------------------------!
! ***calcConfigEnergies*** 
!------------------------------------------------------------------------!

  Subroutine calcConfigEnergies()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Real(kind=SingleReal) ::  startTime, endTime	
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray
!Start Time
    startTime = ProgramTime()
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	write(999,"(A32,F8.4)") "Calculating configuration energy",ProgramTime()
!Allocate array
    !Allocate(configurationEnergy(1:configCount))
!zero out energies
    do i=1,configCount
      configurationEnergy(i) = 0.0D0	
	enddo
!Loop through configurations, and split over MPI processes
	Do i=1,configCount
	  selectProcess = mod(i-1,mpiProcessCount)
	  If(selectProcess.eq.mpiProcessID)Then
	    configurationEnergy(i) = CalcEnergy(i)
	  End If
	End Do
!Collect data from processes
    Allocate(bufferArray(1:configCount))
    If(mpiProcessID.gt.0) Then
!send buffers from all worker processes
      processTo = 0
      Do i=1,(mpiProcessCount-1)
        If(i.eq.mpiProcessID)Then
          tag = 1000 + mpiProcessID	
          Call MPI_send(configurationEnergy,configCount,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End If
      End Do
    End If
    If(mpiProcessID.eq.0) Then
!read buffers from worker processes
      Do j=1,(mpiProcessCount-1)
		processFrom = j
        tag = 1000 + j
        Call MPI_recv(bufferArray,configCount,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
        Do i=1,ceiling(configCount/(1.0D0*mpiProcessCount))
	      k = (i-1)*mpiProcessCount+j+1
	      If(k.le.configCount)Then
	        configurationEnergy(k) = bufferArray(k)
	      End If
	    End Do  
      End Do
    End If
!send array from master to worker processes
	Call MPI_sendout(configurationEnergy,"double")

!Write to file, if master process
	If(mpiProcessID.eq.0)Then
	  write(999,"(A44,I4)") "--------------------------------------------"
	  write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	  write(999,"(A22,I4)") "Configuration Details "
	  write(999,"(A44,I4)") "--------------------------------------------"
	  Do i=1,configCount
	    write(999,"(I4,A1,F12.6,A4,I8,A6)") i," ",&
	      (1.0*configurationEnergy(i))," eV ",&
	    configAtoms(i)," atoms "
	    write(999,"(I4,F12.6,A12)") i,(configurationEnergy(i)/configAtoms(i))," eV per atom"
	    write(999,"(A1)") " "
	  End Do   
	End If  
!End Time
    endTime = ProgramTime() 
	write(999,"(A18,F16.8)") "Energy calc time: ",(endTime-startTime)
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
	  
	  bmFit = PolyFit(dataPoints,polyOrder)
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
	Integer(kind=StandardInteger) :: densityCount
	Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, energy, time
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: density  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayV  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayP  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayE  
!Info	
!neighbourListI(j,1) A type, neighbourListI(j,2) B type, 
!neighbourListI(j,3) A id, neighbourListI(j,4) B id
!Set Variables
    configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)
    atoms = configAtoms(configurationID)
	densityCount = atoms*elementCount
	functionOutOfRangeCounter = 0
!Allocate density array
    if(Allocated(density))then
	  Deallocate(density)
	endif
    Allocate(density(1:atoms))
!zero density
    do j=1,atoms
	  density(j) = 0.0D0
	enddo
!pair energy contribution and density calculation
	pairEnergy = 0.0D0
	embeddingEnergy = 0.0D0
!sum pair energy and density	
	do j=configStart,configStart+configLength-1 
	  yArrayV = SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),1,neighbourListR(j))
	  pairEnergy = pairEnergy + yArrayV(1)
	  !k = neighbourListI(j,3) + (atoms * (neighbourListI(j,2) - 1))	
	  !yArrayP = SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),2,neighbourListR(j))
	  !density(k) = density(k) + yArrayP(1)
	  yArrayP = SearchPotentialPoint(neighbourListI(j,2),0,2,neighbourListR(j))
	  density(neighbourListI(j,3)) = density(neighbourListI(j,3)) + yArrayP(1)
	enddo	
	pairEnergy = 0.5D0 * pairEnergy
!sum embedding energy
    do j=1,atoms
	  !do l=1,elementCount
	    !k = atoms + (atoms*(l-1))
		!yArrayE = SearchPotentialPoint(neighbourListI(j,1),neighbourListI(j,2),3,density(k))
	    !embeddingEnergy = embeddingEnergy + yArrayE(1)
	  !enddo
	  yArrayE = SearchPotentialPoint(neighbourListI(j,1),0,3,density(j))
	  embeddingEnergy = embeddingEnergy + yArrayE(1)
	enddo
!total energy of configuration
    energy = pairEnergy + embeddingEnergy   
  End function CalcEnergy  

  
  
  
  
 Function CalcForces (configurationID) RESULT (energy)  
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, configurationID  
	Integer(kind=StandardInteger) :: configStart, configLength
	Integer(kind=StandardInteger) :: atoms
	Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, energy, time
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: density  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayV  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayP  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayE  
!Set Variables
    configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)
    atoms = configAtoms(configurationID)
!total energy of configuration







  End function CalcForces    
  
  
  
  Function SearchPotentialPoint (atomA, atomB, potType, x, calcForce) RESULT (yArray)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: output
    Integer(kind=StandardInteger) :: outOfRange, inRange
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger) :: atomA, atomB, potType, potKey
	Integer(kind=StandardInteger) :: atomMin, atomMax
	Integer(kind=StandardInteger) :: potStart, potLength, dataPos, dataPosNext
	Integer(kind=StandardInteger) :: dataPosMin, dataPosMax, dataPosOffset
	Integer(kind=StandardInteger) :: interpPointCount
	Real(kind=DoubleReal) :: x, y, xLower, xUpper, yLower, yUpper, xPos, yPos
	Real(kind=DoubleReal) :: xPosUpper, xPosLower
	Real(kind=DoubleReal) :: xDouble, yDouble
	Real(kind=DoubleReal) :: xA,xB,xC,xD,yA,yB,yC,yD,aA,aB,aC,A,B,C
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: interpPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray  	
!optional function inputs
	Integer(kind=StandardInteger),optional :: calcForce
	Integer(kind=StandardInteger) :: calcForceOption
!set optional variables	
	If(present(calcForce))Then
	  calcForceOption = calcForce
	Else
	  calcForceOption = 0
	End If
	!potType = 1	Pair
	!potType = 2	Dens or S-Band Dens
	!potType = 3	Embe or S-Band Embe
	!potType = 4	D-Band Dens (dend)
	!potType = 5	D-Band Embe (embd)	
	!eamType = 1	Standard EAM
	!eamType = 2	2BMEAM	
	!pairCount, densCount, dendCount, embeCount, embdCount, elementCount
!Set result array 	
	Allocate(yArray(1:2))
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
!upper and lower values
    xLower = eamData(potStart,1)
    xUpper = eamData(potStart+potLength-1,1)
!if variable is not in range, zero
    If(x.lt.xLower.or.x.gt.xUpper)Then
	  functionOutOfRangeCounter = functionOutOfRangeCounter + 1
	  yArray(1) = 0.0D0	!Energy to 0 if out of range
	  yArray(2) = 0.0D0	!Energy to 0 if out of range
	Else
!estimate location
	  dataPos = (potStart-1)+floor(potLength*(x/(1.0D0*(xUpper-xLower))))
	  dataPosMin = potStart
	  dataPosMax = potStart+potLength-1
!check if estimate is in range
      outOfRange = 0
	  inRange = 0	  
	  If(dataPos.eq.dataPosMax)Then
	    xPosLower = 1.0D0*eamData(dataPos-1,1)
	    xPosUpper = 1.0D0*eamData(dataPos,1)
	  Else
	    xPosLower = 1.0D0*eamData(dataPos,1)
	    xPosUpper = 1.0D0*eamData(dataPos+1,1)
	  End If
	  inRange = 1 !temporary override
      !If(x.ge.xPosLower.and.x.le.xPosUpper)Then
      !  inRange = 1
	  !Else
	  !  Do k=1,1
	  !	  outOfRange = outOfRange + 1
!slot before
      !    dataPosNext = dataPos - outOfRange
	!	  If(dataPosNext.eq.dataPosMax)Then
	!        xPosLower = 1.0D0*eamData(dataPos-1,1)
	!        xPosUpper = 1.0D0*eamData(dataPos,1)
	!      Else
	!        xPosLower = 1.0D0*eamData(dataPos,1)
	!        xPosUpper = 1.0D0*eamData(dataPos+1,1)
	!      End If
	!	  If(x.ge.xPosLower.and.x.le.xPosUpper)Then
    !        inRange = 1
	!		dataPos = dataPosNext
	!		Exit
	 !     End If	
!slot after
     !     dataPosNext = dataPos + outOfRange
	!	  If(dataPosNext.eq.dataPosMax)Then
	!        xPosLower = 1.0D0*eamData(dataPos-1,1)
	!        xPosUpper = 1.0D0*eamData(dataPos,1)
	!      Else
	!        xPosLower = 1.0D0*eamData(dataPos,1)
	!        xPosUpper = 1.0D0*eamData(dataPos+1,1)
	!      End If
!		  If(x.ge.xPosLower.and.x.le.xPosUpper)Then
     !       inRange = 1
	!		dataPos = dataPosNext
	!		Exit
	!      End If	  
	!    End Do
	!  End If	  
      If(inRange.eq.0)Then !If out of range
	    functionOutOfRangeCounter = functionOutOfRangeCounter + 1
	    yArray(1) = 0.0D0	!Energy to 0 if out of range
	    yArray(2) = 0.0D0	!Energy to 0 if out of range
      Else
!in range, set 
	    xPos = eamData(dataPos,1)
	    yPos = eamData(dataPos,2)	
!Point interpolation using lagrange's formula
        If(eamInterpType.lt.2.or.eamInterpType.gt.7)Then
		  eamInterpType = 4		!Set a default for 5 point interpolation 
		End If
        interpPointCount = eamInterpType + 1
        Allocate(interpPoints(1:interpPointCount,1:2))
		!dataPosMin dataPosMax
!prepare set of data points used in interp
		dataPosOffset = -1*floor(interpPointCount/2.0D0)
		If((dataPos+dataPosOffset).lt.dataPosMin)Then
		  dataPosOffset = 0  
		End If
		If((dataPos-dataPosOffset).gt.dataPosMax)Then
		  dataPosOffset = -1*interpPointCount  
		End If
		Do i=1,interpPointCount
		  interpPoints(i,1) = eamData(dataPos+dataPosOffset+i,1)
		  interpPoints(i,2) = eamData(dataPos+dataPosOffset+i,2)
		End Do
!calculate y and y' from interp points
	    yArray = PointInterpolationFull(interpPoints,x)
	  End If !in data slots 
	End If !in range of whole potential
  End function SearchPotentialPoint    
  

End Module calc