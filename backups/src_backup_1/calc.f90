Module calc

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use units
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
  Public :: calcEval		
  Public :: calcConfigEnergies
  Public :: calcDifference
!functions
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!------------------------------------------------------------------------!
! ***calcEval*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!

  Subroutine calcEval()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
!Calculate energies, forces and stresses
    Call calcConfigEnergies(.true.,.false.)
!Bulk modulus
    Call calcDifference()


  
    print *,configurationEnergy(1),trialResidualSquareSum
  End Subroutine calcEval
  
!------------------------------------------------------------------------!
! ***calcConfigEnergies*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!

  Subroutine calcConfigEnergies(printOut, saveForces)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k, configCounter, atomCounter
	Real(kind=SingleReal) ::  startTime, endTime	
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray
	Logical :: printOut, saveForces
!Start Time
    startTime = ProgramTime()
    If(mpiProcessID.gt.0) Then
!open output file	
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  write(999,"(A32,F8.4)") "Calculating configuration energy",ProgramTime()
	End If
!zero out energies
    Do i=1,configCount
      configurationEnergy(i) = 0.0D0	
	End Do
!zero out forces
    If(calcForcesOnOff.eq.1)Then
      Do i=1,configAtomsTotal
        configurationForceX(i) = 0.0D0
        configurationForceY(i) = 0.0D0
        configurationForceZ(i) = 0.0D0	
	  End Do
	End If  
!-----------------------------------------------------
! Split over MPI processes and calc energies, forces
!-----------------------------------------------------
!Loop through configurations, and split over MPI processes
	Do i=1,configCount
	  selectProcess = mod(i-1,mpiProcessCount)
	  configAtomsMap(i,3) = selectProcess
	  If(selectProcess.eq.mpiProcessID)Then
	    configurationEnergy(i) = CalcEnergy(i)
	  End If
	End Do
!-----------------------------------------------------------
! Collect energy results from all processes then distribute
!-----------------------------------------------------------	
!Collect data from processes
    If(Allocated(bufferArray))Then
	  Deallocate(bufferArray)
	End If
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
!-----------------------------------------------------------
! Collect force results from all processes then distribute
!-----------------------------------------------------------	
	!Collect data from processes
    If(Allocated(bufferArray))Then
	  Deallocate(bufferArray)
	End If
    Allocate(bufferArray(1:configAtomsTotal))
    If(mpiProcessID.gt.0) Then
!send buffers from all worker processes
      Do i=1,(mpiProcessCount-1)
        If(i.eq.mpiProcessID)Then
          processTo = 0
          tag = 1100 + mpiProcessID	
          Call MPI_send(configurationForceX,configAtomsTotal,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
          processTo = 0
          tag = 1200 + mpiProcessID	
          Call MPI_send(configurationForceY,configAtomsTotal,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
          processTo = 0
          tag = 1300 + mpiProcessID	
          Call MPI_send(configurationForceZ,configAtomsTotal,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End If
      End Do
    End If
    If(mpiProcessID.eq.0) Then
!read buffers from worker processes
      Do j=1,(mpiProcessCount-1)
!receive arrays - X force component
		processFrom = j
        tag = 1100 + j
        Call MPI_recv(bufferArray,configAtomsTotal,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
!loop and transfer data
        Do i=1,configCount
		  If(configAtomsMap(i,3).eq.j)Then !Check if this process j calculated forces fr config i
!Transfer data if calculated for those atoms
            Do k=configAtomsMap(i,1),configAtomsMap(i,2)
			  configurationForceX(k) = bufferArray(k)
			End Do
		  End If
		End Do
!receive arrays - Y force component
		processFrom = j
        tag = 1200 + j
        Call MPI_recv(bufferArray,configAtomsTotal,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
!loop and transfer data
        Do i=1,configCount
		  If(configAtomsMap(i,3).eq.j)Then !Check if this process j calculated forces fr config i
!Transfer data if calculated for those atoms
            Do k=configAtomsMap(i,1),configAtomsMap(i,2)
			  configurationForceY(k) = bufferArray(k)
			End Do
		  End If
		End Do
!receive arrays - Z force component
		processFrom = j
        tag = 1300 + j
        Call MPI_recv(bufferArray,configAtomsTotal,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
!loop and transfer data
        Do i=1,configCount
		  If(configAtomsMap(i,3).eq.j)Then !Check if this process j calculated forces fr config i
!Transfer data if calculated for those atoms
            Do k=configAtomsMap(i,1),configAtomsMap(i,2)
			  configurationForceZ(k) = bufferArray(k)
			End Do
		  End If
		End Do	
      End Do
    End If	
!send arrays from master to worker processes
	Call MPI_sendout(configurationForceX,"double")
	Call MPI_sendout(configurationForceY,"double")
	Call MPI_sendout(configurationForceZ,"double")
!-----------------------------------------------------------
! Collect Stress results from all processes then distribute
!-----------------------------------------------------------	
    If(mpiProcessID.eq.0)Then
    Do i=1,configCount
	  Do j=1,9
        !print "(I8,I8,E20.10)",i,j,configurationStress((i-1)*9+j)
      End Do
	End Do
    End If
	

!Write energies to output file, if master process
	If(mpiProcessID.eq.0.and.printOut.eqv..true.)Then
	  write(999,"(A44,I4)") "--------------------------------------------"
	  write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	  write(999,"(A22,I4)") "Configuration Details "
	  write(999,"(A44,I4)") "--------------------------------------------"
	  Do i=1,configCount
	    write(999,"(I4,A1,F12.6,A4,I8,A10,I8)") i," ",&
	      (1.0*configurationEnergy(i))," eV ",&
	    configAtoms(i)," atoms     ",configAtomsMap(i,3)
	    write(999,"(I4,F12.6,A14,F12.6,A8)") i,&
		(configurationEnergy(i)/configAtoms(i))," eV per atom (",&
		configurationRefEnergy(i)," eV Ref)"
		If(calcStressOnOff.eq.1)Then
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Sxx:",configurationStress((i-1)*9+1)," Sxy:",&
		  configurationStress((i-1)*9+2)," Sxz:",configurationStress((i-1)*9+3)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Syx:",configurationStress((i-1)*9+4)," Syy:",&
		  configurationStress((i-1)*9+5)," Syz:",configurationStress((i-1)*9+6)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Szx:",configurationStress((i-1)*9+7)," Szy:",&
		  configurationStress((i-1)*9+8)," Szz:",configurationStress((i-1)*9+9)
		End If
	    write(999,"(A1)") " "
	  End Do   
	End If  
!End Time
    endTime = ProgramTime() 
	If(mpiProcessID.eq.0.and.printOut.eqv..true.)Then
	  write(999,"(A18,F16.8)") "Energy calc time: ",(endTime-startTime)
	  close(999)	
	End If	
!Write forces to output file, if master process
	If(mpiProcessID.eq.0.and.saveForces.eqv..true.)Then	
!open forces output file	
	  open(unit=989,file=trim(outputFileForces),status="old",position="append",action="write")	  
	  configCounter = 1
	  atomCounter = 0
	  Do i=1,configAtomsTotal  
	    atomCounter = atomCounter + 1
	    If(configAtomsMap(configCounter,4).gt.0)Then	!Reference forces
	      write(989,"(I4,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,A10,I8,I8,I8)") atomCounter,&
		  configurationForceX(i),configurationForceY(i),configurationForceZ(i),&
		  configurationRefForceX(i),configurationRefForceY(i),&
		  configurationRefForceZ(i),&
		  "          ",i,configCounter,configAtomsMap(configCounter,3)
		Else
	      write(989,"(I4,E20.10,E20.10,E20.10,A10,I8,I8,I8)") atomCounter,&
		  configurationForceX(i),configurationForceY(i),configurationForceZ(i),&
		  "          ",i,configCounter,configAtomsMap(configCounter,3)
		End If 
		If(atomCounter.eq.configAtoms(configCounter))Then
		  atomCounter = 0
		  configCounter = configCounter + 1
		End If 
	  End Do	  
	  close(989)	
	End If
	
	
  End Subroutine calcConfigEnergies
  

!------------------------------------------------------------------------!
! ***calcConfigBulkProperties*** 
! Calculates the bulk properties of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!

  Subroutine calcConfigBulkProperties(printOut)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k, configCounter, atomCounter
	Real(kind=SingleReal) ::  startTime, endTime	
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray
	Logical :: printOut, saveForces
!Start Time
    startTime = ProgramTime()
    If(mpiProcessID.gt.0) Then
!open output file	
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  write(999,"(A44,F8.4)") "Calculating configuration bulk properties ",ProgramTime()
	End If
  
  
  End Subroutine calcConfigBulkProperties  
    
  
!------------------------------------------------------------------------!
! ***calcBulkModulus*** 
!
!------------------------------------------------------------------------!

Subroutine calcBulkModulusOld()
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
  End Subroutine calcBulkModulusOld
  
  
  
  
  
  
!------------------------------------------------------------------------!
! ***calcDifference*** 
!------------------------------------------------------------------------!
    
  Subroutine calcDifference(printOutIn)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: configCounter, atomCounter
	Real(kind=DoubleReal) :: energyDifference, forceDifference
	Logical, optional :: printOutIn
	Logical :: printOut
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag
!Set optional variables
	If(present(printOutIn))Then
      printOut = printOutIn
	Else
      printOut = .false.	
    End If
!Total residual square sum
    trialResidualSquareSum = 0.0D0  
!Energy Difference
    energyDifference = 0.0D0
    Do i=1,configCount
	  energyDifference = energyDifference + 1.0D0*(configurationEnergy(i)/configAtoms(i)-&
	  configurationRefEnergy(i))**2
	End Do  
	trialResidualSquareSum = energyDifference
!Force Difference
    forceDifference = 0.0D0
    If(calcForcesOnOff.eq.1)Then
	  configCounter = 1
	  atomCounter = 0
  	  Do i=1,configAtomsTotal 
	    atomCounter = atomCounter + 1
		If(configAtomsMap(configCounter,4).gt.0)Then  
          forceDifference=forceDifference+&
		    1.0D0*(configurationForceX(i)-configurationRefForceX(i))**2+&
			   1.0D0*(configurationForceY(i)-configurationRefForceY(i))**2+&
			      1.0D0*(configurationForceZ(i)-configurationRefForceZ(i))**2
		End If		  
		If(atomCounter.eq.configAtoms(configCounter))Then
		  atomCounter = 0
		  configCounter = configCounter + 1
		End If 		  
	  End Do
    End If
	
	!Write forces to output file, if master process
	If(mpiProcessID.eq.0.and.printOut.eqv..true.)Then	
!open forces output file	
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	  
	  write(999,"(A25,E20.10)") "Energy difference:       ",energyDifference	  
	  write(999,"(A25,E20.10)") "Force difference:        ",forceDifference
	  close(999)	
	End If
  
  End Subroutine calcDifference
  
  
  

  
  
  
  
  
  


  
  
  
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
	Integer(kind=StandardInteger) :: i, j, k, l, n, configurationID, keyIJ, forceKey  
	Integer(kind=StandardInteger) :: ii, jj
	Integer(kind=StandardInteger) :: configStart, configLength
	Integer(kind=StandardInteger) :: atoms, atomI, atomJ, atomIType, atomJType
	Integer(kind=StandardInteger) :: densityCount
	Real(kind=DoubleReal) :: pairEnergy, embeddingEnergy, energy, time
	Real(kind=DoubleReal) :: embedForce
	Real(kind=DoubleReal) :: pairEnergyDifferential
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: density  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayV  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayP  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayE 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayEi 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArrayEj 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: atomForce
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: atomForceMagnitude
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: atomForceComponent  
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: densityGrad  
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: embeddingGrad  
!Info	
!neighbourListI(j,1) A type, neighbourListI(j,2) B type, 
!neighbourListI(j,3) A id, neighbourListI(j,4) B id
!Set Variables
    configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)
    atoms = configAtoms(configurationID)
	densityCount = atoms*elementCount
	functionOutOfRangeCounter = 0
!Allocate arrays
    if(Allocated(density))Then
	  Deallocate(density)
	endif
    Allocate(density(1:atoms))
	If(calcForcesOnOff.eq.1)Then 
	  If(Allocated(atomForce))Then
	    Deallocate(atomForce)
	  End If
      Allocate(atomForce(1:atoms))
	  If(Allocated(atomForceMagnitude))Then
	    Deallocate(atomForceMagnitude)
	  End If
      Allocate(atomForceMagnitude(1:configLength))
	  If(Allocated(densityGrad))Then
	    Deallocate(densityGrad)
	  End If
      Allocate(densityGrad(1:configLength,1:2))
	  If(Allocated(embeddingGrad))Then
	    Deallocate(embeddingGrad)
	  End If
      Allocate(embeddingGrad(1:atoms))
	  If(Allocated(atomForceComponent))Then
	    Deallocate(atomForceComponent)
	  End If
	  Allocate(atomForceComponent(1:configLength,1:3))
	End If
!configurationForceMPIKey(i,1)
!zero density and atomForce
    Do j=1,atoms
	  density(j) = 0.0D0 	  	
	End Do
!pair energy contribution and density calculation
	pairEnergy = 0.0D0
	pairEnergyDifferential = 0.0D0
	embeddingEnergy = 0.0D0
!sum pair energy and density	
    k = 0
	do j=configStart,configStart+configLength-1 
	  k = k + 1
	  atomIType = neighbourListI(j,1)
	  atomJType = neighbourListI(j,2)
	  atomI = neighbourListI(j,3)
	  atomJ = neighbourListI(j,4)
!Pair potential
	  yArrayV = SearchPotentialPoint(atomIType,atomJType,1,neighbourListR(j))
	  pairEnergy = pairEnergy + yArrayV(1)
!Density function
	  yArrayP = SearchPotentialPoint(atomJType,0,2,neighbourListR(j))
	  density(atomI) = density(atomI) + yArrayP(1)
	  If(calcForcesOnOff.eq.1)Then	  
!Add pair potential to each force component
        atomForceMagnitude(k) = 1.0D0*yArrayV(2)
!Store density gradient for atom i-j pair
		densityGrad(k,2) = yArrayP(2)
		yArrayP = SearchPotentialPoint(atomIType,0,2,neighbourListR(j))
		densityGrad(k,1) = yArrayP(2)
	  End If
	End Do	
!sum embedding energy and embedding grad for each atom
    do j=1,atoms
	  n = configAtomsStart(configurationID) + j
	 !atomTypeKey
	  yArrayE = SearchPotentialPoint(atomTypeKey(n),0,3,density(j))
	  embeddingEnergy = embeddingEnergy + yArrayE(1)	
	enddo
!sum force	
    If(calcForcesOnOff.eq.1)Then	
      k = 0
	  Do j=configStart,configStart+configLength-1 
	    k = k + 1
		atomIType = neighbourListI(j,1)
		atomJType = neighbourListI(j,2)
		atomI = neighbourListI(j,3)
		atomJ = neighbourListI(j,4)
		yArrayEi = SearchPotentialPoint(atomIType,0,3,density(atomI))
		yArrayEj = SearchPotentialPoint(atomJType,0,3,density(atomJ))
		embedForce = 1.0D0*yArrayEi(2)*densityGrad(k,2)+1.0D0*yArrayEj(2)*densityGrad(k,1)
!Add embedding potential to each force component
        atomForceMagnitude(k) = atomForceMagnitude(k) + 1.0D0*embedForce
!sum up forces by atom
        forceKey = configAtomsStart(configurationID)+atomI
		configurationForceX(forceKey) = configurationForceX(forceKey) + &
		atomForceMagnitude(k) * neighbourListCoords(j,10)
		configurationForceY(forceKey) = configurationForceY(forceKey) + &
		atomForceMagnitude(k) * neighbourListCoords(j,11)
		configurationForceZ(forceKey) = configurationForceZ(forceKey) + &
		atomForceMagnitude(k) * neighbourListCoords(j,12)
	  End Do
	End If	
!calculate stress	
    If(calcStressOnOff.eq.1)Then
	  k = 0
	  Do j=configStart,configStart+configLength-1 
	    k = k + 1
		Do ii=1,3
	      Do jj=1,3
	        n = ((configurationID-1)*9)+3*(ii-1)+jj
			configurationStress(n) = configurationStress(n)+&
			neighbourListCoords(j,6+ii)*atomForceComponent(k,jj)
		  End Do  
		End Do 
	  End Do	
	End If
!total energy of configuration
    energy = 0.5D0 * pairEnergy + embeddingEnergy   
	
	
	
	
!Convert stress to GPa
	If(calcStressOnOff.eq.1)Then
      Do ii=1,3
	    Do jj=1,3
          n = ((configurationID-1)*9)+3*(ii-1)+jj
          configurationStress(n) = UnitConvert(configurationStress(n),"EVAN3","GPA")
          configurationStress(n) = &
	      (1.0D0/(2.0D0*configurationVolume(configurationID))) * configurationStress(n)
        End Do
      End Do		
	End If
  End function CalcEnergy  

  
  Function CalcBulkModulus (configurationID) RESULT (bulkProperties) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, n, configurationID, keyIJ, forceKey  
	Real(kind=DoubleReal) :: bulkProperties
    
  End function CalcBulkModulus  
  
  
  
  
  
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
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: interpPointsTemp
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
	eamDataSet = 2
	If(eamDataSet.eq.1)Then
!get the start point and length of potential data
      potStart = eamKey(potKey,4)
      potLength = eamKey(potKey,5)
      interpPointCount = eamInterpType + 1
!interpolate
	  yArray = PointInterpolationArr(eamData,x,interpPointCount,potStart,potLength)
	End If
	If(eamDataSet.eq.2)Then
!get the start point and length of potential data
      potStart = eamKeyReduced(potKey,4)
      potLength = eamKeyReduced(potKey,5)
      interpPointCount = eamInterpType + 1
!interpolate
	  yArray = PointInterpolationArr(eamDataReduced,x,interpPointCount,potStart,potLength)
	End If
  End function SearchPotentialPoint    
  

End Module calc