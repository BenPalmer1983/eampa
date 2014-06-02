Module calc

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use general
  Use strings		!string functions
  Use maths
  Use units
  Use initialise
  Use input
  Use output
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
  Public :: calcEvalFull		
  Public :: calcConfigEnergies
  Public :: calcDifference
  !Public :: calcOutput
  Public :: calcOutputForces
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
!Clear values
    trialResidualSquareSum = 0.0D0	
!Count energy calculation
    globalCounter(2) = globalCounter(2) + 1
!Calculate energies, forces and stresses
    Call calcConfigEnergies()
!Bulk modulus
    !Call calcBulkModulii()

!Synch MPI processes
	Call synchMpiProcesses()
!Calc difference 
    Call calcDifference()
!Save difference to file
	Call calcOutputDifference() 
  End Subroutine calcEval
  
  
  
  Subroutine calcEvalFull()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
!Clear values
    trialResidualSquareSum = 0.0D0	
!Count energy calculation
    globalCounter(2) = globalCounter(2) + 1
!Calculate energies, forces and stresses
    Call calcConfigEnergies()
!Equilibrium Volume	
	Call calcEquilibriumVolume()
!Bulk modulus
    Call calcBulkModulii()
!Bulk modulus
    !Call calcCubicElasticConstants()
!Synch MPI processes
	Call synchMpiProcesses()
!Calc difference 
    Call calcDifference()
!Save difference to file
	Call calcOutputDifference() 
  End Subroutine calcEvalFull
  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
! ***calcConfigEnergies*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!

  Subroutine calcConfigEnergies()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n,configCounter,atomCounter
	Real(kind=SingleReal) :: startTime, endTime	
	Logical :: printOut, saveForces
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,bufferArraySize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray
!Count energy calculation
    globalCounter(1) = globalCounter(1) + 1
!Start Time
    startTime = ProgramTime()
    If(mpiProcessID.eq.0.and.printOut.eqv..true.)Then
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  write(999,"(A33,F8.4,A3,I4)") "Calculating configuration energy ",&
	  ProgramTime(),"   ",globalCounter(1)
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
    Call synchMpiProcesses() !Synchronise mpi processes at this point
!Loop through configurations, and split over MPI processes
	Do i=1,configCount
	  selectProcess = mod(i-1,mpiProcessCount)
	  configAtomsMap(i,3) = selectProcess
	  If(selectProcess.eq.mpiProcessID)Then
	    configurationEnergy(i) = CalcEnergy(i)
	  End If
	End Do
	Call synchMpiProcesses() !Synchronise mpi processes at this point
!-----------------------------------------------------------
! Collect energy results from all processes then distribute
!-----------------------------------------------------------
	Call distributeData1DDP(configurationEnergy,size(configurationEnergy,1),21)
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
	Call sendData1DDP(configurationForceX,size(configurationForceX))
	Call sendData1DDP(configurationForceY,size(configurationForceY))
	Call sendData1DDP(configurationForceZ,size(configurationForceZ))
!-----------------------------------------------------------
! Collect Stress results from all processes then distribute
!-----------------------------------------------------------	

    If(mpiProcessID.gt.0) Then
!send buffers from all worker processes
      Do i=1,(mpiProcessCount-1)
        If(i.eq.mpiProcessID)Then
          processTo = 0
          tag = 1400 + mpiProcessID	
		  bufferArraySize = (9*configCount)
          Call MPI_send(configurationStress,bufferArraySize,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End If
      End Do
    End If
!read buffers from worker processes
	If(mpiProcessID.eq.0) Then
!read buffers from worker processes
      Do j=1,(mpiProcessCount-1)
!receive arrays - X force component
		processFrom = j
        tag = 1400 + j
		bufferArraySize = (9*configCount)
        Call MPI_recv(bufferArray,bufferArraySize,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
		n = 0
		Do k=1,configCount
		  Do i=1,9
		    n = n + 1
			If(configAtomsMap(k,21).eq.j)Then
			  configurationStress(n) = bufferArray(n) 
			End If
		  End Do
		End Do
      End Do
	End If  
!send array from master to worker processes
	Call sendData1DDP(configurationStress,size(configurationStress))
!Set end time	
    endTime = ProgramTime() 	
!-----------------------------
! Deallocate arrays
!-----------------------------
	If(Allocated(bufferArray))Then
	  Deallocate(bufferArray)
	End If
  End Subroutine calcConfigEnergies


!------------------------------------------------------------------------!
! ***calcEquilibriumVolume*** 
!
!------------------------------------------------------------------------!

Subroutine calcEquilibriumVolume()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n	
	Real(kind=DoubleReal) :: equilibriumVolume  
	Real(kind=DoubleReal), Dimension(1:2) :: equArray
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: outputArray
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
	Integer(kind=StandardInteger) :: tagA,tagB
    Integer(kind=StandardInteger) :: processTo,processFrom,bufferArraySize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArrayA, bufferArrayB
!Loop through the configurations set in the confiAtomsMap array
!Only calculate bulk modulus where there is a reference BM
    Call synchMpiProcesses() !Synchronise mpi processes at this point
    Do i=1,configCount
      If(mpiProcessID.eq.configAtomsMap(i,24))Then
		equArray = CalcEquilibrium(i)
		configurationEquVolume(i) = equArray(1)
		configurationEquLat(i) = equArray(2)
	  End If
    End Do
	Call synchMpiProcesses() !Synchronise mpi processes at this point
	Call distributeData1DDP(configurationEquVolume,size(configurationEquVolume,1),24)
	Call distributeData1DDP(configurationEquLat,size(configurationEquLat,1),24)
  End Subroutine calcEquilibriumVolume    
  
  
!------------------------------------------------------------------------!
! ***calcBulkModulii*** 
!
!------------------------------------------------------------------------!

Subroutine calcBulkModulii()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n	
	Real(kind=DoubleReal) :: bulkModulus  
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,bufferArraySize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray	
!Loop through the configurations set in the confiAtomsMap array
!Only calculate bulk modulus where there is a reference BM
    Call synchMpiProcesses() !Synchronise mpi processes at this point
    Do i=1,configCount
      If(mpiProcessID.eq.configAtomsMap(i,22))Then
		configurationBM(i) = CalcBulkModulus(i)
	  End If
    End Do
	Call synchMpiProcesses() !Synchronise mpi processes at this point	
	Call distributeData1DDP(configurationBM,size(configurationBM,1),22)
  End Subroutine calcBulkModulii  
  
  
  
  
  
!------------------------------------------------------------------------!
! ***calcCubicElasticConstants*** 
!
!------------------------------------------------------------------------!

  Subroutine calcCubicElasticConstants()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n	
	Real(kind=DoubleReal), Dimension(1:3) :: cubicElasticConstants  
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom,bufferArraySize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray	
	!Real(kind=DoubleReal), Dimension( : ), Allocatable :: sendArray	
	!Real(kind=DoubleReal), Dimension(1 : 12) :: sendArray
!Loop through the configurations set in the confiAtomsMap array
!Only calculate bulk modulus where there is a reference BM
    Call synchMpiProcesses() !Synchronise mpi processes at this point
    Do i=1,configCount
      If(mpiProcessID.eq.configAtomsMap(i,23))Then
		cubicElasticConstants = CalcCubicElasticConstant(i)
	    Do j=1,3
	      configurationEC(i,j) = cubicElasticConstants(j)
	    End Do
	  End If
    End Do
	Call synchMpiProcesses() !Synchronise mpi processes at this point
!Distribute data in array, 2D array, double precision
    Call distributeData2DDP(configurationEC,size(configurationEC,1),size(configurationEC,2),23)	
  End Subroutine calcCubicElasticConstants  
  
  
  
  
  

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
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  write(999,"(A44,F8.4)") "Calculating configuration bulk properties ",ProgramTime()
	End If
  
  
  End Subroutine calcConfigBulkProperties  
    
	
	
	
	
  
  
  
  
  
!------------------------------------------------------------------------!
! ***calcDifference*** 
!------------------------------------------------------------------------!
    
  Subroutine calcDifference()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: configCounter, atomCounter
	Real(kind=DoubleReal) :: energyDifference, forceDifference, trialRSS
	Real(kind=DoubleReal) :: volumeDifference
	Real(kind=DoubleReal) :: bmDifference, tempRSS, configForceDifference
	Real(kind=DoubleReal) :: calcEnergy
	Real(kind=DoubleReal) :: maxValue
	Logical :: printOut
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag
!Total residual square sum
    trialResidualSquareSum = 0.0D0  
!zero rss array
    Do i=1,configCount
	  Do j=1,10
	    configurationRSS(i,j) = 0.0D0
	  End Do	
	End Do
!-------------------------------------
!Energy Difference
!-------------------------------------
    energyDifference = 0.0D0
    Do i=1,configCount
!Check if there is a reference energy
      If(configAtomsMap(i,5).gt.0)Then
		calcEnergy = 1.0D0*(configurationEnergy(i)/configAtoms(i))
	    tempRSS = 1.0D0*(configurationRefEnergy(i)-calcEnergy)**2
		configurationRSS(i,1) = configurationRSS(i,1) + tempRSS
	  End If	
	End Do  
!Apply weighting and store
	trialResidualSquareSum = energyDifference
	rssEnergyDifference = energyDifference
!-------------------------------------
!Stress Difference
!-------------------------------------

!-------------------------------------
!Force Difference
!-------------------------------------
    forceDifference = 0.0D0
    If(calcForcesOnOff.eq.1)Then
	  configCounter = 1
	  atomCounter = 0
	  configForceDifference = 0.0D0
  	  Do i=1,configAtomsTotal 
	    atomCounter = atomCounter + 1
		If(configAtomsMap(configCounter,4).gt.0)Then  
		  tempRSS = 1.0D0*(configurationForceX(i)-configurationRefForceX(i))**2+&
			  1.0D0*(configurationForceY(i)-configurationRefForceY(i))**2+&
			  1.0D0*(configurationForceZ(i)-configurationRefForceZ(i))**2	
		  forceDifference=forceDifference+tempRSS  
		  configForceDifference = configForceDifference + tempRSS
		End If		  
		If(atomCounter.eq.configAtoms(configCounter))Then
		  atomCounter = 0
!Store config ev rss
		  configurationRSS(configCounter,3) = configurationRSS(configCounter,3) + configForceDifference
		  configCounter = configCounter + 1
		  configForceDifference = 0.0D0
		End If 		  
	  End Do
	  trialResidualSquareSum = trialResidualSquareSum + forceDifference
	  rssForceDifference = forceDifference
    End If
!-------------------------------------
!Equilibrium Volume
!-------------------------------------
    volumeDifference = 0.0D0
	Do i=1,configCount
	  If(configurationRefEquVolume(i).gt.-2.1D20.and.configurationEquVolume(i).gt.-2.1D20)Then
		tempRSS = (configurationEquVolume(i)-configurationRefEquVolume(i))**2
!Store config ev rss
		configurationRSS(i,4) = configurationRSS(i,4) + tempRSS
!sum total ev rss
		volumeDifference = volumeDifference + tempRSS
	  End If
    End Do	
	trialResidualSquareSum = trialResidualSquareSum + volumeDifference
!-------------------------------------
!Bulk Modulus Difference
!-------------------------------------
    bmDifference = 0.0D0
    Do i=1,configCount 
	  If(configurationBM(i).gt.(-2.0D20).and.configurationRefBM(i).gt.(-2.0D20))Then
	    bmDifference = bmDifference + 1.0D0*(configurationRefBM(i)-configurationBM(i))
		!print *,i,configurationRefBM(i),configurationBM(i)
	  End If
	End Do
!------------------------------------
! Normalise
!------------------------------------    
    Do j=1,4
	  maxValue = 0.0D0
	  Do i=1,configCount
        If(configurationRSS(i,j).ne.0.0D0)Then
          If(configurationRSS(i,j).gt.maxValue)Then
		    maxValue = configurationRSS(i,j)
		  End If
		End If
      End Do
	  If(maxValue.gt.0.0D0)Then
	    Do i=1,configCount
	      configurationRSS(i,j) = configurationRSS(i,j)/maxValue
        End Do
	  End If
	End Do  
!------------------------------------
! Apply config and global weights
!------------------------------------   
	Do j=1,4
	  Do i=1,configCount 	
	    configurationRSS(i,j) = 1.0D0*configurationOptWeights(i,j)*&
		  eamOptWeights(j)*configurationRSS(i,j)
	  End Do
    End Do
!------------------------------------
! Sum
!------------------------------------	
	trialResidualSquareSum = 0.0D0
	Do j=1,4
	  Do i=1,configCount
	    trialResidualSquareSum = trialResidualSquareSum + configurationRSS(i,j)
      End Do
	End Do 
  
  End Subroutine calcDifference
  
  
  

  

  
  
!------------------------------------------------------------------------!
! ***calcOutput*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!
  
  Subroutine calcOutputForces() 
!Outputs configuration forces to file
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k 
	Integer(kind=StandardInteger) :: configCounter, atomCounter
!If root/master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!Write forces to output file, if master process
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
!----------------------------------
  	End If  
  End Subroutine calcOutputForces
  
   
!------------------------------------------------------------------------!
! ***calcOutputDifference*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!
  
  Subroutine calcOutputDifference() 
!Outputs configuration calculation details
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k 
!Write energies to output file, if master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	 
!Write to file	  
	  write(999,"(A19,I4,A7,E20.10,A14,E20.10,A13,E20.10)") "Evaluation count: ",&
	  globalCounter(1),", RSS: ",trialResidualSquareSum,&
	  ", Energy RSS: ",rssEnergyDifference,&
	  ", Force RSS: ",rssForceDifference
	  
!Close output file
      close(999)
	End If  
  End Subroutine calcOutputDifference

!================================================================================================================================================  
!================================================================================================================================================
!------------------------------------------------------------------------!
! ***calcOutputTimes*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!  
  Subroutine calcOutputTimes() 
!Outputs configuration calculation details
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k 
!Write energies to output file, if master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	 
!Write to file	
	  If(globalTimer(1).gt.0.0D0)Then	    
	    write(999,"(A6,A41,E20.10,E20.10)") "      ",&
		"Calculation preparation:                 ",globalTimer(1),globalTimer(2)
	  End If
	  If(globalTimer(3).gt.0.0D0)Then	    
	    write(999,"(A6,A41,E20.10,E20.10)") "      ",&
		"Energy calculations:                     ",globalTimer(3),globalTimer(4)
	  End If
!Close output file
      close(999)
	End If  
  End Subroutine calcOutputTimes
  
!================================================================================================================================================  
!================================================================================================================================================
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
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: refStress  
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
	  If(neighbourListR(j).le.configurationRadiusCutoff(configurationID))Then	
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
	  Else
	    If(calcForcesOnOff.eq.1)Then
		  atomForceMagnitude(k) = 0.0D0
		  densityGrad(k,2) = 0.0D0 
		  densityGrad(k,1) = 0.0D0
		End If
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
		embedForce = 0.0D0
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
!Stress
        If(calcStressOnOff.eq.1)Then 
		  Do ii=1,3
		    Do jj=1,3 
              n = ((configurationID-1)*9)+3*(ii-1)+jj
			  configurationStress(n) = configurationStress(n)+&
			  neighbourListCoords(j,3+ii)*atomForceMagnitude(k)*&
			  neighbourListCoords(j,9+jj)
			End Do  
          End Do			
		End If
	  End Do
	End If	
!total energy of configuration
    energy = 0.5D0 * pairEnergy + embeddingEnergy   
!Calculate volume	
	!configurationVolume(configurationID) = calcConfigurationVolume(configurationID)
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

!================================================================================================================================================  
!================================================================================================================================================
  Function CalcEquilibrium (configurationID) RESULT (resultArray) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, n, configurationID, keyIJ, forceKey  
	Real(kind=DoubleReal) :: strainAmount, strainIncrement, configEnergyVal	
	Real(kind=DoubleReal) :: volume, equilibriumVolume, trialEquVolume, trialLat
	Real(kind=DoubleReal), Dimension(1:2) :: secondDerivative
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: strainArray, volEnergy, strainVolArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: volEnergyTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bmFit	
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: deriveData
	Real(kind=DoubleReal), Dimension(1:2) :: resultArray
!Allocate arrays
	Allocate(StrainArray(1:3,1:3))
	Allocate(volEnergy(1:7,1:2))
	Allocate(strainVolArray(1:7,1:2))
!-----------------------------------------------------
! Step 1 - save starting neighbourlist
!-----------------------------------------------------
!Store original neighbour list
    Call storeNeighboutList()
!-----------------------------------------------------
! Step 2 - first set of strains
!-----------------------------------------------------
!set increment
	strainIncrement = 0.01D0
!Loop over strains
    j = 0
	Do i=-3,3
	  j = j + 1
	  strainAmount = 1.0D0+strainIncrement*i
!Set up strain array
      StrainArray(1,1) = strainAmount
      StrainArray(1,2) = 0.0D0
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.0D0
      StrainArray(2,2) = strainAmount
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = strainAmount
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal	  
	  strainVolArray(j,1) = volume
	  strainVolArray(j,2) = strainAmount	  
    End Do
!Load original neighbour list
    Call loadNeighboutList()
!Fit polynomial
	coefficients = PolyFit(volEnergy,4)
!Calculate second derivative and minimum value
	deriveData = CalcSecondDerivativeMinimum(volEnergy, 3)
    trialEquVolume = deriveData(2)
	coefficients = PolyFit(strainVolArray,3)
	trialLat = ComputePolynomial(coefficients,trialEquVolume)
	volEnergyTemp = volEnergy
	Deallocate(volEnergy)
!-----------------------------------------------------
! Step 3 - second set of strains
!-----------------------------------------------------
!Allocate volEnergy array
	Allocate(volEnergy(1:14,1:2))
!Transfer data from first strain set
	Do i=1,7
	  volEnergy(i,1) = volEnergyTemp(i,1)
	  volEnergy(i,2) = volEnergyTemp(i,2)
	End Do
!set increment
	strainIncrement = 0.007D0
!Loop over strains
    j = 7
	Do i=-3,3
	  j = j + 1
	  strainAmount = trialLat+strainIncrement*i
!Set up strain array
      StrainArray(1,1) = strainAmount
      StrainArray(1,2) = 0.0D0
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.0D0
      StrainArray(2,2) = strainAmount
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = strainAmount
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal	
	End Do
!Load original neighbour list
    Call loadNeighboutList()
!Fit polynomial
	coefficients = PolyFit(volEnergy,4)
!Calculate second derivative and minimum value
	deriveData = CalcSecondDerivativeMinimum(volEnergy, 3)
    trialEquVolume = deriveData(2)
	coefficients = PolyFit(strainVolArray,3)
	trialLat = ComputePolynomial(coefficients,trialEquVolume)
	volEnergyTemp = volEnergy
	Deallocate(volEnergy)	
!-----------------------------------------------------
! Step 4 - third set of strains
!-----------------------------------------------------
!Allocate volEnergy array
	Allocate(volEnergy(1:21,1:2))
!Transfer data from first strain set
	Do i=1,14
	  volEnergy(i,1) = volEnergyTemp(i,1)
	  volEnergy(i,2) = volEnergyTemp(i,2)
	End Do
!set increment
	strainIncrement = 0.003D0
!Loop over strains
    j = 14
	Do i=-3,3
	  j = j + 1
	  strainAmount = trialLat+strainIncrement*i
!Set up strain array
      StrainArray(1,1) = strainAmount
      StrainArray(1,2) = 0.0D0
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.0D0
      StrainArray(2,2) = strainAmount
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = strainAmount
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal	
	End Do
!Load original neighbour list
    Call loadNeighboutList()
!Fit polynomial
	coefficients = PolyFit(volEnergy,4)
!Calculate second derivative and minimum value
	deriveData = CalcSecondDerivativeMinimum(volEnergy, 3)
    trialEquVolume = deriveData(2)
	coefficients = PolyFit(strainVolArray,3)
	trialLat = ComputePolynomial(coefficients,trialEquVolume)
!-----------------------------------------------------
! Step 5 - store equib volume and lattice multiplier
!-----------------------------------------------------	
	!configurationEquVolume(configurationID) = 1.0D0 * trialEquVolume
	!configurationEquLat(configurationID) = 1.0D0 * trialLat	
	resultArray(1) = 1.0D0 * trialEquVolume	
	resultArray(2) = 1.0D0 * trialLat
  End function CalcEquilibrium    

!================================================================================================================================================  
!================================================================================================================================================
  Function CalcBulkModulus (configurationID) RESULT (bulkModulus) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, n, configurationID, keyIJ, forceKey  
	Real(kind=DoubleReal) :: strainAmount, strainIncrement, configEnergyVal	
	Real(kind=DoubleReal) :: volume, bulkModulus, equVolume, equLat
	Real(kind=DoubleReal), Dimension(1:2) :: equArray
	Real(kind=DoubleReal), Dimension(1:2) :: secondDerivative
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: strainArray, volEnergy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bmFit	
!Allocate arrays
	Allocate(StrainArray(1:3,1:3))
	Allocate(volEnergy(1:9,1:2))	
!Store original neighbour list
    Call storeNeighboutList()
!Check if equilibrium has been calculated	
    If(configurationEquVolume(configurationID).lt.-2.0D20)Then
	  equArray = CalcEquilibrium(configurationID)
      configurationEquVolume(configurationID) = equArray(1)
      configurationEquLat(configurationID) = equArray(2)
	End If
!Load equilibrium volume and lattice parameter multiplier
	equVolume = 1.0D0 * configurationEquVolume(configurationID)
	equLat = 1.0D0 * configurationEquLat(configurationID)	
	If(equLat.eq.0.0D0)Then
	  equLat = 1.0D0
	End If
!set strain increment
	strainIncrement = 0.01D0
!Loop over strains
    j = 0
	Do i=-4,4
	  j = j + 1
	  strainAmount = equLat+strainIncrement*i
!Set up strain array
      StrainArray(1,1) = strainAmount
      StrainArray(1,2) = 0.0D0
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.0D0
      StrainArray(2,2) = strainAmount
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = strainAmount
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  If(i.eq.0.and.equVolume.eq.0.0D0)Then
	    equVolume = volume
	  End If
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal
    End Do
!Load original neighbour list
    Call loadNeighboutList()
!calculate second derivative d2E/dV2	
	secondDerivative = CalcSecondDerivativeMinimum(volEnergy)
	bulkModulus = secondDerivative(2)*secondDerivative(1)
	bulkModulus = bulkModulus*elementaryCharge*1.0D30*1.0D-9	
  End function CalcBulkModulus
!================================================================================================================================================  
!================================================================================================================================================
  Function CalcCubicElasticConstant (configurationID) RESULT (elasticConstants) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, n, configurationID
	Real(kind=DoubleReal) :: strainAmount, strainIncrement, configEnergyVal
	Real(kind=DoubleReal) :: strainAmountA, strainAmountB, strainAmountC
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: strainArray, volEnergy, strainEnergy	
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients	
	Real(kind=DoubleReal) :: orthorhombicStrainModulus, monoclinicStrainModulus
	Real(kind=DoubleReal) :: tetragonalStrainModulus
	
	Real(kind=DoubleReal) :: volume, bulkModulus, equVolume, equLat
	Real(kind=DoubleReal), Dimension(1:2) :: equArray
	Real(kind=DoubleReal), Dimension(1:2) :: secondDerivative
	Real(kind=DoubleReal), Dimension(1:3) :: elasticConstants
!Allocate arrays
	Allocate(StrainArray(1:3,1:3))
	Allocate(volEnergy(1:6,1:2))	
	Allocate(strainEnergy(1:6,1:2))		
!Check if equilibrium has been calculated	
    If(configurationEquVolume(configurationID).lt.-2.0D20)Then
	  equArray = CalcEquilibrium(configurationID)
      configurationEquVolume(configurationID) = equArray(1)
      configurationEquLat(configurationID) = equArray(2)
	End If	
!Check if bulk modulus has been calculated	
    If(configurationBM(configurationID).lt.-2.0D20)Then
      configurationBM(configurationID) = CalcBulkModulus(configurationID)
	End If	
!Load equilibrium volume and lattice parameter multiplier
	equVolume = 1.0D0 * configurationEquVolume(configurationID)
	equLat = 1.0D0 * configurationEquLat(configurationID)	
	If(equLat.eq.0.0D0)Then
	  equLat = 1.0D0
	End If
	bulkModulus = 1.0D0 * configurationBM(configurationID)	
!---------------------------------------------------------------------	
!Tetragonal Strain
!---------------------------------------------------------------------	
!set strain increment
	strainIncrement = 0.01D0
!Loop over strains
    j = 0
	Do i=0,5
	  j = j + 1
	  !strainAmount = equLat+strainIncrement*i
	  strainAmountA = equLat+strainIncrement*i
	  strainAmountB = equLat+strainIncrement*i
	  strainAmountC = equLat+(1+(strainIncrement*i))**(-2)-1
!Set up strain array
      StrainArray(1,1) = strainAmountA
      StrainArray(1,2) = 0.0D0
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.0D0
      StrainArray(2,2) = strainAmountB
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = strainAmountC
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  If(i.eq.0.and.equVolume.eq.0.0D0)Then
	    equVolume = volume
	  End If
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal
	  strainEnergy(j,1) = strainAmountA
	  strainEnergy(j,2) = configEnergyVal
    End Do
!Fit polynomial
	coefficients = PolyFit(strainEnergy,2)
	tetragonalStrainModulus = coefficients(3)/(3.0D0*equVolume)
	tetragonalStrainModulus = tetragonalStrainModulus*elementaryCharge*1.0D30*1.0D-9	
!---------------------------------------------------------------------	
!Orthorhombic Strain
!---------------------------------------------------------------------	
!set strain increment
	strainIncrement = 0.01D0
!Loop over strains
    j = 0
	Do i=0,5
	  j = j + 1
	  !strainAmount = equLat+strainIncrement*i
	  strainAmountA = equLat+strainIncrement*i
	  strainAmountB = equLat-strainIncrement*i
	  strainAmountC = equLat+((strainIncrement*i)**2/(1-(strainIncrement*i)**2))
!Set up strain array
      StrainArray(1,1) = strainAmountA
      StrainArray(1,2) = 0.0D0
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.0D0
      StrainArray(2,2) = strainAmountB
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = strainAmountC
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  If(i.eq.0.and.equVolume.eq.0.0D0)Then
	    equVolume = volume
	  End If
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal
	  strainEnergy(j,1) = strainAmountA
	  strainEnergy(j,2) = configEnergyVal
    End Do
!Fit polynomial
	coefficients = PolyFit(strainEnergy,2)
	orthorhombicStrainModulus = coefficients(3)/equVolume
	orthorhombicStrainModulus = orthorhombicStrainModulus*elementaryCharge*1.0D30*1.0D-9
!---------------------------------------------------------------------	
!Monoclinic Strain
!---------------------------------------------------------------------	
!set strain increment
	strainIncrement = 0.002D0
!Loop over strains
    j = 0
	Do i=0,5
	  j = j + 1
	  strainAmountA = strainIncrement*i
	  strainAmountB = equLat+((strainIncrement*i)**2/(1-(strainIncrement*i)**2))
!Set up strain array
      StrainArray(1,1) = 1.0D0
      StrainArray(1,2) = 0.5D0*strainIncrement*i
      StrainArray(1,3) = 0.0D0
      StrainArray(2,1) = 0.5D0*strainIncrement*i
      StrainArray(2,2) = 1.0D0
      StrainArray(2,3) = 0.0D0
      StrainArray(3,1) = 0.0D0
      StrainArray(3,2) = 0.0D0
      StrainArray(3,3) = 1.0D0+((strainIncrement*i)**2)/(4-(strainIncrement*i)**2)
!Load original neighbour list
      Call loadNeighboutList()
!Multiply apply distortion to the neighbour list
	  Call applyDistortionVector(configurationID,StrainArray)
!Calculate energy
      configEnergyVal = CalcEnergy(configurationID)
	  volume = calcConfigurationVolume(configurationID,StrainArray)
	  If(i.eq.0.and.equVolume.eq.0.0D0)Then
	    equVolume = volume
	  End If
	  volEnergy(j,1) = volume
	  volEnergy(j,2) = configEnergyVal
	  strainEnergy(j,1) = strainAmountA
	  strainEnergy(j,2) = configEnergyVal
    End Do
!Fit polynomial
	coefficients = PolyFit(strainEnergy,2)
	monoclinicStrainModulus = 2.0D0*coefficients(3)/equVolume
	monoclinicStrainModulus = monoclinicStrainModulus*elementaryCharge*1.0D30*1.0D-9	
!---------------------------------------------------------------------	
!Set elastic constants
!---------------------------------------------------------------------		
	elasticConstants(2) = (3.0D0*bulkModulus-1.0D0*orthorhombicStrainModulus)/3.0D0	!C12
	elasticConstants(1) = 1.0D0*orthorhombicStrainModulus+1.0D0*elasticConstants(2) !C11
	elasticConstants(3) = 1.0D0*monoclinicStrainModulus
  End function CalcCubicElasticConstant    
!================================================================================================================================================  
!================================================================================================================================================
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
	If(eamDataSet.eq.1)Then		!Original input EAM potential
!get the start point and length of potential data
      potStart = eamKey(potKey,4)
      potLength = eamKey(potKey,5)
!Return 0.0D0 if out of range
      If(x.ge.eamData(potStart,1).and.x.le.eamData((potStart+potLength-1),1))Then
        interpPointCount = eamInterpType + 1
!interpolate
	    yArray = PointInterpolationArr(eamData,x,interpPointCount,potStart,potLength)
	  Else
	    yArray(1) = 0.0D0
	    yArray(2) = 0.0D0
	  End If
	End If
	If(eamDataSet.eq.2)Then		!Reduced data set
!get the start point and length of potential data
      potStart = eamKeyReduced(potKey,4)
      potLength = eamKeyReduced(potKey,5)
	  !Return 0.0D0 if out of range
      If(x.ge.eamDataReduced(potStart,1).and.&
	  x.le.eamDataReduced((potStart+potLength-1),1))Then
        interpPointCount = eamInterpType + 1
!interpolate
	    yArray = PointInterpolationArr(eamDataReduced,x,interpPointCount,&
		potStart,potLength)
	  Else
	    yArray(1) = 0.0D0
	    yArray(2) = 0.0D0
	  End If
	End If
	If(eamDataSet.eq.3)Then		!Trial new potential
!get the start point and length of potential data
      potStart = eamKeyTrial(potKey,4)
      potLength = eamKeyTrial(potKey,5)
!Return 0.0D0 if out of range
      If(x.ge.eamDataTrial(potStart,1).and.&
	  x.le.eamDataTrial((potStart+potLength-1),1))Then
        interpPointCount = eamInterpType + 1
!interpolate
	    yArray = PointInterpolationArr(eamDataTrial,x,interpPointCount,&
		potStart,potLength)
	  Else
	    yArray(1) = 0.0D0
	    yArray(2) = 0.0D0
	  End If
	End If
  End function SearchPotentialPoint    
  
  
  

!================================================================================================================================================  
!================================================================================================================================================  
! MPI Data Distribution Subroutines
!================================================================================================================================================  
!================================================================================================================================================  
!------------------------------------------------------------------------!
! ***distributeData1D*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!   
  Subroutine distributeData1DDP(dataArray,arraySize,mapColumn) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: arraySize
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount,mapColumn
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension(1:arraySize) :: dataArray, sendArray, recvArray	
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
	Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)	
!-----------------------------------------	
! Collect array from workers to master
!-----------------------------------------    
    If(processID.eq.0)Then
!RECV by master process	  
	  Do i=1,(processCount-1)
	    processFrom = i
		tag = 1000 + i
        Call MPI_recv(recvArray,arraySize,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
!Read specific worker slots into master array
        Do j=1,configCount
		  If(configAtomsMap(j,mapColumn).eq.i)Then
		    dataArray(j) = recvArray(j)
		  End If
		End Do
	  End Do
	Else
!SEND by worker processes
	  processTo = 0
	  sendArray = dataArray
      Do i=1,(processCount-1)
        If(i.eq.processID)Then
          tag = 1000 + i	
          Call MPI_send(sendArray,arraySize,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End If
      End Do
	End If
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
	Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)	
!-----------------------------------------	
! Send array from master to workers
!-----------------------------------------   
    If(processID.eq.0)Then
!SEND by master process	 
	  sendArray = dataArray
      Do i=1,(processCount-1)
	    processTo = i
        tag = 2000 + i
		Call MPI_send(sendArray,arraySize,&
		MPI_double_precision,processTo,tag,MPI_comm_world,error)
	  End Do	
	Else
!RECV by worker processes
      processFrom = 0
	  tag = 2000 + processID
      Call MPI_recv(recvArray,arraySize,&
	  MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
	  dataArray = recvArray
	End If    
  End Subroutine distributeData1DDP  
!------------------------------------------------------------------------!
! ***distributeData1D*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine distributeData2DDP(dataArray,arraySizeH,arraySizeW,mapColumn) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n
	Integer(kind=StandardInteger) :: arraySizeH,arraySizeW
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount,mapColumn
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension(1:arraySizeH,1:arraySizeW) :: dataArray
	Real(kind=DoubleReal), Dimension(1:arraySizeH) :: dataArrayTemp, sendArray, recvArray
!Loop array columns
	Do n=1,arraySizeW
!-----------------------------------------
! Make temp 1D array
!-----------------------------------------	
	  Do i=1,arraySizeH
	    dataArrayTemp(i) = dataArray(i,n)
	  End Do
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
	  Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
      Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)	
!-----------------------------------------	
! Collect array from workers to master
!-----------------------------------------    
      If(processID.eq.0)Then
!RECV by master process	  
	    Do i=1,(processCount-1)
	      processFrom = i
		  tag = 1000 + i
          Call MPI_recv(recvArray,arraySizeH,&
		  MPI_double_precision,processFrom,tag,MPI_comm_world,status,error )
!Read specific worker slots into master array
          Do j=1,configCount
		    If(configAtomsMap(j,mapColumn).eq.i)Then
		      dataArrayTemp(j) = recvArray(j)
		    End If
		  End Do
	    End Do
	  Else
!SEND by worker processes
	    processTo = 0
	    sendArray = dataArrayTemp
        Do i=1,(processCount-1)
          If(i.eq.processID)Then
            tag = 1000 + i	
            Call MPI_send(sendArray,arraySizeH,&
		    MPI_double_precision,processTo,tag,MPI_comm_world,error)
	      End If
        End Do
	  End If
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
	  Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
      Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)	
!-----------------------------------------	
! Send array from master to workers
!-----------------------------------------   
      If(processID.eq.0)Then
!SEND by master process	 
	    sendArray = dataArrayTemp
        Do i=1,(processCount-1)
	      processTo = i
          tag = 2000 + i
		  Call MPI_send(sendArray,arraySizeH,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End Do	
	  Else
!RECV by worker processes
        processFrom = 0
	    tag = 2000 + processID
        Call MPI_recv(recvArray,arraySizeH,&
	    MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
	    dataArrayTemp = recvArray
	  End If    
!-----------------------------------------
! Merge 1D array with 2D array
!-----------------------------------------		  
      Do i=1,arraySizeH
	    dataArray(i,n) = dataArrayTemp(i)
	  End Do
	End Do
  End Subroutine distributeData2DDP
!------------------------------------------------------------------------!
! ***sendData1D*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine sendData1DDP(dataArray,arraySize) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: arraySize
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension(1:arraySize) :: dataArray, sendArray, recvArray	
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
	Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)	
!-----------------------------------------	
! Send array from master to workers
!-----------------------------------------   
    If(processID.eq.0)Then
!SEND by master process	 
	  sendArray = dataArray
      Do i=1,(processCount-1)
	    processTo = i
        tag = 2000 + i
		Call MPI_send(sendArray,arraySize,&
		MPI_double_precision,processTo,tag,MPI_comm_world,error)
	  End Do	
	Else
!RECV by worker processes
      processFrom = 0
	  tag = 2000 + processID
      Call MPI_recv(recvArray,arraySize,&
	  MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
	  dataArray = recvArray
	End If    
  End Subroutine sendData1DDP
  

End Module calc