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
	Real(kind=DoubleReal) :: rss
	Real(kind=DoubleReal), Dimension(1:1024,1:10) :: rssTable
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
    Call calcDifference(rss,rssTable)
!Save difference to file
	Call calcOutputDifference() 
  End Subroutine calcEval
  
  
  
  Subroutine calcEvalFull()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Real(kind=DoubleReal) :: rss
	Real(kind=DoubleReal), Dimension(1:1024,1:10) :: rssTable
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
    Call calcDifference(rss,rssTable)
	!print *,configurationRSS(1,4)
	!print *,configurationRSS(1,5)
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
	  !outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  !open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
	  !write(999,"(A33,F8.4,A3,I4)") "Calculating configuration energy ",&
	  !ProgramTime(),"   ",globalCounter(1)
	End If
!zero out arrays
	configurationEnergy = 0.0D0
	configurationForce = 0.0D0
	configurationStresses = 0.0D0
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
    Call MPI_sumData2DDP(configurationForce,&
	size(configurationForce,1),size(configurationForce,2))
    Call MPI_sumData2DDP(configurationStresses,&
	size(configurationStresses,1),size(configurationStresses,2))
!Save forces to file
    If(saveFileForces.eq."Y")Then
	  Call outputForces()
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
    
  Subroutine calcDifference(rss, rssTable)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n
	Integer(kind=StandardInteger) :: configCounter, atomCounter, configAtomsCount
	Real(kind=DoubleReal) :: energyDifference, forceDifference, trialRSS
	Real(kind=DoubleReal) :: volumeDifference
	Real(kind=DoubleReal) :: bmDifference, tempRSS, configForceDifference
	Real(kind=DoubleReal) :: calcEnergy
	Real(kind=DoubleReal) :: maxValue
	Real(kind=DoubleReal), Dimension(1:10) :: normMultipliers
	Real(kind=DoubleReal) :: rss, configRss
	Real(kind=DoubleReal), Dimension(1:1024,1:10) :: rssTable
	Logical :: printOut
!mpi variables
    Integer(kind=StandardInteger) :: selectProcess,status,error,tag
!Total residual square sum
    trialResidualSquareSum = 0.0D0  
!zero rss array
    rssTable = 0.0D0
!-------------------------------------
!Energy Difference
!-------------------------------------
    Do i=1,configCount
!Check if there is a reference energy
      If(configAtomsMap(i,5).gt.0)Then
		calcEnergy = 1.0D0*(configurationEnergy(i)/configAtoms(i))
		rssTable(i,1) = 1.0D0*(configurationRefEnergy(i)-calcEnergy)**2
		!rss energy per atom
	  End If	
	End Do  
!-------------------------------------
!Stress Difference
!-------------------------------------
    Do i=1,configCount	
	  If(configurationRefStresses(i,1).gt.-2.1D20)Then
	    Do j=1,9
		  rssTable(i,2) = rssTable(i,2) + &
		    1.0D0*((configurationRefStresses(i,j)-configurationStresses(i,j))&
			    /9.0D0)**2
		  !per stress tensor component
		End Do
      End If
	End Do
!-------------------------------------
!Force Difference
!-------------------------------------
    If(calcForcesOnOff.eq.1)Then
	  configCounter = 1
	  atomCounter = 0
	  configForceDifference = 0.0D0
!Loop through atoms
	  Do i=1,configCount
	    configAtomsCount = configAtoms(i)
	    Do j=1,configAtomsCount
	      n = configAtomsStart(i) + j		  
		  If(configurationRefForce(n,1).gt.-2.0D20)Then
		    Do k=1,3
		      rssTable(i,3) = rssTable(i,3) + &
		      1.0D0*((configurationRefForce(n,k)-configurationForce(n,k))/&
			  (configAtomsCount*3.0D0))**2	!rss avg force per atom per direction
		    End Do
		  End If
	    End Do
	  End Do	  
    End If
!-------------------------------------
!Equilibrium Volume
!-------------------------------------
	Do i=1,configCount
	  configAtomsCount = configAtoms(i)
	  If(configurationRefEquVolume(i).gt.-2.1D20.and.configurationEquVolume(i).gt.-2.1D20)Then
!Store config ev rss
		rssTable(i,4) = rssTable(i,4) + &
		    1.0D0*((configurationEquVolume(i)-configurationRefEquVolume(i))/&
			configAtomsCount)**2 !rss volume per atom
	  End If
    End Do	
!-------------------------------------
!Bulk Modulus Difference
!-------------------------------------
    bmDifference = 0.0D0
    Do i=1,configCount 
	  If(configurationBM(i).gt.(-2.0D20).and.configurationRefBM(i).gt.(-2.0D20))Then
	    rssTable(i,5) = rssTable(i,5) + &
		    1.0D0*(configurationRefBM(i)-configurationBM(i))**2
	  End If
	End Do
!------------------------------------
! Normalise
!------------------------------------  
	normMultipliers = 1.0D0
	normMultipliers(1) = (1.0D-1)**2		!energy normal per atom eV
	normMultipliers(2) = (1.0D-2)**2		!stress normal per direction gpa
	normMultipliers(3) = (1.0D-4)**2		!force normal per atom eV/ang
	normMultipliers(4) = (0.1D0)**2		    !vol normal per atom ang3
	normMultipliers(5) = (5.0D0)**2		    !bm normal per atom gpa
    Do i=1,configCount 
	  Do j=1,10 
	    rssTable(i,j) = rssTable(i,j) / normMultipliers(j)
	  End Do
	End Do 
!------------------------------------
! Apply config and global weights
!------------------------------------   
	Do j=1,10
	  Do i=1,configCount 
        If(configurationOptWeights(i,j).gt.0)Then
	      rssTable(i,j) = 1.0D0*configurationOptWeights(i,j)*&
		  eamOptWeights(j)*rssTable(i,j)
		End If
	  End Do
    End Do
!------------------------------------
! Sum
!------------------------------------	
	rss = 0.0D0
	Do j=1,10
	  Do i=1,configCount
	    rss = rss + rssTable(i,j)
      End Do
	End Do 
    configurationRSS = rssTable
	trialResidualSquareSum = rss
!------------------------------------
! Record rss
!------------------------------------	
	
!Create output file
    If(mpiProcessID.eq.0)Then
      globalCounter(3) = globalCounter(3) + 1
	  open(unit=979,file=trim(outputDirectory)//"/"//"rssLog.dat",&
	  status="old",position="append",action="write")	
	  If(globalCounter(3).eq.1)Then
	    write(979,"(A88)") &
		"========================================================================================"
	    write(979,"(A88)") &
		"n     Config    RSS         Energy      Stress      Forces      Eq.Vol      BM          "
	    write(979,"(A88)") &
		"========================================================================================"
	  End If
	  Do i=1,configCount
	    configRss = 0.0D0
	    Do j=1,10
	      configRss = configRss + rssTable(i,j)
	    End Do
	    write(979,"(I5,A1,I9,A1,F11.5,A1,F11.5,A1,F11.5,A1,F11.5,A1,F11.5,A1,F11.5,A1)")&
		globalCounter(2)," ",i," ",configRss," ",&
		rssTable(i,1)," ",rssTable(i,2)," ",rssTable(i,3)," ",rssTable(i,4)," ",rssTable(i,5)," "		
	  End Do 
      write(979,"(A17,F11.5)") "Calculation RSS: ",rss
	  write(979,"(A88)") &
	  "----------------------------------------------------------------------------------------"
!close output file
	  close(979)
	End If
	
	
  End Subroutine calcDifference
  
  
  

  

  
  
!------------------------------------------------------------------------!
! ***calcOutput*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!
  

  
   
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
	  !write(999,"(A19,I4,A7,E20.10,A14,E20.10,A13,E20.10)") "Evaluation count: ",&
	  !globalCounter(1),", RSS: ",trialResidualSquareSum,&
	  !", Energy RSS: ",rssEnergyDifference,&
	  !", Force RSS: ",rssForceDifference
	  
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
    Integer(kind=StandardInteger) :: configurationID
	Real(kind=DoubleReal) :: energy
	Integer(kind=StandardInteger) :: i, j, k, l, n, m
	Integer(kind=StandardInteger) :: configStart, configLength, configAtomsCount
	Integer(kind=StandardInteger) :: densityCount, forceKeyA, forceKeyB
!Data arrays
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
    Real(kind=DoubleReal), Dimension(1:500000) :: density
!Data variables	
	Real(kind=DoubleReal) :: r, pairEnergy, embeddingEnergy
	Real(kind=DoubleReal) :: forceM, forceS, forceX, forceY, forceZ
	Real(kind=DoubleReal), Dimension(1:3) :: forceArr, forceA
	Real(kind=DoubleReal) :: densDerivAB, densDerivBA, embeDerivA, embeDerivB
!config values
	configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)
    configAtomsCount = configAtoms(configurationID)
!Set starting values
    energy = 0.0D0	
    pairEnergy = 0.0D0
    density	= 0.0D0
	embeddingEnergy	= 0.0D0	
!Zero forces for this configuration
    Do i=1,configAtomsCount
	  n = configAtomsStart(configurationID) + i
	  configurationForce(n,1) = 0.0D0
	  configurationForce(n,2) = 0.0D0
	  configurationForce(n,3) = 0.0D0
	End Do
!
! Neighbour list information
! neighbourListR(n)    Atom separation
! neighbourListI(n,1)  Atom A type
! neighbourListI(n,2)  Atom B type
! neighbourListI(n,3)  Atom A id
! neighbourListI(n,4)  Atom B id
!
!--------------------------------------------------
! Loop 1 - sum pair energy and density	
!--------------------------------------------------
    k = 0
	Do n=configStart,configStart+configLength-1 		!loop through neighbour list for configuration
	  k = k + 1
!If within cutoff radius
      If(neighbourListR(n).le.configurationRadiusCutoff(configurationID))Then
!Pot Energy: Sum pair potential energy between A and all the Bs around it
	    yArray = SearchPotentialPoint(neighbourListI(n,1),&			!Returns yArray(1) value of function at r [neighbourListR(n)]
		         neighbourListI(n,2),1,neighbourListR(n))           !yArray(2) derivative of function at r [neighbourListR(n)]
	    pairEnergy = pairEnergy+yArray(1)
!Forces: pair potential component of force	 
        forceArr(1) = yArray(2)*neighbourListCoords(n,10)	!x component direction and force from atom A to B
        forceArr(2) = yArray(2)*neighbourListCoords(n,11)	!y component direction and force from atom A to B
        forceArr(3) = yArray(2)*neighbourListCoords(n,12)	!z component direction and force from atom A to B
!Force key
        forceKeyA = configAtomsStart(configurationID) + neighbourListI(n,3)
        forceKeyB = configAtomsStart(configurationID) + neighbourListI(n,4)
!Force on atom A		
		configurationForce(forceKeyA,1)=configurationForce(forceKeyA,1)-1.0D0*forceArr(1)
		configurationForce(forceKeyA,2)=configurationForce(forceKeyA,2)-1.0D0*forceArr(2)
		configurationForce(forceKeyA,3)=configurationForce(forceKeyA,3)-1.0D0*forceArr(3)
!Force on atom B		
		configurationForce(forceKeyB,1)=configurationForce(forceKeyB,1)+1.0D0*forceArr(1)
		configurationForce(forceKeyB,2)=configurationForce(forceKeyB,2)+1.0D0*forceArr(2)
		configurationForce(forceKeyB,3)=configurationForce(forceKeyB,3)+1.0D0*forceArr(3)
!Pot Energy: Electron density of each A due to the electrons of the Bs around it	
        yArray = SearchPotentialPoint(neighbourListI(n,2),0,2,neighbourListR(n))	
		density(neighbourListI(n,3)) = density(neighbourListI(n,3)) + yArray(1)
!Pot Energy: Electron density of each B due to the electrons of the As around it			
        yArray = SearchPotentialPoint(neighbourListI(n,1),0,2,neighbourListR(n))	
		density(neighbourListI(n,4)) = density(neighbourListI(n,4)) + yArray(1)
!Stress from pair contribution
		If(neighbourListI(n,6).eq.1)Then
	      Do i=1,3   !coord i j k
		    Do j=1,3 !force on k due to l
		      k = (i-1)*3+j
			  configurationStresses(configurationID,k) = &
			  configurationStresses(configurationID,k) + &
			  (neighbourListCoords(n,i+3)-neighbourListCoords(n,i)) * forceArr(j)
			End Do
          End Do
        End If	
	  End If
	End Do
!--------------------------------------------------
! Loop 2 - embedding energy	
!--------------------------------------------------
    Do n=1,configAtomsCount		!loop through atoms in the configuration
!Pot Energy: Embedding energy	
	  m = configAtomsStart(configurationID) + n			!atom id
      yArray = SearchPotentialPoint(atomTypeKey(m),0,3,density(n))		!get energy to embed atom at that density
	  embeddingEnergy = embeddingEnergy+yArray(1)	
	End Do  
!--------------------------------------------------
! Loop 3 - remaining derivative values for force	
!--------------------------------------------------
    k = 0
	Do n=configStart,configStart+configLength-1 
	  k = k + 1	
!If within cutoff radius
      If(neighbourListR(n).le.configurationRadiusCutoff(configurationID))Then  
!Density derivatives
	    yArray = SearchPotentialPoint(neighbourListI(n,2),0,2,neighbourListR(n))
	    densDerivAB = yArray(2)		!derivative of the density
	    yArray = SearchPotentialPoint(neighbourListI(n,1),0,2,neighbourListR(n))
	    densDerivBA = yArray(2)
!Embedding derivatives
        yArray = SearchPotentialPoint(neighbourListI(n,1),0,3,density(neighbourListI(n,3)))
	    embeDerivA = yArray(2)
        yArray = SearchPotentialPoint(neighbourListI(n,2),0,3,density(neighbourListI(n,4)))
	    embeDerivB = yArray(2)
!Force, embedding dens A
	    forceM = (embeDerivA*densDerivAB+embeDerivB*densDerivBA)
	    forceArr(1) = forceM*neighbourListCoords(n,10)
	    forceArr(2) = forceM*neighbourListCoords(n,11)
	    forceArr(3) = forceM*neighbourListCoords(n,12)
!Force key
        forceKeyA = configAtomsStart(configurationID) + neighbourListI(n,3)
        forceKeyB = configAtomsStart(configurationID) + neighbourListI(n,4)
!Force on atom A		
	    configurationForce(forceKeyA,1)=configurationForce(forceKeyA,1)-1.0D0*forceA(1)
	    configurationForce(forceKeyA,2)=configurationForce(forceKeyA,2)-1.0D0*forceA(2)
	    configurationForce(forceKeyA,3)=configurationForce(forceKeyA,3)-1.0D0*forceA(3)
!Force on atom B		
	    configurationForce(forceKeyB,1)=configurationForce(forceKeyB,1)+1.0D0*forceA(1)
	    configurationForce(forceKeyB,2)=configurationForce(forceKeyB,2)+1.0D0*forceA(2)
	    configurationForce(forceKeyB,3)=configurationForce(forceKeyB,3)+1.0D0*forceA(3)
!Stress from embedding contribution
		If(neighbourListI(n,6).eq.1)Then
	      Do i=1,3   !coord i j k
		    Do j=1,3 !force on k due to l
		      k = (i-1)*3+j
			  configurationStresses(configurationID,k) = &
			  configurationStresses(configurationID,k) + &
			  (neighbourListCoords(n,i+3)-neighbourListCoords(n,i)) * forceArr(j)
			End Do
          End Do
        End If	
	  End If	
	End Do
!total energy of configuration
    energy = pairEnergy+embeddingEnergy   	
!convert stress to gpa
    Do i=1,9
      configurationStresses(configurationID,i) = &
	    UnitConvert(configurationStresses(configurationID,i),"EVAN3","GPA")
	  configurationStresses(configurationID,i) = &
	    (1.0D0/(2.0D0*configurationVolume(configurationID)))*configurationStresses(configurationID,i)
	End Do 
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
	Real(kind=DoubleReal), Dimension(1:50,1:3) :: bmData
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bmFit	
!Allocate arrays
	Allocate(StrainArray(1:3,1:3))
	Allocate(volEnergy(1:9,1:2))	
    If(saveFileBM.eq."Y")Then	
	  bmData = -2.1D20
	End If
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
!bm data array
      If(saveFileBM.eq."Y")Then
        bmData(j,1) = 1.0D0*strainAmount
        bmData(j,2) = 1.0D0*volume
        bmData(j,3) = 1.0D0*configEnergyVal
	  End If
    End Do
!Load original neighbour list
    Call loadNeighboutList()
!calculate second derivative d2E/dV2	
	secondDerivative = CalcSecondDerivativeMinimum(volEnergy)
	bulkModulus = secondDerivative(2)*secondDerivative(1)
	bulkModulus = bulkModulus*elementaryCharge*1.0D30*1.0D-9		
!Output BM details
    If(saveFileBM.eq."Y")Then
	  Call outputBM(configurationID,bmData)	
	End If
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