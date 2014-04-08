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
  Public :: calcEvalFull		
  Public :: calcConfigEnergies
  Public :: calcDifference
  Public :: calcOutput
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
	Real(kind=DoubleReal) :: trialRSS
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
	Real(kind=DoubleReal) :: trialRSS
!Count energy calculation
    globalCounter(2) = globalCounter(2) + 1
!Calculate energies, forces and stresses
    Call calcConfigEnergies()
!Bulk modulus
    Call calcBulkModulii()
!Synch MPI processes
	Call synchMpiProcesses()
!Calc difference 
    Call calcDifference()
!Save difference to file
	Call calcOutputDifference() 
	eamFunctionWobbliness = eamWobbliness(eamKey, eamData, 10000)
	!wobbliness(potKey)
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
	Call MPI_sendout(configurationStress,"double")		
!Set end time	
    endTime = ProgramTime() 	
	!startTime
  End Subroutine calcConfigEnergies
  
  
  
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
    Do i=1,size(configAtomsMap,1)
      If(mpiProcessID.eq.configAtomsMap(i,22))Then
		configurationBM(i) = CalcBulkModulus(i)
	  End If
    End Do
	Allocate(bufferArray(1:configCount))
    If(mpiProcessID.gt.0) Then
!send buffers from all worker processes
      Do i=1,(mpiProcessCount-1)
        If(i.eq.mpiProcessID)Then
          processTo = 0
          tag = 1500 + mpiProcessID	
		  bufferArraySize = configCount
          Call MPI_send(configurationBM,bufferArraySize,&
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
        tag = 1500 + j
		bufferArraySize = configCount
        Call MPI_recv(bufferArray,bufferArraySize,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
		processFrom = j
		Do i=1,configCount
		  If(configAtomsMap(i,22).eq.processFrom)Then
		    configurationBM(i) = bufferArray(i)
		  End If
		End Do
      End Do
	End If  
!send array from master to worker processes
	Call MPI_sendout(configurationBM,"double")		
  End Subroutine calcBulkModulii  
  
  
  

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
	Real(kind=DoubleReal) :: bmDifference, tempRSS, configForceDifference
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
	    tempRSS = 1.0D0*(configurationEnergy(i)/configAtoms(i)-&
	    configurationRefEnergy(i))**2
	    energyDifference = energyDifference + tempRSS
		configurationRSS(i,1) = configurationRSS(i,1) + tempRSS
	  End If	
	End Do  
	trialResidualSquareSum = energyDifference*100.0D0
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
		  configurationRSS(configCounter,2) = &
		    configurationRSS(configCounter,2) + configForceDifference
		  configCounter = configCounter + 1
		  configForceDifference = 0.0D0
		End If 		  
	  End Do
	  trialResidualSquareSum = trialResidualSquareSum + forceDifference
    End If
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
	
	
	!Write forces to output file, if master process
	!If(mpiProcessID.eq.0.and.printOut.eqv..true.)Then	
!open forces output file	
	  !open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
      !write(999,"(A40,E20.10)") "      ----------------------------------"		  
      !write(999,"(A40,E20.10)") "      Potential vs Reference Data RSS"		  
      !write(999,"(A40,E20.10)") "      ----------------------------------"		  
	  !write(999,"(A31,E20.10)") "      RSS              :       ",trialRSS	  
	  !write(999,"(A31,E20.10)") "      Energy difference:       ",energyDifference	  
	  !write(999,"(A31,E20.10)") "      Force difference:        ",forceDifference
	  !write(999,"(A31,E20.10)") "      Bulk Modulus difference: ",bmDifference
	  !close(999)	
	!End If
	
	
  
  End Subroutine calcDifference
  
  
  

  
!------------------------------------------------------------------------!
! ***calcOutput*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!
  
  Subroutine calcOutput () 
!Outputs configuration calculation details
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k 
	Integer(kind=StandardInteger) :: configCounter, atomCounter  
!Write energies to output file, if master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	 
!Write to file	  
	  write(999,"(A44,I4)") "--------------------------------------------"
	  write(999,"(A19,I4)") "Calculation count: ",globalCounter(1)
	  write(999,"(A44,I4)") "--------------------------------------------"
	  Do i=1,configCount
	    write(999,"(I4,A1,F12.6,A4,I8,A10,I8)") i," ",&
	      (1.0*configurationEnergy(i))," eV ",&
	    configAtoms(i)," atoms     ",configAtomsMap(i,3)
		If(configAtomsMap(i,5).gt.0)Then
	      write(999,"(I4,F12.6,A14,F12.6,A8)") i,&
		  (configurationEnergy(i)/configAtoms(i))," eV per atom (",&
		  configurationRefEnergy(i)," eV Ref)"
		Else
		  write(999,"(I4,F12.6,A21)") i,&
		  (configurationEnergy(i)/configAtoms(i))," eV per atom (no ref)"
		End If
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
		If(configurationBM(i).gt.-2.0D20)Then
		  If(configurationRefBM(i).gt.-2.0D20)Then
		    write(999,"(I4,A20,F18.10,A9,F18.10)") i,&
		    "  Bulk Modulus/GPa: ",configurationBM(i)," Ref BM: ",&
		    configurationRefBM(i)
		  Else
		    write(999,"(I4,A20,F18.10,A13)") i,&
		    "  Bulk Modulus/GPa: ",configurationBM(i)," (No Ref BM) "
		  End If
		End If
		If(configurationRSS(i,1).gt.0.0D0)Then
		  write(999,"(I4,A14,F18.10)") i,&
		  "  Energy RSS: ",configurationRSS(i,1)
		End If
		If(configurationRSS(i,2).gt.0.0D0)Then
		  write(999,"(I4,A13,F18.10)") i,&
		  "  Force RSS: ",configurationRSS(i,2)
		End If
	    write(999,"(A1)") " "
	  End Do 
	  If(trialResidualSquareSum.gt.-2.0D20)Then
		write(999,"(A24,E20.10)") "RSS all configurations: ",trialResidualSquareSum
	  End If
	  !If(eamFunctionWobbliness(1).gt.0.0D0)Then	    
		!write(999,"(A24)") "EAM Functions Wobbliness: "
		!Do i=1,size(eamFunctionWobbliness,1)
		!  write(999,"(A4,I4,E20.10)") "Pot ",i,eamFunctionWobbliness(i)
		!End Do
	  !End If
	  write(999,"(A1)") " "
!Close output file
      close(999)
	End If  
  End Subroutine calcOutput
  
  
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
	  write(999,"(A19,I4,A7,E20.10)") "Evaluation count: ",&
	  globalCounter(1),", RSS: ",trialResidualSquareSum
!Close output file
      close(999)
	End If  
  End Subroutine calcOutputDifference
  
    
   
   
   
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

  
  Function CalcBulkModulus (configurationID) RESULT (bulkModulus) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, l, n, configurationID, keyIJ, forceKey  
	Real(kind=DoubleReal) :: strainAmount, strainIncrement, configEnergyVal	
	Real(kind=DoubleReal) :: volume, bulkModulus
	Real(kind=DoubleReal), Dimension(1:2) :: secondDerivative
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: strainArray, volEnergy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bmFit	
!Allocate arrays
	Allocate(StrainArray(1:3,1:3))
	Allocate(volEnergy(1:9,1:2))
!Store original neighbour list
    Call storeNeighboutList()
!set increment
	strainIncrement = 0.01D0
!Loop over strains
    j = 0
	Do i=-4,4
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
	 ! bulkModulus = minimumVolume * calcPolynomial (bmFit, minimumVolume, 2)
	 ! bulkModulus = bulkModulus * elementaryCharge * 1.0D30 * 1.0D-9
    End Do
	secondDerivative = CalcSecondDerivativeMinimum(volEnergy)
	bulkModulus = secondDerivative(2)*secondDerivative(1)
	bulkModulus = bulkModulus*elementaryCharge*1.0D30*1.0D-9
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
  

End Module calc