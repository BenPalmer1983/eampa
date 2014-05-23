Module mpif

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!declare global variables 
  Real(kind=DoubleReal) :: mpiProcessTimeStart
!Privacy of functions/subroutines/variables
  Private
!variables
  Public :: mpiProcessTimeStart
!Subroutines
  Public :: MPI_distributeArray1D
  Public :: MPI_sendout
  Public :: MPI_sendout2D
  Public :: MPI_timer  
  Public :: MPI_sendData1DDP
  Public :: MPI_sendData2DInt
  Public :: MPI_sendData2DDP
!Functions
  !Public :: 
  
  
  Contains  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!   



  Subroutine MPI_distributeArray1D(distributedArray,processMap)
! Takes an array, with values on different processes, combines and distributes
!
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: distributedArray	
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: processMap
!mpi variables
    Integer(kind=StandardInteger) :: processID, processCount
    Integer(kind=StandardInteger) :: status,error,tag,processTemp
    Integer(kind=StandardInteger) :: processTo,processFrom
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferArray		
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
	!Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
    !Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)


  End Subroutine MPI_distributeArray1D




 

!send array - just doubles at the moment
  Subroutine MPI_sendout(buffer,dataType,bufferDimensionOption)
 !force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,processID,processCount,tag
	Integer(kind=StandardInteger) :: status,error,bufferSize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: buffer
	Character(*), INTENT(IN) :: dataType
	Integer(kind=StandardInteger) :: bufferDimension
	Integer(kind=StandardInteger), optional :: bufferDimensionOption
!get process id and count
	Call MPI_comm_rank( MPI_COMM_WORLD,processID,error )
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!optional variables	
	If(Present(bufferDimensionOption))Then
	  bufferDimension = bufferDimensionOption
	Else
	  bufferDimension = 1
	End If	
!declare array
    	
!-----------
! 1D Array 
!-----------
    If(bufferDimension.eq.1)Then
!set variables
      bufferSize = size(buffer,1)
!Send out data to all processes
      If(processID.eq.0) Then
        Do i=1,(processCount-1)
          tag = 1000 + i	
		  If(dataType(1:7).eq."integer")Then
            Call MPI_send(buffer,bufferSize,MPI_integer,i,tag,&
		    MPI_comm_world,error)
		  End If
		  If(dataType(1:6).eq."double")Then
            Call MPI_send(buffer,bufferSize,MPI_double_precision,i,tag,&
		    MPI_comm_world,error)
		  End If
        End Do
      Else
        tag = 1000 + processID
	    If(dataType(1:7).eq."integer")Then
          Call MPI_recv(buffer,bufferSize,MPI_integer,0,tag,&
          MPI_comm_world,status,error)  
	    End If
	    If(dataType(1:6).eq."double")Then
	      Call MPI_recv(buffer,bufferSize,MPI_double_precision,0,tag,&
          MPI_comm_world,status,error) 
	    End If 
      End If
	End If
  End Subroutine MPI_sendout
  
  
  
  Subroutine MPI_sendcollect(buffer)
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,processID,processCount,tag
	Integer(kind=StandardInteger) :: status,error,bufferSize
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: buffer
!get process id and count
	Call MPI_comm_rank( MPI_COMM_WORLD,processID,error )
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)    
  
  
  
  
  End Subroutine MPI_sendcollect
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
!Send 2D array - not working yet
  !send array
  Subroutine MPI_sendout2D(buffer,dataType,bufferDimensionOption)
 !force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,j,processID,processCount,tag
	Integer(kind=StandardInteger) :: status,error,bufferSize,bufferWidth
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: buffer
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: bufferTemp
	Character(*), INTENT(IN) :: dataType
	Integer(kind=StandardInteger) :: bufferDimension
	Integer(kind=StandardInteger), optional :: bufferDimensionOption
!get process id and count
	Call MPI_comm_rank( MPI_COMM_WORLD,processID,error )
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!optional variables	
	If(Present(bufferDimensionOption))Then
	  bufferDimension = bufferDimensionOption
	Else
	  bufferDimension = 2
	End If	
!declare array
    	
!-----------
! 2D Array 
!-----------
    If(bufferDimension.eq.2)Then
!set variables
      bufferSize = size(buffer,1)
	  bufferWidth = size(buffer,2)
!loop through columns in array	  
	  Do j=1,bufferWidth
!make 1D array	    
		If(Allocated(bufferTemp))Then
		Else
		  Allocate(bufferTemp(1:bufferSize))
		End If
!Send out data to all processes
        If(processID.eq.0) Then
		  Do i=1,bufferSize
		    bufferTemp(i) = buffer(i,j)
		  End Do
          Do i=1,(processCount-1)
            tag = 1000 + i + 10 * j	
		    If(dataType(1:7).eq."integer")Then
              Call MPI_send(buffer,bufferSize,MPI_integer,i,tag,&
		      MPI_comm_world,error)
		    End If
		    If(dataType(1:6).eq."double")Then
              Call MPI_send(bufferTemp,bufferSize,MPI_double_precision,i,tag,&
		      MPI_comm_world,error)
		    End If
          End Do
        Else
          tag = 1000 + processID + 10 * j
	      If(dataType(1:7).eq."integer")Then
            Call MPI_recv(bufferTemp,bufferSize,MPI_integer,0,tag,&
            MPI_comm_world,status,error)  
	      End If
	      If(dataType(1:6).eq."double")Then
	        Call MPI_recv(bufferTemp,bufferSize,MPI_double_precision,0,tag,&
            MPI_comm_world,status,error) 
	      End If 
		  Do i=1,bufferSize
		    buffer(i,j) = bufferTemp(i)
		  End Do
        End If
	  End Do
	End If
	
	
  End Subroutine MPI_sendout2D
  
  
  
  
  Subroutine MPI_timer(mpiProcessTime,mpiProcessTimeReal,mpiProcessTimeFlag)
!mpiProcessTimeFlag = 0 - sets start time on all processes
!mpiProcessTimeFlag = 1 - returns total processing time on all processes
!mpiProcessTime total time added over all processes
!mpiProcessTimeReal maximum time taken by 1 of the processes, real time taken
!force declaration of all variables
	Implicit None
!declare variables  
    Integer(kind=StandardInteger) :: i,j,k
	Real(kind=DoubleReal) :: cpuTime
    Real(kind=DoubleReal) :: mpiProcessTime, mpiProcessTimeReal
    Real(kind=DoubleReal) :: mpiProcessTimeDifference
	Real(kind=DoubleReal) :: bufferVal
	Integer(kind=StandardInteger) :: mpiProcessTimeFlag, loopSize
!MPI Variables
    Integer(kind=StandardInteger) :: processTo,processFrom,processID,processCount,tag
	Integer(kind=StandardInteger) :: status,error,bufferSize,processIDTemp   
!get process id and count
	Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!Start 
	loopSize = processCount - 1
    If(mpiProcessTimeFlag.eq.1)Then
!Calculate process time
      mpiProcessTime = 0.0D0
      mpiProcessTimeReal = 0.0D0
	  Call cpu_time(cpuTime) 
	  mpiProcessTimeDifference = 1.0D0*(cpuTime - mpiProcessTimeStart)	!Time elapsed on process  
	  If(processID.gt.0) Then
!send buffers from all worker processes
        processTo = 0
        Do i=1,loopSize
          If(i.eq.processID)Then
            tag = 7231 + processID	
            Call MPI_send(mpiProcessTimeDifference,1,&
		    MPI_double_precision,processTo,tag,MPI_comm_world,error)
	      End If
        End Do
	  ElseIf(processID.eq.0)Then
        mpiProcessTime = mpiProcessTimeDifference	
        mpiProcessTimeReal	= mpiProcessTimeDifference	
		Do j=1,loopSize
		  processFrom = j
          tag = 7231 + j
          Call MPI_recv(bufferVal,1,&
		  MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
		  mpiProcessTime = mpiProcessTime + bufferVal
		  If(bufferVal.gt.mpiProcessTimeReal)Then
		    mpiProcessTimeReal = bufferVal
		  End If
	    End Do 		
      End If
!Send out to all workers
      If(processID.eq.0)Then
	    Do j=1,loopSize		  
		  tag = 8231 + j
		  processTo = j	
          Call MPI_send(mpiProcessTime,1,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
		  tag = 8431 + j
		  processTo = j	
          Call MPI_send(mpiProcessTimeReal,1,&
		  MPI_double_precision,processTo,tag,MPI_comm_world,error)
	    End Do 
	  Elseif(processID.gt.0)Then
	    processIDTemp = processID
        tag = 8231 + processIDTemp
	    processFrom = 0
        Call MPI_recv(bufferVal,1,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
		mpiProcessTime = bufferVal
        tag = 8431 + processIDTemp
	    processFrom = 0
        Call MPI_recv(bufferVal,1,&
		MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
		mpiProcessTimeReal = bufferVal
	  End If
	Else
!Set start time for each process
	  Call cpu_time(cpuTime)
	  mpiProcessTimeStart = cpuTime
	End If
  End Subroutine MPI_timer
  
  
  
  
  
!--------------------------------------------------------------------------------------------------- 
! Distribute array subroutines
!--------------------------------------------------------------------------------------------------- 
 
!------------------------------------------------------------------------!
! ***MPI_sendData1DInt*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine MPI_sendData1DInt(dataArray,arraySize) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k
	Integer(kind=StandardInteger) :: arraySize
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Integer(kind=StandardInteger), Dimension(1:arraySize) :: dataArray, sendArray, recvArray	
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
  End Subroutine MPI_sendData1DInt 
!------------------------------------------------------------------------!
! ***MPI_sendData1DDP*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine MPI_sendData1DDP(dataArray,arraySize) 
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
  End Subroutine MPI_sendData1DDP
!------------------------------------------------------------------------!
! ***MPI_sendData2DInt*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine MPI_sendData2DInt(dataArray,arraySizeH,arraySizeW) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n
	Integer(kind=StandardInteger) :: arraySizeH,arraySizeW
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
	Integer(kind=StandardInteger), Dimension(1:arraySizeH,1:arraySizeW) :: dataArray
	Integer(kind=StandardInteger), Dimension(1:arraySizeH) :: dataArrayTemp, sendArray, recvArray
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
  End Subroutine MPI_sendData2DInt 
!------------------------------------------------------------------------!
! ***MPI_sendData2DDP*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine MPI_sendData2DDP(dataArray,arraySizeH,arraySizeW) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n
	Integer(kind=StandardInteger) :: arraySizeH,arraySizeW
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount
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
  End Subroutine MPI_sendData2DDP 

  
  
  
  

End Module mpif