Module msubs

! Setup Modules
  Use kinds

! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Declare global variables 
! Privacy of functions/subroutines/variables
  Private
! Variables


! Subroutines
  Public :: M_synchProcesses
! Send data from Master to Workers - fixed arrays
  Public :: M_sendData1DDP
  Public :: M_sendData1DInt
  Public :: M_sendData2DInt
  Public :: M_sendData2DDP

! Distribute from Master to Workers  
  Public :: M_dist_char
  
! Sum data from master and workers
  Public :: M_sumData2DDP
  Public :: M_sumData1DDP
  
  
!Functions
  
  
  Contains  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!   



  Subroutine M_synchProcesses()
!force declaration of all variables
    Implicit None   
!Internal subroutine variables    
    Integer(kind=StandardInteger) :: i, send, receive, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
!call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
! Send out data to all processes from root
    If(processID.eq.0) Then
      send = 123
      Do i=1,(processCount-1)
        processTo = i
        tag = 1000 + i  
        Call MPI_SEND(send,1,MPI_INTEGER,processTo,tag,&
        MPI_COMM_WORLD,error)
      End Do
    End If
! Collect from root by workers
    If(processID.gt.0) Then      
      processFrom = 0
      tag = 1000 + processID
      Call MPI_RECV(receive,1,MPI_INTEGER,processFrom,tag,&
      MPI_COMM_WORLD,status,error)  
    End If     
    
  End Subroutine M_synchProcesses 



  
  
  
!--------------------------------------------------------------------------------------------------- 
! Distribute array subroutines
!--------------------------------------------------------------------------------------------------- 
 
!------------------------------------------------------------------------!
! ***MPI_sendData1DInt*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine M_sendData1DInt(dataArray,arraySize) 
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i
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
  End Subroutine M_sendData1DInt 
!------------------------------------------------------------------------!
! ***MPI_sendData1DDP*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine M_sendData1DDP(dataArray,arraySize) 
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i
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
    MPI_integer,processTo,tag,MPI_comm_world,error)
    End Do  
  Else
!RECV by worker processes
      processFrom = 0
    tag = 2000 + processID
      Call MPI_recv(recvArray,arraySize,&
    MPI_integer,processFrom,tag,MPI_comm_world,status,error)
    dataArray = recvArray
  End If    
  End Subroutine M_sendData1DDP
!------------------------------------------------------------------------!
! ***MPI_sendData2DInt*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine M_sendData2DInt(dataArray,arraySizeH,arraySizeW) 
!force declaration of all variables
    Implicit None  
!declare private variables
    Integer(kind=StandardInteger) :: i,n
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
      MPI_integer,processTo,tag,MPI_comm_world,error)
      End Do  
    Else
!RECV by worker processes
        processFrom = 0
      tag = 2000 + processID
        Call MPI_recv(recvArray,arraySizeH,&
      MPI_integer,processFrom,tag,MPI_comm_world,status,error)
      dataArrayTemp = recvArray
    End If    
!-----------------------------------------
! Merge 1D array with 2D array
!-----------------------------------------      
      Do i=1,arraySizeH
      dataArray(i,n) = dataArrayTemp(i)
    End Do
  End Do
  End Subroutine M_sendData2DInt 
!------------------------------------------------------------------------!
! ***MPI_sendData2DDP*** 
! Takes an array and process map, merges array values and distributes
!------------------------------------------------------------------------!  
  Subroutine M_sendData2DDP(dataArray,arraySizeH,arraySizeW) 
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i,n
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
  End Subroutine M_sendData2DDP 
  
  
  
  
  
  Subroutine M_dist_char(sendString)
    Implicit None   ! Force declaration of all variables
! Private variables    
    Character(*) :: sendString
    Integer(kind=StandardInteger) :: i, sendInt, strLen
    Real(kind=DoubleReal) :: sendD, recvD
! Mpi variables
    Integer(kind=StandardInteger) :: processID,processCount
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
! Init variables
    strLen = 0
! Set mpi variable values
    Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)      
!-----------------------------------------
! Send/Recv 
!-----------------------------------------
    print *,".....",processID,strLen
    If(processID.eq.0)Then
!SEND from master process
      sendInt = Len(sendString)
      strLen = sendInt
      sendD = 1.0D0
      recvD = sendD
      Do i=1,(processCount-1)
        processTo = i
        tag = 1000 + i
        print *,"send to ",i,sendD,tag
        Call MPI_send(sendD,1,&
        MPI_DOUBLE_PRECISION,processTo,tag,MPI_comm_world,error)
      End Do  
    Else  
!RECV by master process
      processFrom = 0
      tag = 1000 + processID
      !Call MPI_recv(recvInt,1,&
      !MPI_INTEGER,processFrom,tag,MPI_comm_world,status,error)
      Call MPI_recv(recvD,1,&
      MPI_DOUBLE_PRECISION,processFrom,tag,MPI_comm_world,status,error)
      print *,"collect on ",processID, recvD, tag
      !strLen = recvInt
    End If
    print *,"..",processID,recvD
  
  
  End Subroutine M_dist_char 
  !Subroutine MPI_dist_char_send(sendString)
  !  Implicit None   ! Force declaration of all variables
! Private variables    
  !  Character(*) :: sendString
  !MPI_CHAR
  
  
  !End Subroutine MPI_dist_char_send 
  
  
!------------------------------------------------------------------------!
! ***MPI_sumData2DDP*** 
! Takes an array and sums values across mpi processes
!------------------------------------------------------------------------!  
  Subroutine M_sumData2DDP(dataArray,arraySizeH,arraySizeW) 
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i,j,n
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
! Sum arrays from workers with master
!-----------------------------------------       
    If(processID.eq.0)Then
!RECV by master process
      Do i=1,(processCount-1)
          processFrom = i
        tag = 1000 + i
          Call MPI_recv(recvArray,arraySizeH,&
        MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
      Do j=1,arraySizeH
          dataArrayTemp(j) = dataArrayTemp(j) + recvArray(j)
      End Do
    End Do   
    Else
!SEND by worker process   
      sendArray = dataArrayTemp
      processTo = 0
        tag = 1000 + processID
    Call MPI_send(sendArray,arraySizeH,&
    MPI_double_precision,processTo,tag,MPI_comm_world,error)
    End If 
!-----------------------------------------
! Merge 1D array with 2D array
!-----------------------------------------      
      Do i=1,arraySizeH
      dataArray(i,n) = dataArrayTemp(i)
    End Do  
  End Do
!-----------------------------------------
! Send out merged/summed array
!-----------------------------------------  
    Call M_sendData2DDP(dataArray,arraySizeH,arraySizeW) 
  End Subroutine M_sumData2DDP 
  
!------------------------------------------------------------------------!
! ***MPI_sumData1DDP*** 
! Takes an array and sums values across mpi processes
!------------------------------------------------------------------------!  
  Subroutine M_sumData1DDP(dataArray,arraySize) 
!force declaration of all variables
  Implicit None  
!declare private variables
  Integer(kind=StandardInteger) :: i,j
  Integer(kind=StandardInteger) :: arraySize
!mpi variables
    Integer(kind=StandardInteger) :: processID,processCount
    Integer(kind=StandardInteger) :: status,error,tag
    Integer(kind=StandardInteger) :: processTo,processFrom
  Real(kind=DoubleReal), Dimension(1:arraySize) :: dataArray
  Real(kind=DoubleReal), Dimension(1:arraySize) :: sendArray, recvArray  
!-----------------------------------------
! Set mpi variable values
!-----------------------------------------
    Call MPI_comm_rank(MPI_COMM_WORLD,processID,error)
      Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)  
!-----------------------------------------  
! Sum arrays from workers with master
!-----------------------------------------       
    If(processID.eq.0)Then
!RECV by master process
      Do i=1,(processCount-1)
          processFrom = i
        tag = 1000 + i
          Call MPI_recv(recvArray,arraySize,&
        MPI_double_precision,processFrom,tag,MPI_comm_world,status,error)
      Do j=1,arraySize
          dataArray(j) = dataArray(j) + recvArray(j)
      End Do
    End Do   
    Else
!SEND by worker process   
      sendArray = dataArray
      processTo = 0
        tag = 1000 + processID
    Call MPI_send(sendArray,arraySize,&
    MPI_double_precision,processTo,tag,MPI_comm_world,error)
    End If   
!-----------------------------------------
! Send out merged/summed array
!-----------------------------------------  
      Call M_sendData1DDP(dataArray,arraySize)   
  End Subroutine M_sumData1DDP 
  
  

End Module msubs