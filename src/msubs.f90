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

! Distribute from Master to Workers
  Public :: M_distChar
  Public :: M_distInt
  Public :: M_distDouble
  Public :: M_distLogical
  Public :: M_distInt1D
  Public :: M_distDouble1D
  Public :: M_distDouble2D
  Public :: M_collDouble1D
  Public :: M_collDouble1DMap
  Public :: M_collDouble2D
  Public :: M_sumDouble
  
! Distribute custom types  

! Functions

  Contains

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

  Subroutine M_synchProcesses()
! Synchronises processes by sending an integer from root to workers
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, send, receive, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
! Send out data to all processes from root
    If(processCount.gt.1)Then
      If(processID.eq.0)Then
        send = 123
        Do i=1,(processCount-1)
          processTo = i
          tag = 97220 + i
          Call MPI_SSEND(send,1,MPI_INTEGER,processTo,tag,&
          MPI_COMM_WORLD,error) ! Mpi_ssend always waits
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 97220 + processID
        Call MPI_RECV(receive,1,MPI_INTEGER,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
      End If
    End If
  End Subroutine M_synchProcesses

! ---------------------------------------------------------------------------------------------------
! Distribute char/dp/int subroutines
! ---------------------------------------------------------------------------------------------------

! ------------------------------
! M_distChar
! ------------------------------
  Subroutine M_distChar(sendString)
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(*) :: sendString
    Integer(kind=StandardInteger) :: stringLen
    Integer(kind=StandardInteger) :: i, send, receive, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        send = Len(sendString)
        Do i=1,(processCount-1)
          stringLen = send
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(send,1,MPI_INTEGER,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(receive,1,MPI_INTEGER,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        stringLen = receive
      End If
      Call M_distCharA(sendString,stringLen)
    End If
  End Subroutine M_distChar
  Subroutine M_distCharA(sendStringIn, stringLen)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: stringLen
    Character(Len=stringLen) :: sendStringIn, sendString, recvString
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        Do i=1,(processCount-1)
          sendString = sendStringIn
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(sendString,stringLen,MPI_CHARACTER,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(recvString,stringLen,MPI_CHARACTER,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendStringIn = recvString
      End If
    End If
  End Subroutine M_distCharA

! ------------------------------
! M_distDouble
! ------------------------------
  Subroutine M_distDouble(sendIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: sendIn, send, receive
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        Do i=1,(processCount-1)
          send = sendIn
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(send,1,MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(receive,1,MPI_DOUBLE_PRECISION,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendIn = receive
      End If
    End If
  End Subroutine M_distDouble

! ------------------------------
! M_distInt
! ------------------------------
  Subroutine M_distInt(sendIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: sendIn, send, receive
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        Do i=1,(processCount-1)
          send = sendIn
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(send,1,MPI_INTEGER,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(receive,1,MPI_INTEGER,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendIn = receive
      End If
    End If
  End Subroutine M_distInt
    
! ------------------------------
! M_distInt
! ------------------------------
  Subroutine M_distLogical(sendIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Logical :: sendIn, send, receive
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        Do i=1,(processCount-1)
          send = sendIn
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(send,1,MPI_LOGICAL,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(receive,1,MPI_LOGICAL,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendIn = receive
      End If
    End If
  End Subroutine M_distLogical

! ---------------------------------------------------------------------------------------------------
! Distribute array subroutines
! ---------------------------------------------------------------------------------------------------

! ------------------------------
! M_distInt1D
! ------------------------------
  Subroutine M_distInt1D(sendIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), Dimension(:) :: sendIn
    Integer(kind=StandardInteger), Dimension(1:size(sendIn)) :: send, receive
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        Do i=1,(processCount-1)
          send = sendIn
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(send,size(sendIn),MPI_INTEGER,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(receive,size(sendIn),MPI_INTEGER,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendIn = receive
      End If
    End If
  End Subroutine M_distInt1D

! ------------------------------
! M_distDouble1D
! ------------------------------
  Subroutine M_distDouble1D(sendIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(:) :: sendIn
    Real(kind=DoubleReal), Dimension(1:size(sendIn)) :: send, receive
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send out data to all processes from root
      If(processID.eq.0)Then
        Do i=1,(processCount-1)
          send = sendIn
          processTo = i
          tag = 1000 + i
          Call MPI_SEND(send,size(sendIn),MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,error)
        End Do
      End If
! Collect from root by workers
      If(processID.gt.0)Then
        processFrom = 0
        tag = 1000 + processID
        Call MPI_RECV(receive,size(sendIn),MPI_DOUBLE_PRECISION,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        sendIn = receive
      End If
    End If
  End Subroutine M_distDouble1D

! ------------------------------
! M_distDouble2D
! ------------------------------
  Subroutine M_distDouble2D(sendIn)
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal), Dimension(:,:) :: sendIn
    Real(kind=DoubleReal), Dimension(1:size(sendIn,1)) :: send, receive
    Integer(kind=StandardInteger) :: i, j, k, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer(kind=StandardInteger) :: arrayX, arrayY
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Init variables
      arrayX = size(sendIn,1)
      arrayY = size(sendIn,2)
! Send out data to all processes from root
      Do j=1,arrayY
        If(processID.eq.0)Then
          Do i=1,(processCount-1)
            Do k=1,arrayX
              send(k) = sendIn(k,j)
            End Do
            processTo = i
            tag = 1000 + i
            Call MPI_SEND(send,arrayX,MPI_DOUBLE_PRECISION,processTo,tag,&
            MPI_COMM_WORLD,error)
          End Do
        End If
! Collect from root by workers
        If(processID.gt.0)Then
          processFrom = 0
          tag = 1000 + processID
          Call MPI_RECV(receive,arrayX,MPI_DOUBLE_PRECISION,processFrom,tag,&
          MPI_COMM_WORLD,status,error)
          Do k=1,arrayX
            sendIn(k,j) = receive(k)
          End Do
        End If
      End Do
    End If
  End Subroutine M_distDouble2D

! ---------------------------------------------------------------------------------------------------
! Collect array subroutines
! ---------------------------------------------------------------------------------------------------

! ------------------------------
! M_collDouble1D
! ------------------------------

  Subroutine M_collDouble1D(sendIn, selectedProcessIn, startKeyIn, endKeyIn)
! Subroutine M_collDouble1D(sendIn, startKeyIn, endKeyIn)
! Collects 1D dp array from  workers back to root
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), optional :: selectedProcessIn, startKeyIn, endKeyIn
    Integer(kind=StandardInteger) :: selectedProcess, startKey, endKey
    Integer(kind=StandardInteger) :: arrayX
    Real(kind=DoubleReal), Dimension(:) :: sendIn
    Real(kind=DoubleReal), Dimension(1:size(sendIn)) :: send, receive
    Integer(kind=StandardInteger) :: n, i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Optional variables
    selectedProcess = -1
    startKey = 1
    endKey = size(sendIn,1)
    If(present(startKeyIn))Then
      startKey = startKeyIn
    End If
    If(present(endKeyIn))Then
      endKey = endKeyIn
    End If
    If(present(selectedProcessIn))Then
      selectedProcess = selectedProcessIn
    End If
! Init variables
    n = 0
    send = 0.0D0
    receive = 0.0D0
    arrayX = size(sendIn,1)
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! All processes back to root
      If(selectedProcess.eq.-1)Then
! Send from workers
        If(processID.gt.0)Then
          send = sendIn
          processTo = 0
          tag = 5000 + processID
          Call MPI_SEND(send,arrayX,MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,status,error)
        End If
! Collect from workers by root
        If(processID.eq.0)Then
          Do n=1,processCount-1
            processFrom = n
            tag = 5000 + n
            Call MPI_RECV(receive,arrayX,MPI_DOUBLE_PRECISION,processFrom,tag,&
            MPI_COMM_WORLD,status,error)
            Do i=startKey,endKey
              If(n.eq.mod(i-1,processCount))Then
                sendIn(i) = receive(i)
              End If
            End Do
          End Do
        End If
! Selected process back to root
      Else
! Send from selected process
        If(processID.eq.selectedProcess)Then
          send = sendIn
          processTo = 0
          tag = 5000 + processID
          Call MPI_SEND(send,arrayX,MPI_DOUBLE_PRECISION,processTo,tag,&
          MPI_COMM_WORLD,status,error)
        End If
! Collect from selected process by root
        If(processID.eq.0)Then
          processFrom = selectedProcess
          tag = 5000 + selectedProcess
          Call MPI_RECV(receive,arrayX,MPI_DOUBLE_PRECISION,processFrom,tag,&
          MPI_COMM_WORLD,status,error)
          Do i=startKey,endKey
            sendIn(i) = receive(i)
          End Do
        End If
      End If
    End If
  End Subroutine M_collDouble1D

! ------------------------------
! M_collDouble1D
! ------------------------------

  Subroutine M_collDouble1DMap(sendIn, mapIn, mapCol)
! Collects 1D dp array from  workers back to root
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: arrayX, mapCol
    Real(kind=DoubleReal), Dimension(:) :: sendIn
    Integer(kind=StandardInteger), Dimension(:,:) :: mapIn
    Real(kind=DoubleReal), Dimension(1:size(sendIn)) :: send, receive
    Integer(kind=StandardInteger) :: n, i, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Init variables
    n = 0
    send = 0.0D0
    receive = 0.0D0
    arrayX = size(sendIn,1)
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! All processes back to root
      If(processID.gt.0)Then
        send = sendIn
        processTo = 0
        tag = 1000 + processID
        Call MPI_SEND(send,arrayX,MPI_DOUBLE_PRECISION,processTo,tag,&
        MPI_COMM_WORLD,status,error)
      End If
! Collect from processes by root
      If(processID.eq.0)Then
        Do n=1,processCount-1
          processFrom = n
          tag = 1000 + n
          Call MPI_RECV(receive,arrayX,MPI_DOUBLE_PRECISION,processFrom,tag,&
          MPI_COMM_WORLD,status,error)
          Do i=1,size(receive)
            If(n.eq.mapIn(i,mapCol))Then
              sendIn(i) = receive(i)
            End If
          End Do
        End Do
      End If
    End If
  End Subroutine M_collDouble1DMap

! ------------------------------
! M_collDouble2D
! ------------------------------

  Subroutine M_collDouble2D(sendIn, selectedProcessIn, startKeyIn, endKeyIn)
! Collects 1D dp array from  workers back to root
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger), optional :: selectedProcessIn, startKeyIn, endKeyIn
    Integer(kind=StandardInteger) :: selectedProcess, startKey, endKey
    Integer(kind=StandardInteger) :: arrayX, arrayY
    Real(kind=DoubleReal), Dimension(:,:) :: sendIn
    Real(kind=DoubleReal), Dimension(1:size(sendIn)) :: send, receive
    Integer(kind=StandardInteger) :: n, i, j, processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Optional variables
    selectedProcess = -1
    startKey = 1
    endKey = size(sendIn,1)
    If(present(startKeyIn))Then
      startKey = startKeyIn
    End If
    If(present(endKeyIn))Then
      endKey = endKeyIn
    End If
    If(present(selectedProcessIn))Then
      selectedProcess = selectedProcessIn
    End If
! Init variables
    n = 0
    send = 0.0D0
    receive = 0.0D0
    arrayX = size(sendIn,1)
    arrayY = size(sendIn,2)
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! All processes back to root
      If(selectedProcess.eq.-1)Then
        Do j=1,arrayY
! Send from workers
          If(processID.gt.0)Then
! Prep send array
            Do i=1,arrayX
              send(i) = sendIn(i,j)
            End Do
            processTo = 0
            tag = 1000 + processID
            Call MPI_SEND(send,arrayX,MPI_DOUBLE_PRECISION,processTo,tag,&
            MPI_COMM_WORLD,status,error)
          End If
! Collect from workers by root
          If(processID.eq.0)Then
            Do n=1,processCount-1
              processFrom = n
              tag = 1000 + n
              Call MPI_RECV(receive,arrayX,MPI_DOUBLE_PRECISION,processFrom,tag,&
              MPI_COMM_WORLD,status,error)
              Do i=startKey,endKey
                If(n.eq.mod(i-1,processCount))Then
                  sendIn(i,j) = receive(i)
                End If
              End Do
            End Do
          End If
        End Do
! Selected process back to root
      ElseIf(selectedProcess.gt.0)Then   ! Except 0
        Do j=1,arrayY
! Send from selected process
          If(processID.eq.selectedProcess)Then
            send(i) = 0.0D0
            Do i=1,arrayX
              send(i) = sendIn(i,j)
            End Do
            processTo = 0
            tag = 11200 + 100*j + processID
! print *,"send ",processID,processTo,tag
            Call MPI_SEND(send,arrayX,MPI_DOUBLE_PRECISION,processTo,tag,&
            MPI_COMM_WORLD,status,error)
          End If
! Collect from selected process by root
          If(processID.eq.0)Then
            processFrom = selectedProcess
            tag = 11200 + 100*j + selectedProcess
! print *,"recv ",processID,processFrom,tag
            Call MPI_RECV(receive,arrayX,MPI_DOUBLE_PRECISION,processFrom,tag,&
            MPI_COMM_WORLD,status,error)
            Do i=startKey,endKey
              sendIn(i,j) = receive(i)
            End Do
          End If
        End Do
      End If
    End If
  End Subroutine M_collDouble2D

! ---------------------------------------------------------------------------------------------------
! Sum variable/array subroutines
! ---------------------------------------------------------------------------------------------------

! ------------------------------
! M_collDouble1D
! ------------------------------

  Subroutine M_sumDouble(procA, procB, doubleIn)
! Adds doubleIn on procB to doubleIn on procA
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: procA, procB
    Real(kind=DoubleReal) :: doubleIn, send, receive
    Integer(kind=StandardInteger) :: processTo, processFrom, tag
    Integer(kind=StandardInteger) :: processID,processCount,error
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Init variables
    send = 0.0D0
    receive = 0.0D0
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
    If(processCount.gt.1)Then
! Send from proc B
      If(processID.eq.procB)Then
        processTo = procA
        tag = 3211+procB
        send = doubleIn
        Call MPI_SEND(send,1,MPI_DOUBLE_PRECISION,processTo,tag,&
        MPI_COMM_WORLD,status,error)
      End If
! Recv by proc A and add to existing value
      If(processID.eq.procA)Then
        processFrom = procB
        tag = 3211+procB
        Call MPI_RECV(receive,4,MPI_DOUBLE_PRECISION,processFrom,tag,&
        MPI_COMM_WORLD,status,error)
        doubleIn = doubleIn + receive
      End If
    End If
  End Subroutine M_sumDouble

! UNDER DEVELOPMENT

End Module msubs
