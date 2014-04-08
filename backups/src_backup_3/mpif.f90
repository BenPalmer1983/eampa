Module mpif

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
  
!Variables

!Subroutines

!Functions

  
!Privacy of functions/subroutines/variables
  Private
  Public :: MPI_sendout
  Public :: MPI_sendout2D
  
  
  
  Contains

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
  
  
  
  
  
  
  
  
  

End Module mpif