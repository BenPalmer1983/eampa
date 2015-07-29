Module typesM
! --------------------------------------------------------------!
! Data typesMPI
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Define data typesMPI
! ----------------------------------------
! Updated: 24th June 2015
! ----------------------------------------
! Setup Modules  
  Use kinds
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Declare global variables
! Privacy of functions/subroutines/variables
  Private
! Public subroutines  
  Public :: SetUpMPITypes

  Contains
  
  Subroutine SetUpMPITypes()  
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: error
    Integer(kind=StandardInteger) :: iExtent, rExtent, dpExtent
    Integer(kind=StandardInteger) :: elementTypes(0:1), blockCounts(0:1), blockOffset(0:1)
    Integer(kind=StandardInteger) :: rssConfigID
    Integer(kind=StandardInteger) :: mpiProcessCount, mpiProcessID
    Integer(kind=StandardInteger) :: i, processTo, processFrom, tag
    Integer, Dimension(MPI_STATUS_SIZE) :: status
! Define my rssConfig type in Fortran
    Type :: rssConfig
      sequence
      Real(kind=DoubleReal) :: total=0.0D0
      Real(kind=DoubleReal) :: energy=0.0D0
      Real(kind=DoubleReal) :: force=0.0D0
      Real(kind=DoubleReal) :: stress=0.0D0
      Integer(kind=StandardInteger) :: n
    End Type rssConfig
! Declare rssTest as rssConfig type
    Type (rssConfig) :: rssTest  
! Extent of data types
    Call MPI_TYPE_EXTENT(MPI_INTEGER, iExtent, error)
    Call MPI_TYPE_EXTENT(MPI_REAL, rExtent, error)
    Call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dpExtent, error)    
! Block 1
    blockCounts(0) = 4
    blockOffset(0) = 0
    elementTypes(0) = MPI_DOUBLE_PRECISION
! Block 2
    blockCounts(1) = 1 
    blockOffset(1) = 4 * dpExtent
    elementTypes(1) = MPI_INTEGER
! Make mpi structure
    Call MPI_TYPE_STRUCT(2, blockCounts, blockOffset, elementTypes, rssConfigID, error)
    Call MPI_TYPE_COMMIT(rssConfigID, error)
! Get process ID and total mpi process count
    Call MPI_Comm_size(MPI_COMM_WORLD ,mpiProcessCount,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiProcessID,error)
! Set values on the root process    
    If(mpiProcessID.eq.0)Then
      rssTest%total = 8.1D0
    End If
! wait and print out  
    !print *,mpiProcessID,"px",rssTest%total
! Send from root to workers using rssConfig type
    If(mpiProcessID.eq.0)Then
      Do i=1,(mpiProcessCount-1)
        processTo = i
        tag = 114 + i
        Call MPI_SEND(rssTest,1,rssConfigID,processTo,tag,MPI_COMM_WORLD,error)
      End Do
    End If
! Recieve by worker processes from root process
    If(mpiProcessID.gt.0)Then
      processFrom = 0
      tag = 114 + mpiProcessID
      Call MPI_RECV(rssTest,1,rssConfigID,processFrom,tag,MPI_COMM_WORLD,status,error)
    End If  
! wait and printout
    !print *,mpiProcessID,"px",rssTest%total  
  End Subroutine SetUpMPITypes

End Module typesM