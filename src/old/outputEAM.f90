Module outputeam

! Setup Modules
  Use kinds
  Use constants
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

!Privacy of functions/subroutines/variables
  Private
!Variables  
!Subroutines  
  Public :: outputPreparedEAMFile
!Functions
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

  Subroutine outputPreparedEAMFile(eamKeyArray, eamDataArray, fileName, numberOfPointsIn, processIn)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints, processFlag
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	Character(*) :: fileName
!optional variables	
	Integer(kind=StandardInteger), optional :: numberOfPointsIn
	Integer(kind=StandardInteger), optional :: processIn
	  If(Present(numberOfPointsIn))Then
	    numberOfPoints = numberOfPointsIn
	  Else
        numberOfPoints = 0	  
	  End If
	  If(Present(processIn))Then
	    processFlag = processIn
	  Else
	    processFlag = -1
	  End If	
!Store reduced potential to file
	If(mpiProcessID.eq.0)Then
	  open(unit=24,file=trim(fileName))
	  If(numberOfPoints.lt.10)Then
	    Do potKey=1,size(eamKeyArray,1)
!Get number of points to reduce to
	      If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	        write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyArray(potKey,1))," ",&
		    elements(eamKeyArray(potKey,2))," "
		  Else  
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		      eamTypeText = "DENS "
		    End If
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		      eamTypeText = "SDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		      eamTypeText = "EMBE "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		      eamTypeText = "SEMB "
		    End If
		    If(eamKeyArray(potKey,3).eq.4)Then
		      eamTypeText = "DDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.5)Then
		      eamTypeText = "DEMB "
		    End If
		    write(24,"(A5,A2)") eamTypeText,elements(eamKeyArray(potKey,1))
	      End If	
!pot positions
	      potStart = eamKeyArray(potKey,4)
          potLength = eamKeyArray(potKey,5)
	      potEnd = potStart + potLength - 1
!loop over data points
          k = 0
	      Do i=potStart,potEnd
		    k = k + 1
	        x = eamDataArray(i,1)
		    y = eamDataArray(i,2)
		    dy = eamDataArray(i,3)
!write to file
            write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		    y,"  ",dy,"    ",potKey,k,i
	      End Do	
        End Do
	  ElseIf(numberOfPoints.ge.10)Then
!Interpolate between points
!Loop through potential functions
        dataPointCounter = 1
	    Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	      If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	        write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyArray(potKey,1))," ",&
		    elements(eamKeyArray(potKey,2))," "
		  Else  
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		      eamTypeText = "DENS "
		    End If
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		      eamTypeText = "SDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		      eamTypeText = "EMBE "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		      eamTypeText = "SEMB "
		    End If
		    If(eamKeyArray(potKey,3).eq.4)Then
		      eamTypeText = "DDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.5)Then
		      eamTypeText = "DEMB "
		    End If
		    write(24,"(A5,A2)") eamTypeText,elements(eamKeyArray(potKey,1))
	      End If	
	      potStart = eamKeyArray(potKey,4)
          potLength = eamKeyArray(potKey,5)
	      potEnd = potStart + potLength - 1
!loop through reduced data points
!Interpolate 4th order polynomial (5 data points)
	      Do i=1,numberOfPoints
	        x = eamDataArray(potStart,1)+&
		    1.0D0*(i-1)*((eamDataArray(potEnd,1)-&
			eamDataArray(potStart,1))/(numberOfPoints-1))
		    !yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength,"N")
		    yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength)
		    y = yArray(1)
		    dy = yArray(2)
!write to file
            write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		    y,"  ",dy,"    ",potKey,i,dataPointCounter	
!increment counter
		    dataPointCounter = dataPointCounter + 1		
	      End Do	
        End Do
	  End If
!Close file
	  close(24)
	End If
!Wait for all processes to catch up	 
    If(processFlag.lt.0)Then
	  Call synchMpiProcesses()
	End If  
  End Subroutine outputPreparedEAMFile 

  

  
  

End Module outputeam