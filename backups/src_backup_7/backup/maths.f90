Module maths

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
  Private  
  !functions
  Public :: solvePolynomial
  !subroutines
  Public :: polyFit
  Public :: splineFit
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  

  function solvePolynomial (coefficients, lower, upper) RESULT (output)
    	
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: upper, lower, output
	Real(kind=DoubleReal) :: x,y,dydx
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget, factor, difference
	Integer(kind=StandardInteger) i,j,k
	
	!Set values
	convergenceTarget = 0
	convergence = 1000
	convergenceThreshold = 0.001
	
	!set start value for x
	factor = 0.5
	difference = upper - lower
	x = lower + factor * difference	
	do while(convergence.gt.convergenceThreshold)
	  difference = factor * difference
	  y = 0
	  do i=0,size(coefficients)-1
	    y = y + x**(i) * coefficients(i)			
	  enddo
	  dydx = 0
	  do i=1,size(coefficients)-1
	    dydx = dydx + i * x**(i-1) * coefficients(i)			
	  enddo		  
	  convergence = abs(convergenceTarget - y)
	  if(convergence.gt.convergenceThreshold)then
	    if((dydx.lt.0.and.y.ge.0).or.(dydx.ge.0.and.y.lt.0))then
	      x = x + difference	
	    else
	      x = x - difference
	    endif
	  endif
	enddo
	
	output = x

  end function solvePolynomial 


!Subroutines

  Subroutine splineFit(inputPoints)
  
    !In variables
	!double precision, intent(IN), Dimension( : , : ), Allocatable :: inputPoints
	double precision, Dimension( : , : ), Allocatable :: inputPoints
	
	Allocate(inputPoints(1:17,1:2))
	
	inputPoints(1,1) = 1.5000000000000000e+00 
	inputPoints(1,2) = 4.43594
	inputPoints(2,1) = 1.7500000000000000e+00 
	inputPoints(2,2) = 3.65848
	inputPoints(3,1) = 2.0000000000000000e+00 
	inputPoints(3,2) = 1.83126
	inputPoints(4,1) = 2.2500000000000000e+00 
	inputPoints(4,2) = 0.60767
	inputPoints(5,1) = 2.5000000000000000e+00 
	inputPoints(5,2) = 0.09596
	inputPoints(6,1) = 2.7500000000000000e+00 
	inputPoints(6,2) = -0.07819
	inputPoints(7,1) = 3.0000000000000000e+00 
	inputPoints(7,2) = -0.10924
	inputPoints(8,1) = 3.2500000000000000e+00 
	inputPoints(8,2) = -0.10353
	inputPoints(9,1) = 3.3799999999999999e+00 
	inputPoints(9,2) = -0.09831
	inputPoints(10,1) = 3.5000000000000000e+00 
	inputPoints(10,2) = -0.09320
	inputPoints(11,1) = 3.7500000000000000e+00 
	inputPoints(11,2) = -0.08425
	inputPoints(12,1) = 4.0000000000000000e+00 
	inputPoints(12,2) = -0.07807
	inputPoints(13,1) = 4.5000000000000000e+00 
	inputPoints(13,2) = -0.04693
	inputPoints(14,1) = 5.0000000000000000e+00 
	inputPoints(14,2) = -0.00878
	inputPoints(15,1) = 5.5000000000000000e+00 
	inputPoints(15,2) = 0.00515
	inputPoints(16,1) = 6.0000000000000000e+00 
	inputPoints(16,2) = 0.00095
	inputPoints(17,1) = 6.5000000000000000e+00 
	inputPoints(17,2) = 0.000000
	
  
!End of spline subroutine
  End Subroutine splineFit
  
  
    
    
  

  Subroutine polyFit(inputPoints,order,polyCoefficients)
 
    
 
!force declaration of all variables
	Implicit None
	
    Integer, Parameter :: sp = Selected_Real_Kind(6,37)    ! single real
    Integer, Parameter :: dp = Selected_Real_Kind(15,307)  ! double real
	Integer, Parameter :: qp = Selected_Real_Kind(15,307) ! temporary
    Integer, Parameter :: wp = dp                          ! working real
    Integer, Parameter :: ip = Selected_Int_Kind(12)       ! long integer
	
!declare variables
    !in-out vars
	real, Dimension( : , : ), Allocatable :: inputPoints
	Real( Kind = qp ), Dimension( : ), Allocatable :: polyCoefficients
	integer :: order
	
	!internal
	Real( Kind = qp) :: matrixentry, matrixtempa, matrixtempb, idfactor, tempreal, trifactor
    Real( Kind = qp ), dimension(10,20) :: xmatrix
    Real( Kind = qp ), dimension(10,10) :: xinvmatrix
	!Real( Kind = qp ), dimension( : ), Allocatable :: xmatrix
    !Real( Kind = qp ), dimension( : ), Allocatable :: xinvmatrix
    Real( Kind = qp ), dimension(10,1) :: ymatrix, coefficients
    Real( Kind = qp ), dimension(10,1) :: rowvalue
    Real( Kind = qp ) :: randnum, lms, x, y, ycalc
    Integer :: filerow, maxrows, totalrows, ios, row, column, exponentVal, matrixsize
    Integer :: rowi, rowj, columnk, i, j
    Integer :: matrixarows, matrixacolumns, matrixbrows, matrixbcolumns, strstart, strend
    Logical :: sorting
    CHARACTER(len=20) :: temp
    CHARACTER(len=255) :: output

    !Start matrix operations
    matrixsize = order + 1
	
	!Allocate arrays
	!Allocate(xmatrix(1:matrixsize,1:2*matrixsize))
	!Allocate(xinvmatrix(1:matrixsize,1:matrixsize))
	
	!do rowi=1,size(inputPoints)/2
    !  print *,inputPoints(rowi,1),inputPoints(rowi,2)
    !end do	

    !make xmatrix
	!print *,"X Matrix"
    do row=1,matrixsize
      do column=1,matrixsize
        exponentVal = (row - 1) + (column - 1)
        matrixentry = 1.0D0
        do rowi=1,size(inputPoints)/2
          matrixentry = 1.0D0 * (matrixentry + inputPoints(rowi,1)**(1.0D0*exponentVal))
        end do
        xmatrix(row,column) = 1.0D0 * matrixentry
      end do
    end do
	
	!ymatrix
	!print *,"Y Matrix"
    do row=1,matrixsize
      matrixentry = 1.0D0
      do rowi=1,size(inputPoints)/2
        matrixentry = 1.0D0 * &
		(matrixentry + (1.0D0 * inputPoints(rowi,1)**(1.0D0 * (row-1)))*&
		(1.0D0 * inputPoints(rowi,2)))
      end do
	  !print *,row,matrixentry
      ymatrix(row,1) = matrixentry
    end do
  
    !optimise matrix order
	Do i=1,size(rowvalue)
	  rowvalue(i,1) = 1.0D0
	Enddo
    DO row=1,matrixsize
      DO column=1,matrixsize
        IF(xmatrix(row,column).EQ.0)THEN
          rowvalue(row,1) = 1.0D0 * (rowvalue(row,1) &
		  + 1.0D0*(10**(1.0D0*(matrixsize-row))))
		ELSE
		  rowvalue(row,1) = 1.0D0
        END IF
      END DO
    END DO
    sorting = .true.
    DO WHILE(sorting)
      sorting = .false.
      DO row=1,(matrixsize-1)
        IF(rowvalue(row,1).GT.rowvalue(row+1,1))THEN
          sorting = .true.
          matrixtempa = 1.0D0 * rowvalue(row,1)
          matrixtempb = 1.0D0 * rowvalue(row+1,1)
          rowvalue(row,1) = 1.0D0 * matrixtempb
          rowvalue(row+1,1) = 1.0D0 * matrixtempa
          DO column=1,matrixsize
            matrixtempa = 1.0D0 * xmatrix(row,column)
            matrixtempb = 1.0D0 * xmatrix(row+1,column)
            xmatrix(row,column) = 1.0D0 * matrixtempb
            xmatrix(row+1,column) = 1.0D0 * matrixtempa
          END DO
        END IF
      END DO
    END DO
		
	!append identity to xmatrix
    do row=1,matrixsize
      do column=matrixsize+1,2 * matrixsize
        IF((column-matrixsize).EQ.row)THEN
          xmatrix(row,column) = 1.0D0
        ELSE
          xmatrix(row,column) = 0.0D0
        END IF
      END DO
    END DO	
	
	!print *,"X Matrix"
    !do row=1,matrixsize
    !  do column=1,2 *matrixsize
	!	print *,row,column,xmatrix(row,column)
    !  end do
    !end do
  
    !make lower triangle of 0s
    DO rowi=1,matrixsize-1
      DO rowj=rowi+1,matrixsize
	    trifactor = 1.0D0 * (xmatrix(rowj,rowi)/xmatrix(rowi,rowi))
        DO columnk=1,2*matrixsize
          IF(columnk.LE.rowi)THEN
            xmatrix(rowj,columnk) = 1.0D0
          ELSE
			xmatrix(rowj,columnk) = 1.0D0 * xmatrix(rowj,columnk) - &
			1.0D0 * (trifactor) * xmatrix(rowi,columnk)
          END IF
        END DO
      END DO
    END DO
	
	
	!print *,"X Matrix"
    !do row=1,matrixsize
    !  do column=1,2 *matrixsize
	!	print *,row,column,xmatrix(row,column)
    !  end do
    !end do
	
	!make upper triangle of zeros
    DO rowi=matrixsize,2,-1
      DO rowj=rowi-1,1,-1
	    trifactor = 1.0D0 * (xmatrix(rowj,rowi)/xmatrix(rowi,rowi))
        DO columnk=1,2*matrixsize
          xmatrix(rowj,columnk) = 1.0D0 * xmatrix(rowj,columnk) - &
		  1.0D0 * (trifactor) * xmatrix(rowi,columnk)
        END DO
      END DO
    END DO
	
	!print *,"X Matrix"
    !do row=1,matrixsize
    !  do column=1,2 *matrixsize
	!	print *,row,column,xmatrix(row,column)
    !  end do
    !end do
  
    !recreate identity on left
    DO row=1,matrixsize
      IF(xmatrix(row,row).NE.1)THEN
        DO column=1,(2*matrixsize)
          IF(row.EQ.column)THEN
            idfactor = xmatrix(row,column)
            xmatrix(row,column) = 1.0D0
          ELSE
            IF(xmatrix(row,column).EQ.0)THEN
              xmatrix(row,column) = 0.0D0
            ELSE
              xmatrix(row,column) = (1.0D0 * xmatrix(row,column)) / &
			  (1.0D0 * idfactor)
            END IF
          END IF
        END DO
      END IF
    END DO

    !transfer inverted matrix
    DO row=1,matrixsize
      DO column=1,matrixsize
        xinvmatrix(row,column) = xmatrix(row,column+matrixsize)
      END DO
    END DO
	
	!mult ymatrix by xinvmatrix
    matrixarows = matrixsize
    matrixacolumns = matrixsize
    matrixbrows = matrixsize
    matrixbcolumns = 1

    DO column=1,matrixbcolumns
      DO row=1,matrixarows
        matrixentry = 0.0D0
        DO columnk=1,matrixacolumns
          matrixentry = 1.0D0 * matrixentry + 1.0D0 * xinvmatrix(row,columnk) * ymatrix(columnk,column)
        END DO
        coefficients(row,column) = matrixentry
      END DO
    END DO

	If(Allocated(polyCoefficients).eqv..true.)Then
	  Deallocate(polyCoefficients)
	End If
	Allocate(polyCoefficients(0:order))
	do row=1,matrixsize
	  polyCoefficients(row-1) = coefficients(row,1)
	  !print *,(row-1),coefficients(row,1)
	enddo
	
  
  !End of polyFit subroutine
  End Subroutine polyFit

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
End Module maths