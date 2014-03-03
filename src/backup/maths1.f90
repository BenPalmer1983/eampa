Module maths

! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
  Private  
  !functions
  Public :: solvePolynomial
  Public :: SolvePolynomialRand
  Public :: calcPolynomial
  Public :: polyFit
  Public :: polynomialFit
  Public :: QuadraticInterpolation
  Public :: QuadraticInterpolationCalc
  Public :: CubicInterpolationCalc
  
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
	Integer(kind=StandardInteger) i,j,k,maxLoops	
!Set values
	convergenceTarget = 0
	convergence = 1000
	convergenceThreshold = 0.00001
	maxLoops = 0
!set start value for x
	factor = 0.5
	difference = upper - lower
	x = lower + factor * difference	
	do while(convergence.gt.convergenceThreshold.and.maxLoops.le.10000)
	  maxLoops = maxLoops + 1
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
    
  
  
  
  Function SolvePolynomialRand (coefficients, lower, upper) RESULT (output)    	
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: upper, lower, output
	Real(kind=DoubleReal) :: x,y,dydx,yBest,xBest
    Real(kind=DoubleReal) :: convergence, convergenceThreshold, convergenceTarget
	Real(kind=DoubleReal) :: factor, difference, randNumber
	Integer(kind=StandardInteger) i,j,k,maxLoops  
!polynomial of the form ax^2+bx+c=0
    do i=1,1000000
!get random number	  
	  Call RANDOM_NUMBER(randNumber)
!set x value
	    x=lower+randNumber*(upper-lower)
!calculate y
        y=0.0D0
	    do j=0,size(coefficients)-1
	      y=y+1.0D0*x**(j)*coefficients(j)			
	    enddo
		y=abs(y)
!store best answer
      if(i.eq.1.or.y.lt.ybest)then
	    yBest = y
		xBest = x
	  endif
	enddo
    output = xBest
  End Function SolvePolynomialRand 
  
  
  
  function calcPolynomial (polyCoefficients, x, derivative) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k, derivative
	Integer(kind=StandardInteger) :: coeffCount 
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyTemp
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyWorking
    Real(kind=DoubleReal) :: x, y
!count coefficients
	coeffCount = size(polyCoefficients)	
!force no or positive derivative
	if(derivative.lt.0)then
	  derivative = 0
	endif
!run calc
	if(derivative.gt.coeffCount)then
!if order of derivative higher than coeffs/order polynomial, then answer is zero	
	  y = 0.0D0
	  Allocate(polyWorking(0:0))
!transfer to polyWorking
	  polyWorking(0) = 0.0D0
	  coeffCount = 1
	else
	  if(derivative.eq.0)then
!if no derivative, use input coeffs	
        Allocate(polyWorking(0:(coeffCount-1)))
!transfer to polyWorking
	    do i=0,(coeffCount-1)
		  polyWorking(i) = polyCoefficients(i)
		enddo
	  else
!loop through each derivative
	    do i=1,derivative
		  coeffCount = coeffCount - 1
		  do j=0,(coeffCount-1)
		    polyCoefficients(j) = (j+1)*polyCoefficients(j+1)
		  enddo
		enddo	
        Allocate(polyWorking(0:(coeffCount-1)))
!transfer to polyWorking
	    do i=0,(coeffCount-1)
		  polyWorking(i) = polyCoefficients(i)
		enddo
	  endif
	  y = 0.0D0
	  do i=0,(coeffCount-1)
	    y = y+1.0D0*polyWorking(i)*x**(1.0D0*i)
	  enddo 
	endif 
!returns y
  End Function calcPolynomial
  
  
  
  
  Function polynomialFit(inputPoints,order) RESULT (polyCoefficients) 
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i,j,k,l,m,x,y,z,verbose
	Integer(kind=StandardInteger) :: order, exponentInt, optimiseSum
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputPoints
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: polyCoefficients
	Integer(kind=StandardInteger) :: dataPoints, matrixSize
	Integer(kind=StandardInteger) :: xColMax, xRowMax, yRowMax
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yMatrix
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: xMatrixRow
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: xMatrixInverse
	Real(kind=DoubleReal) :: exponentValue, matrixSumEntry
	Real(kind=DoubleReal) :: xA, xB, triangleFactor, tempDouble, inverseFactor
	Logical :: optimising
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: optimiseArray	
!i=column, j=row
!set values
    verbose = 0
	dataPoints = size(inputPoints)/2
	matrixSize = order + 1	
!Allocate arrays
    Allocate(polyCoefficients(0:order))
    Allocate(xMatrix(1:2*matrixSize,1:matrixSize))
    Allocate(xMatrixRow(1:2*matrixSize))	
    Allocate(yMatrix(1:matrixSize))	
    Allocate(optimiseArray(1:matrixSize))
    Allocate(xMatrixInverse(1:matrixSize,1:matrixSize))	
!Fill arrays
	do i=1,(2*matrixSize)
	  do j=1,matrixSize
	    xMatrix(i,j) = 0.0D0
	  enddo
	enddo
	do i=1,matrixSize
	  yMatrix(i) = 0.0D0
	enddo	
	do i=1,matrixSize
	  xMatrixRow(i) = 0.0D0
	enddo	
	do i=1,matrixSize
	  do j=1,matrixSize
	    xMatrixInverse(i,j) = 0.0D0
	  enddo
	enddo		
!Build Least Squares Fitting Vandermonde matrix
    do i=1,matrixSize
	  do j=1,matrixSize
	    exponentValue = 1.0D0*i+1.0D0*j-2.0D0
		do k=1,dataPoints
	      xMatrix(i,j) = xMatrix(i,j)+1.0D0*inputPoints(k,1)**exponentValue
		enddo
	  enddo
	enddo
	do j=1,matrixSize
	  exponentValue = 1.0D0*j-1.0D0
	  do k=1,dataPoints
	    yMatrix(j) = yMatrix(j)+1.0D0*inputPoints(k,2)*inputPoints(k,1)**exponentValue
	  enddo
	enddo	
!print start
    if(verbose.eq.1)then
	  do j=1,matrixSize
	    print *,xMatrix(1,j),xMatrix(2,j),xMatrix(3,j),xMatrix(4,j),xMatrix(5,j),&
	    xMatrix(6,j),xMatrix(7,j),xMatrix(8,j)
	  enddo
	  print *,""
	  print *,""
	endif  
!print end
!optimise matrix order
    do j=1,matrixSize		!row at a time
	  optimiseSum = 0
	  do i=1,matrixSize 
	    if(xMatrix(i,j).ne.0.0D0)then
	      exponentInt = matrixSize - i
		  optimiseSum = optimiseSum + 2**exponentInt
		endif
	  enddo
      optimiseArray(j) = optimiseSum 
	enddo	
	optimising = .true.
    do while(optimising)
      optimising = .false.
      do j=1,(matrixSize-1)        !loop through rows
        if(optimiseArray(j).lt.optimiseArray(j+1))then
		  optimising = .true.
!reorder optimising array
          xA = optimiseArray(j)
          xB = optimiseArray(j+1)
		  optimiseArray(j) = xB
		  optimiseArray(j+1) = xA
!reorder xMatrix
          do i=1,(2*matrixSize)   
		    xA = 1.0D0*xMatrix(i,j)
            xB = 1.0D0*xMatrix(i,j+1)
		    xMatrix(i,j) = 1.0D0*xB
		    xMatrix(i,j+1) = 1.0D0*xA
		  enddo
		endif
	  enddo
    enddo	  	
!append identity array to xMatrix
	do i=1,matrixSize
      do j=1,matrixSize
	    if(i.eq.j)then
	      xMatrix(i+matrixSize,j) = 1.0D0
		endif
	  enddo
	enddo
!print start
    if(verbose.eq.1)then
	  do j=1,matrixSize
	    print *,xMatrix(1,j),xMatrix(2,j),xMatrix(3,j),xMatrix(4,j),xMatrix(5,j),&
	    xMatrix(6,j),xMatrix(7,j),xMatrix(8,j)
	  enddo
	  print *,""
	  print *,""
	endif  
!print end	
!make lower triangle of zeros
    do j=1,matrixSize-1
	  do k=j+1,matrixSize
        triangleFactor = 1.0D0*((1.0D0*xMatrix(j,j))/(1.0D0*xMatrix(j,k)))
		!calculate new row values
		do i=1,(2*matrixSize) !loop over all columns
		  if(xMatrix(i,k).eq.0.0D0.and.i.le.matrixSize)then
		    xMatrixRow(i) = 0.0D0
		  else	
		    xMatrixRow(i) = 1.0D0*&
		    ((1.0D0*xMatrix(j,j))/(1.0D0*xMatrix(j,k)))*xMatrix(i,k)-1.0D0*xMatrix(i,j)
		  endif	
		enddo
!replace row values
		do i=1,(2*matrixSize) !loop over all columns
		  xMatrix(i,k) = xMatrixRow(i)
		enddo
	  enddo
!force zeros in the lower triangle
      do k=j+1,matrixSize
	    xMatrix(j,k) = 0.0D0
	  enddo
	enddo
!print start
    if(verbose.eq.1)then
	  do j=1,matrixSize
	    print *,xMatrix(1,j),xMatrix(2,j),xMatrix(3,j),xMatrix(4,j),xMatrix(5,j),&
	    xMatrix(6,j),xMatrix(7,j),xMatrix(8,j)
	  enddo
	  print *,""
	  print *,""
	endif  
!print end
!make upper triangle of zeros
    do j=matrixSize,2,-1
	  do k=j-1,1,-1
		do i=1,(2*matrixSize) !loop over all columns
		  if(xMatrix(i,k).eq.0.0D0.and.i.le.matrixSize)then
		    xMatrixRow(i) = 0.0D0
		  else	
            xMatrixRow(i) = 1.0D0*&
			((1.0D0*xMatrix(j,j))/(1.0D0*xMatrix(j,k)))*xMatrix(i,k)-1.0D0*xMatrix(i,j)			
		  endif	
		enddo
!replace row values
		do i=1,(2*matrixSize) !loop over all columns
		  xMatrix(i,k) = xMatrixRow(i)
		enddo
	  enddo
!force zeros in the lower triangle
      do k=j-1,1,-1
	    xMatrix(j,k) = 0.0D0
	  enddo	
	enddo	
!print start
    if(verbose.eq.1)then
	  do j=1,matrixSize
	    print *,xMatrix(1,j),xMatrix(2,j),xMatrix(3,j),xMatrix(4,j),xMatrix(5,j),&
	    xMatrix(6,j),xMatrix(7,j),xMatrix(8,j)
	  enddo
	  print *,""
	  print *,""
	endif  
!print end
!make identity on the left
    do j=1,matrixSize		!loop through rows
	  inverseFactor = 1.0D0*xMatrix(j,j)
	  do i=1,2*matrixSize
		xMatrix(i,j) = (1.0D0*xMatrix(i,j))/inverseFactor
	  enddo
	  !xMatrix(j,j) = 1.0D0
	enddo
!print start
    if(verbose.eq.1)then
	  do j=1,matrixSize
	    print *,xMatrix(1,j),xMatrix(2,j),xMatrix(3,j),xMatrix(4,j),xMatrix(5,j),&
	    xMatrix(6,j),xMatrix(7,j),xMatrix(8,j)
	  enddo
	  print *,""
	  print *,""
	endif  
!print end
!transfer right side of xMatrix to xMatrixInverse
  do j=1,matrixSize		!loop through rows
    do i=1,matrixSize
	  xMatrixInverse(i,j) = 1.0D0*xMatrix(i+matrixSize,j)
    enddo
  enddo
!print start
    if(verbose.eq.1)then
	do j=1,matrixSize
	  print *,xMatrixInverse(1,j),xMatrixInverse(2,j),xMatrixInverse(3,j),&
	  xMatrixInverse(4,j)
	enddo
	print *,""
	print *,""
	do j=1,matrixSize
	  print *,yMatrix(j)
	enddo
	print *,""
	print *,""
	endif
!print end

!multiply inverse with yMatrix	
	polyCoefficients = matMul(xMatrixInverse,yMatrix)
!End of polyFit subroutine
  End Function polynomialFit
  
  
  
  
  
  
  
  
  
    
    
      
  Function printMatrix(matrixIn) RESULT (x)    
	Integer(kind=StandardInteger) :: x,i,j
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: matrixIn	
	x = 1	
	do j=1,size(matrixIn,2)	!row
	  do i=1,size(matrixIn,1) !column
	    print *,matrixIn(i,j)
	  enddo
	enddo
  End Function printMatrix
  
  
  Function QuadraticInterpolation(xA,yA,xB,yB,xC,yC) RESULT (coefficients)
!force declaration of all variables
	Implicit None
!declare variables
	Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC
	Real(kind=DoubleReal) :: aA,aB,aC
	Allocate(coefficients(0:2))  	
	aA = 1.0D0*(yA/((xA-xB)*(xA-xC)))
    aB = 1.0D0*(yB/((xB-xA)*(xB-xC)))
	aC = 1.0D0*(yC/((xC-xA)*(xC-xB)))
	coefficients(2) = 1.0D0*(aA+aB+aC)
	coefficients(1) = -1.0D0*(aA*(xB+xC)+aB*(xA+xC)+aC*(xA+xB))
	coefficients(0) = 1.0D0*(aA*xB*xC+aB*xA*xC+aC*xA*xB)
  End Function QuadraticInterpolation
  !------------------------------
  Function QuadraticInterpolationCalc(xA,yA,xB,yB,xC,yC,x) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC,x,y
	Allocate(coefficients(0:2)) 
	coefficients = QuadraticInterpolation(xA,yA,xB,yB,xC,yC)
	y = 1.0D0*coefficients(0)+1.0D0*coefficients(1)*x+1.0D0*coefficients(2)*x**2
  End Function QuadraticInterpolationCalc
  
    
  
  Function CubicInterpolationCalc(xA,yA,xB,yB,xC,yC,xD,yD,x) RESULT (y)
!force declaration of all variables
	Implicit None
!declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: coefficients
	Real(kind=DoubleReal) :: xA,yA,xB,yB,xC,yC,xD,yD,x,y
	Real(kind=DoubleReal) :: aA,aB,aC,aD
	
	aA = 1.0D0*yA*(((x-xB)*(x-xC)*(x-xD))/((xA-xB)*(xA-xC)*(xA-xD)))
	aB = 1.0D0*yB*(((x-xA)*(x-xC)*(x-xD))/((xB-xA)*(xB-xC)*(xB-xD)))
	aC = 1.0D0*yC*(((x-xA)*(x-xB)*(x-xD))/((xC-xA)*(xC-xB)*(xC-xD)))
	aD = 1.0D0*yD*(((x-xA)*(x-xB)*(x-xC))/((xD-xA)*(xD-xB)*(xD-xC)))
		
	y = 1.0D0*(aA+aB+aC+aD)
  End Function CubicInterpolationCalc
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
  

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