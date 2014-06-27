Module general

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 1st May 2014
!----------------------------------------


! Setup Modules
  Use kinds

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!Privacy of functions/subroutines/variables
  Private    
!Public subroutines
  Public :: swapArrayRows1D, swapArrayRows2D
  Public :: extractArrayColumnDP, extractArrayColumnInt
  Public :: makeDir
!Public functions
  Public :: dpToString, intToString
  Public :: GetClockTime
  
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Subroutines
!-------------------------------------------------------------------------  
! 


  Subroutine extractArrayColumnDP(inputArray,outputArray,column) 
!Extract one column of a 2D dp array array(row,col)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, column
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
!Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
!Copy column
    Do i=1,size(inputArray,1)
	  outputArray(i) = inputArray(i,column)
	End Do
  End Subroutine extractArrayColumnDP
  
  Subroutine extractArrayColumnInt(inputArray,outputArray,column) 
  !Subroutine extractArrayColumnInt(inputArray,outputArray,column) 
!Extract one column of a 2D int array array(row,col)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, column
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: inputArray
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: outputArray
!Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
!Copy column
    !Do i=1,size(inputArray,1)
	!  outputArray(i) = inputArray(i,column)
	!End Do
  End Subroutine extractArrayColumnInt
  
  
    
  Subroutine swapArrayRows1D(matrix,rowA,rowB) 
!Swap rows of square dp matrix
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
!Set variables
	matH = size(matrix,1)
	matW = 1
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
!Allocate arrays
	  Allocate(rowAArr(1:matW))
	  Allocate(rowBArr(1:matW))
!Swap rows
	  Do i=1,matW
	    rowAArr(i) = matrix(rowA)
	    rowBArr(i) = matrix(rowB)
	  End Do
	  Do i=1,matW
	    matrix(rowA) = rowBArr(i)
	    matrix(rowB) = rowAArr(i)
	  End Do
    End If
  End Subroutine swapArrayRows1D
  
  
  
  Subroutine swapArrayRows2D(matrix,rowA,rowB) 
!Swap rows of square dp matrix
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
	Real(kind=DoubleReal), Dimension( : , :), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
!Set variables
	matH = size(matrix,1)
	matW = size(matrix,2)
!Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
!Allocate arrays
	  Allocate(rowAArr(1:matW))
	  Allocate(rowBArr(1:matW))
!Swap rows
	  Do i=1,matW
	    rowAArr(i) = matrix(rowA,i)
	    rowBArr(i) = matrix(rowB,i)
	  End Do
	  Do i=1,matW
	    matrix(rowA,i) = rowBArr(i)
	    matrix(rowB,i) = rowAArr(i)
	  End Do
    End If
  End Subroutine swapArrayRows2D
 
  
  Subroutine makeDir(path) 
!Swap rows of square dp matrix
!force declaration of all variables
	Implicit None	
!declare private variables
    Character(len=128) :: path
  End Subroutine makeDir
  
  
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 
 
  Function dpToString(inputDP) RESULT (outputString)
!force declaration of all variables
	Implicit None	
!declare private variables
    Real(kind=DoubleReal) :: inputDP
	Integer(kind=StandardInteger) :: stringLength
	Character(len=32) :: outputString
!Read dp to string
    inputDP = 1.0D0 * inputDP
	Write(outputString,"(ES16.8E3)") inputDP
  End Function dpToString
  
  Function intToString(inputInt) RESULT (outputString)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: inputInt
	Integer(kind=StandardInteger) :: stringLength
	Character(len=32) :: outputString
!Read int to string
	Write(outputString,"(I16)") inputInt
  End Function intToString
  
  
  
  
  Function GetClockTime () RESULT (outputTime)
    ! -- Argument and result
	Real(kind=DoubleReal) :: outputTime 
	Call cpu_time(outputTime)
  End Function GetClockTime    
  
  
End Module general  