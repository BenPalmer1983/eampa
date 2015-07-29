Module general

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! ----------------------------------------
! Updated: 1st May 2014
! ----------------------------------------

! Setup Modules
  Use kinds

! force declaration of all variables
  Implicit None
! Include MPI header
! Include 'mpif.h'
! Privacy of functions/subroutines/variables
  Private
! Public subroutines
  Public :: swapArrayRows1D, swapArrayRows2D
  Public :: extractArrayColumnDP, extractArrayColumnInt
  Public :: makeDir, rmFile, rmDir, randFileName, tempFileName
  Public :: strToIntArr, strToDPArr, strToStrArr
  Public :: timeAcc
  Public :: readFile
! Public functions
  Public :: dpToString, intToString
  Public :: GetClockTime
  Public :: StrToUpper
  Public :: NumericOnly
  Public :: RemoveSpaces
  Public :: CorrectFilePath
  Public :: TrimSpaces
  Public :: BlankString
  Public :: BlankStringArray
  Public :: BlankString2DArray
  Public :: SpacesRight
  Public :: RemoveComments
  Public :: RemoveQuotes

  character( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  Contains

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! List of Subroutines
! -------------------------------------------------------------------------
!

  Subroutine extractArrayColumnDP(inputArray,outputArray,column)
! Extract one column of a 2D dp array array(row,col)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, column
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
! Copy column
    Do i=1,size(inputArray,1)
      outputArray(i) = inputArray(i,column)
    End Do
  End Subroutine extractArrayColumnDP

  Subroutine extractArrayColumnInt(inputArray,outputArray)
! Subroutine extractArrayColumnInt(inputArray,outputArray,column)
! Extract one column of a 2D int array array(row,col)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: inputArray
    Integer(kind=StandardInteger), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
! Copy column
! Do i=1,size(inputArray,1)
!  outputArray(i) = inputArray(i,column)
! End Do
  End Subroutine extractArrayColumnInt

  Subroutine swapArrayRows1D(matrix,rowA,rowB)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
! Set variables
    matH = size(matrix,1)
    matW = 1
! Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
! Allocate arrays
      Allocate(rowAArr(1:matW))
      Allocate(rowBArr(1:matW))
! Swap rows
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
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
    Real(kind=DoubleReal), Dimension( : , :), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
! Set variables
    matH = size(matrix,1)
    matW = size(matrix,2)
! Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
! Allocate arrays
      Allocate(rowAArr(1:matW))
      Allocate(rowBArr(1:matW))
! Swap rows
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
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("mkdir -p "//trim(path))
  End Subroutine makeDir

  Subroutine rmFile(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("rm -f "//trim(path))
  End Subroutine rmFile

  Subroutine rmDir(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("rm -fR "//trim(path))
  End Subroutine rmDir

  Subroutine randFileName(fileName)
    Integer(kind=StandardInteger) :: i, characterNum
    Character(len=8) :: fileName
    Real(kind=DoubleReal) :: randNumber
    Character(len=52), Parameter :: alpha = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    fileName = "tmp     "
    Do i=4,8
      Call RANDOM_NUMBER(randNumber)
      characterNum = Ceiling(52.0E0*randNumber+1.0E0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.52)Then
        characterNum = 52
      End If
      fileName(i:i) = alpha(characterNum:characterNum)
    End Do
  End Subroutine randFileName

  Subroutine tempFileName(fileName)
    Integer(kind=StandardInteger) :: i, characterNum
    Character(len=8) :: fileName
    Real(kind=DoubleReal) :: randNumber
    Character(len=52), Parameter :: alpha = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    fileName = "tmp     "
    Do i=4,8
      Call RANDOM_NUMBER(randNumber)
      characterNum = Ceiling(52.0E0*randNumber+1.0E0)
      If(characterNum.lt.1)Then
        characterNum = 1
      End If
      If(characterNum.gt.52)Then
        characterNum = 52
      End If
      fileName(i:i) = alpha(characterNum:characterNum)
    End Do
  End Subroutine tempFileName

  Subroutine strToIntArr(stringIn,intArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger), Dimension(:) :: intArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(32) :: intString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    intString = BlankString(intString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        Read(intString,*) intArr(k)
        intString = BlankString(intString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        intString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToIntArr

  Subroutine strToDPArr(stringIn,dpArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal), Dimension(:) :: dpArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(32) :: intString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    intString = BlankString(intString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        Read(intString,*) dpArr(k)
        intString = BlankString(intString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        intString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToDPArr

  Subroutine strToStrArr(stringIn,strArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Character(*), Dimension(:) :: strArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(Len(strArr)) :: tempString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    tempString = BlankString(tempString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        strArr(k) = tempString
        tempString = BlankString(tempString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        tempString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToStrArr

! Time Accumulator subroutine
  Subroutine timeAcc(time,timeStart,timeEnd)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal) :: time,timeStart,timeEnd
    time = time + timeEnd - timeStart
  End Subroutine timeAcc
  
  
! ---------------------------------------------------------------------------------------------------
  Subroutine readFile(inputFilePath, fileArray, n)
! Subroutine to read file into an array
! Removes comments !.....
! Removes blank lines
! Removes leading spaces
    Implicit None ! Force declaration of all variables
! Private variables
    Character(*) :: inputFilePath
    Integer(kind=StandardInteger) :: i, j, n, ios
    Character(*), Dimension(:) :: fileArray
    Character(len=255) :: fileRow, fileRowTemp
    Logical :: commentsFlag
! open file
    Open(UNIT=4323,FILE=trim(inputFilePath))
    n = 0
    Do i=1,size(fileArray,1)
      Read(4323,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
! remove comments
      commentsFlag = .false.
      Do j=1,255
        If(fileRow(j:j).eq."!")Then
          commentsFlag = .true.
        End If
        If(commentsFlag)Then
          fileRow(j:j) = " "
        End If
      End Do  
! remove blank lines
      fileRowTemp = trim(adjustl(fileRow))
      If(fileRowTemp(1:1).ne." ")Then
        n = n + 1
        fileArray(n) = fileRowTemp
      End If
    End Do  
! Close file    
    close(4323)
  End Subroutine readFile


! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! List of Functions
! -------------------------------------------------------------------------
!

  Function dpToString(inputDP) RESULT (outputString)
! force declaration of all variables
    Implicit None
! declare private variables
    Real(kind=DoubleReal) :: inputDP
    Character(len=32) :: outputString
! Read dp to string
    inputDP = 1.0D0 * inputDP
    Write(outputString,"(ES16.8E3)") inputDP
  End Function dpToString

  Function intToString(inputInt) RESULT (outputString)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: inputInt
    Character(len=32) :: outputString
! Read int to string
    Write(outputString,"(I16)") inputInt
  End Function intToString

  Function GetClockTime () RESULT (outputTime)
! -- Argument and result
    Real(kind=DoubleReal) :: outputTime
    Call cpu_time(outputTime)
  End Function GetClockTime

  Function StrToUpper (input) RESULT (output)
! -- Argument and result
    CHARACTER(*), INTENT(IN) :: input
    CHARACTER(LEN(input)) :: output
! -- Local variables
    Integer(kind=StandardInteger) :: i, n
! -- Copy input string
    output = input
! -- Loop over string elements
    Do i = 1, LEN( output )
! -- Find location of letter in lower case constant string
      n = INDEX( LOWER_CASE, output( i:i ) )
! -- If current substring is a lower case letter, make it upper case
      If( n /= 0 ) output( i:i ) = UPPER_CASE( n:n )
      End Do
    End Function StrToUpper

    Function NumericOnly (input) RESULT (output)
! -- Argument and result
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      Integer(kind=StandardInteger) :: i, n
! -- Copy input string
      outputTemp = input
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
      n = 0
      Do i = 1, LEN( outputTemp )
        If(outputTemp( i:i ).eq.".".or.(iachar(outputTemp( i:i )).ge.48.and.iachar(outputTemp( i:i )).le.57))Then
          n = n + 1
          output( n:n ) = outputTemp( i:i )
        Else
          output( i:i ) = " "
        End If
      End Do
    End Function NumericOnly

    Function RemoveSpaces (input) RESULT (output)
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      Integer(kind=StandardInteger) :: i, j
! -- Copy input string
      outputTemp = input
! Blank output
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
! transfer outputtemp to output without spaces
      j = 0
      Do i = 1, LEN( outputTemp )
        If(outputTemp( i:i ).ne." ")Then
          j = j + 1
          output( j:j ) = outputTemp( i:i )
        End If
      End Do
    End Function RemoveSpaces

    Function TrimSpaces (input) RESULT (output)
      Character(*), INTENT(IN) :: input
      Character(LEN(trim(adjustl(input)))) :: output
      output = trim(adjustl(input))
    End Function TrimSpaces

    Function BlankString (input) RESULT (output)
      Character(*), INTENT(IN) :: input
      Character(Len(input)) :: output
      Integer(kind=StandardInteger) :: i
      Do i=1,Len(input)
        output(i:i) = " "
      End Do
    End Function BlankString

    Function BlankStringArray (input) RESULT (output)
      Character(*), Dimension(:), INTENT(IN) :: input
      Character(Len(input)) :: line
      Character(Len(input)), Dimension(1:size(input,1)) :: output
      Integer(kind=StandardInteger) :: i
      Do i=1,Len(input)
        line(i:i) = " "
      End Do
      Do i=1,size(input,1)
        output(i) = line
      End Do
    End Function BlankStringArray

    Function BlankString2DArray (input) RESULT (output)
      Character(*), Dimension(:,:), INTENT(IN) :: input
      Character(Len(input)) :: line
      Character(Len(input)), Dimension(1:size(input,1),1:size(input,2)) :: output
      Integer(kind=StandardInteger) :: i, j
      Do i=1,Len(input)
        line(i:i) = " "
      End Do
      Do i=1,size(input,1)
        Do j=1,size(input,2)
          output(i,j) = line
        End Do
      End Do
    End Function BlankString2DArray

    Function SpacesRight (input) RESULT (output)
! Adds spaces to right of string
      Character(*), INTENT(IN) :: input
      Character(Len(input)) :: tempStr, output
      Integer(kind=StandardInteger) :: i
      tempStr = trim(adjustl(input))
      Do i=1,Len(tempStr)
        output(i:i) = " "
      End Do
      Do i=1,Len(trim(tempStr))
        output(i:i) = tempStr(i:i)
      End Do
    End Function SpacesRight

    Function RemoveComments (input) RESULT (output)
! Removes comments from
      Character(*), INTENT(IN) :: input
      Character(Len(input)) :: output
      Integer(kind=StandardInteger) :: i
      output = BlankString(output)
! Comment character !
      Do i=1,Len(input)
        If(input(i:i).eq."!")Then
          Exit
        Else
          output(i:i) = input(i:i)
        End If
      End Do
    End Function RemoveComments

    Function RemoveQuotes (input) RESULT (output)
! Removes comments from
      Character(*), INTENT(IN) :: input
      Character(Len(input)) :: output
      Integer(kind=StandardInteger) :: i, j
      output = BlankString(output)
! Comment character !
      j = 0
      Do i=1,Len(input)
        If(ichar(input(i:i)).eq.34.or.ichar(input(i:i)).eq.39)Then
! Do nothing
        Else
          j = j + 1
          output(j:j) = input(i:i)
        End If
      End Do
    End Function RemoveQuotes

    Function CorrectFilePath (input) RESULT (output)
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! -- Local variables
      Integer(kind=StandardInteger) :: i, n
! -- Copy input string
      outputTemp = input
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
      n = 0
      Do i = 1, LEN( outputTemp )
        If((iachar(outputTemp( i:i )).ge.32.and.iachar(outputTemp( i:i )).le.126))Then
          If(outputTemp( i:i ).ne."?".and.outputTemp( i:i ).ne."*".and.&
            outputTemp( i:i ).ne."%".and.outputTemp( i:i ).ne."+")Then
            n = n + 1
            output( n:n ) = outputTemp( i:i )
          End If
        End If
      End Do
      output = trim (output)
    End Function CorrectFilePath

  End Module general
