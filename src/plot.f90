Module plot
! --------------------------------------------------------------!
! Plot
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Uses Matplotlib python library to build charts
! ----------------------------------------
! Updated: 17th June 2015
! ----------------------------------------
! Setup Modules
  Use kinds
  Use types
  Use maths
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
!  
  
! Public Subroutines
  Public :: makePlot
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine makePlot(outputDirectory, outputName, tempDirectory,&
    dataArrayIn, colStart, colEnd, rowStart, rowEnd, chartInput, polyFitFlag_in)
    Implicit None  ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, iA, jA
    Integer(kind=StandardInteger) :: colStart, colEnd, rowStart, rowEnd, colEndA, rowEndA
    Character(*) :: outputDirectory, tempDirectory, outputName
    Real(kind=DoubleReal), Dimension( : , : ) :: dataArrayIn
    Real(kind=DoubleReal), Dimension(1:(rowEnd-rowStart+1),1:(colEnd-colStart+1)) :: dataArray
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: polyFitData
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsA
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsB
    Real(kind=DoubleReal), Dimension(1:5) :: coefficientsC
    Real(kind=DoubleReal), Dimension(1:6) :: coefficientsD
    Character(len=8) :: fileName
    Character(len=14) :: dpStr, dpStrA, dpStrB
    Character(len=32768) :: xLine, yLine
    Character(len=65536) :: outputLine
    Integer(kind=StandardInteger) :: termExitStat
    Type(chart) :: chartInput
    Logical, Optional :: polyFitFlag_in
    Logical :: polyFitFlag
! Optional variables    
    polyFitFlag = .false.
    If(Present(polyFitFlag_in))Then
      polyFitFlag = polyFitFlag_in
    End If
! Prepare data array
    Do i=rowStart,rowEnd
      iA = i-rowStart+1
      Do j=colStart,colEnd
        jA = j-colStart+1
        dataArray(iA,jA) = dataArrayIn(i,j)
      End Do
    End Do 
! Reset end
    rowEndA = rowEnd-rowStart+1   
    colEndA = colEnd-colStart+1
! If polyFit selected  
    If(polyFitFlag)Then
      If(rowEndA.eq.3)Then
        coefficientsA = PolyFit(dataArray,2)
        polyFitData = PolyPoints(coefficientsA,dataArray(1,1),dataArray(rowEndA,1),100) 
      Else If(rowEndA.eq.4)Then
        coefficientsB = PolyFit(dataArray,3)
        polyFitData = PolyPoints(coefficientsB,dataArray(1,1),dataArray(rowEndA,1),100) 
      Else If(rowEndA.eq.5)Then
        coefficientsC = PolyFit(dataArray,4)
        polyFitData = PolyPoints(coefficientsC,dataArray(1,1),dataArray(rowEndA,1),100) 
      Else If(rowEndA.eq.6)Then
        coefficientsD = PolyFit(dataArray,5)
        polyFitData = PolyPoints(coefficientsD,dataArray(1,1),dataArray(rowEndA,1),100) 
      End If
    End If  
    !PolyPoints(coefficients,xStart,xEnd,points)
! Prepare x-values
    xLine = BlankString(xLine)
! Insert x-values
    If(polyFitFlag)Then
      Do i=1,100 
        write(dpStr,"(E14.6)") polyFitData(i,1)
        xLine((i-1)*15+1:(i-1)*15+14) = dpStr
        If(i.lt.100)Then
          xLine((i-1)*15+15:(i-1)*15+15) = ","
        End If  
      End Do 
    Else
      Do i=1,rowEndA 
        write(dpStr,"(E14.6)") dataArray(i,1)
        xLine((i-1)*15+1:(i-1)*15+14) = dpStr
        If(i.lt.rowEndA)Then
         xLine((i-1)*15+15:(i-1)*15+15) = ","
        End If  
      End Do 
    End If  
! Make temp random name
    fileName = RandName()    
! Open file
    open(unit=701,file=(trim(tempDirectory)//"/"//fileName//".py"))
! write python headers
    write(701,"(A)") "#!/usr/bin/env python"
    write(701,"(A)") "import numpy as np"
    write(701,"(A)") "import matplotlib"
    write(701,"(A)") "matplotlib.use('Agg')"
    write(701,"(A)") "import matplotlib.pyplot as plt"
! Set figure sizes
    write(701,"(A)") "plt.figure(figsize=(1792/144, 1008/144), dpi=144)"
! Data
    outputLine = BlankString(outputLine)
    yLine = BlankString(yLine)
! Insert y-values
    If(polyFitFlag)Then
      Do i=1,100 
        write(dpStr,"(E14.6)") polyFitData(i,2)
        yLine((i-1)*15+1:(i-1)*15+14) = dpStr
        If(i.lt.100)Then
          yLine((i-1)*15+15:(i-1)*15+15) = ","
        End If  
      End Do
    Else
      Do i=1,rowEndA 
        write(dpStr,"(E14.6)") dataArray(i,2)
        yLine((i-1)*15+1:(i-1)*15+14) = dpStr
        If(i.lt.rowEndA)Then
          yLine((i-1)*15+15:(i-1)*15+15) = ","
        End If  
      End Do 
    End If  
! Make plot line and write
    outputLine = "plt.plot(["//trim(xLine)//"],["//trim(yLine)//"],label=' ')"  
    write(701,"(A)") trim(RemoveSpaces(outputLine))   
! Write Titles
    write(701,"(A)") "plt.title('"//trim(chartInput%title)//"')"
! Axes Labels
    write(701,"(A)") "plt.xlabel('"//trim(chartInput%xAxis)//"')"
    write(701,"(A)") "plt.ylabel('"//trim(chartInput%yAxis)//"')"
! Resize axis
    If(chartInput%xMin.lt.1.0D99.and.chartInput%xMax.gt.-1.0D99)Then
      write(dpStrA,"(E14.6)") chartInput%xMin      
      write(dpStrB,"(E14.6)") chartInput%xMax      
      write(701,"(A)") "plt.xlim("//dpStrA//","//dpStrB//")"    
    End If
    If(chartInput%yMin.lt.1.0D99.and.chartInput%yMax.gt.-1.0D99)Then
      write(dpStrA,"(E14.6)") chartInput%yMin      
      write(dpStrB,"(E14.6)") chartInput%yMax      
      write(701,"(A)") "plt.ylim("//dpStrA//","//dpStrB//")"    
    End If
! Legend
    !write(701,"(A)") "plt.legend()"
    !write(701,"(A)") "plt.legend()"
! Set output file
    write(701,"(A)") "plt.savefig('"//trim(outputDirectory)//"/"//outputName//"',dpi=144)"
! Close file
    close(701)
! Run python and the file to create the chart
    Call execute_command_line("python "//trim(tempDirectory)//"/"//fileName//".py",&
    exitstat=termExitStat)
! Clean python file    
    If(chartInput%cleanPyFile)Then
      Call system("rm -f "//(trim(tempDirectory)//"/"//fileName//".py"))
    End If
  End Subroutine makePlot
  

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!
  
  Function BlankString(input) RESULT (output)
    Character(*), INTENT(IN) :: input
    Character(Len(input)) :: output
    Integer(kind=StandardInteger) :: i
    Do i=1,Len(input)
      output(i:i) = " "
    End Do
  End Function BlankString
  
  Function RandName() RESULT (fileName)
! force declaration of all variables
    Implicit None  
    Integer(kind=StandardInteger) :: i, characterNum
    Character(len=8) :: fileName
    Real(kind=DoubleReal) :: randNumber
    Character(len=52), Parameter :: alpha = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
! Prepare string
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
  End Function RandName
  
  Function RemoveSpaces(input)RESULT(output)
! force declaration of all variables
    Implicit None    
      CHARACTER(*), INTENT(IN) :: input
      CHARACTER(LEN(input)) :: outputTemp
      CHARACTER(LEN(input)) :: output
! Local variables
      Integer(kind=StandardInteger) :: i, j
! Copy input string
      outputTemp = input
! Blank output
      Do i = 1, LEN( outputTemp )
        output( i:i ) = " "
      End Do
! transfer outputtemp to output without spaces
      j = 0
      Do i = 1, LEN(outputTemp)
        If(outputTemp( i:i ).ne." ")Then
          j = j + 1
          output( j:j ) = outputTemp( i:i )
        End If
      End Do
  End Function RemoveSpaces
End Module plot
