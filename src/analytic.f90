Module analytic
! --------------------------------------------------------------!
! Description
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
!
! ----------------------------------------
! Updated: 15th Dec 2015
! ----------------------------------------
! Setup Modules
  Use kinds
  Use types
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use initialise
  Use loadData
  Use globals
  Use readEAM
  Use output
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: readAnalytic
  Public :: updateAnalytic
  Public :: varyAnalytic
  Public :: varyAnalyticRandom
! Public Functions
  !Public :: publicFunction

  Contains

! --------------------------------------------------------------------------------------------------!
!                                                                                                   !
! MODULE SUBROUTINES                                                                                  !
!                                                                                                   !
!                                                                                                   !
! --------------------------------------------------------------------------------------------------!

  Subroutine readAnalytic()
! Read in parameters
    Implicit None   ! Force declaration of all variables
! Private
    Character(len=128) :: generatedPotFile


    generatedPotFile = trim(outputDirectory)//"/generated.pot"
! Read in analytic potential file
    Call readAnalyticFile()
! Make a pot file to read in
    Call makePotentialFile(generatedPotFile)
! Store input apf
    apfDataIn = apfData

    eamFilePath = trim(outputDirectory)//"/generated.pot"

    eamFileType = 1
    Call readEAMFile()             ! readEAM.f90

  End Subroutine readAnalytic

! --------------------------------------------------------------------------------------------------!

  Subroutine readAnalyticFile()
! Read in parameters
    Implicit None   ! Force declaration of all variables
! Private
    Character(len=128), Dimension(1:50) :: analyticFile
    Character(len=128) :: fileRow
    Character(len=16) :: bufferA, bufferB, bufferC, bufferD
    Integer(kind=StandardInteger) :: fileRows
    Integer(kind=StandardInteger) :: i, j, fCount
    Integer :: functionCounter, paramID
! Read file
    Call readFile(trim(eamFilePath), analyticFile, fileRows)
! convert to upper case
    Do i=1,fileRows
      analyticFile(i) = StrToUpper(analyticFile(i))
    End Do
! read through
    i = 0
    functionCounter = 0
! Start loop
    Do While(i.lt.fileRows)
      i = i + 1
      fileRow = analyticFile(i)
! -----------
! PAIR
! -----------
      If(fileRow(1:4).eq."PAIR")Then
        functionCounter = functionCounter + 1
! Pair
        apfData(functionCounter)%functionType = "PAIR"
        Read(fileRow,*) bufferA, bufferB, bufferC
        apfData(functionCounter)%elementA = bufferB
        apfData(functionCounter)%elementB = bufferC
! point count row
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%dataPointCount
! fmin
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%fMin
! fmax
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%fMax
! pair zbl
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%pairZBL
        If(apfData(functionCounter)%pairZBL.gt.0.0D0)Then
          Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
          Read(bufferC,*) apfData(functionCounter)%pairZBL_A
          Read(bufferD,*) apfData(functionCounter)%pairZBL_B
        End If
! function form
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        bufferA = StrToUpper(bufferA)
        If(bufferA(1:4).eq."TYPE")Then
          !---------------------------------------------------------------------
          If(bufferB.eq."LJ")Then  ! [1] Lennard Jones, 2 parameters
          !---------------------------------------------------------------------
            apfData(functionCounter)%functionForm = 1
            apfData(functionCounter)%functionParameterCount = 2
            ! read sigma value
            i = i + 1
            fileRow = analyticFile(i)
            fCount = CountRowFields(fileRow)
            If(fCount.eq.1)Then
              Read(fileRow,*) bufferA
              Read(bufferA,*) apfData(functionCounter)%functionParameters(1)
            End If
            If(fCount.eq.3)Then
              Read(fileRow,*) bufferA, bufferB, bufferC
              Read(bufferA,*) apfData(functionCounter)%functionParameters(1)
              Read(bufferB,*) apfData(functionCounter)%functionParametersLB(1)
              Read(bufferC,*) apfData(functionCounter)%functionParametersUB(1)
            End If
            ! read epsilon value
            i = i + 1
            fileRow = analyticFile(i)
            fCount = CountRowFields(fileRow)
            If(fCount.eq.1)Then
              Read(fileRow,*) bufferA
              Read(bufferA,*) apfData(functionCounter)%functionParameters(2)
            End If
            If(fCount.eq.3)Then
              Read(fileRow,*) bufferA, bufferB, bufferC
              Read(bufferA,*) apfData(functionCounter)%functionParameters(2)
              Read(bufferB,*) apfData(functionCounter)%functionParametersLB(2)
              Read(bufferC,*) apfData(functionCounter)%functionParametersUB(2)
            End If
          End If
          !---------------------------------------------------------------------
          If(bufferB.eq."MORSE")Then  ! [4] Morse, 3 parameters
          !---------------------------------------------------------------------
            apfData(functionCounter)%functionForm = 4
            apfData(functionCounter)%functionParameterCount = 3
            Do paramID=1,apfData(functionCounter)%functionParameterCount
              ! read sigma value
              i = i + 1
              fileRow = analyticFile(i)
              fCount = CountRowFields(fileRow)
              If(fCount.eq.1)Then
                Read(fileRow,*) bufferA
                Read(bufferA,*) apfData(functionCounter)%functionParameters(paramID)
              End If
              If(fCount.eq.3)Then
                Read(fileRow,*) bufferA, bufferB, bufferC
                Read(bufferA,*) apfData(functionCounter)%functionParameters(paramID)
                Read(bufferB,*) apfData(functionCounter)%functionParametersLB(paramID)
                Read(bufferC,*) apfData(functionCounter)%functionParametersUB(paramID)
              End If
            End Do
          End If
        End If
      End If
! -----------
! DENS + EMBE
! -----------
      If(fileRow(1:4).eq."DENS".or.fileRow(1:4).eq."EMBE")Then
        functionCounter = functionCounter + 1
! DENS
        apfData(functionCounter)%functionType = fileRow(1:4)
        Read(fileRow,*) bufferA, bufferB
        apfData(functionCounter)%elementA = bufferB
! point count row
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%dataPointCount
! fmin
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%fMin
! fmax
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) apfData(functionCounter)%fMax
! function form
        i = i + 1
        fileRow = analyticFile(i)
        Read(fileRow,*) bufferA, bufferB
        !---------------------------------------------------------------------
        If(bufferB.eq."BASICDENSITY")Then  ! [2] Basic Density, 2 parameters
        !---------------------------------------------------------------------
          apfData(functionCounter)%functionForm = 2
          apfData(functionCounter)%functionParameterCount = 2
          i = i + 1
          fileRow = analyticFile(i)
          fCount = CountRowFields(fileRow)
          If(fCount.eq.1)Then
            Read(fileRow,*) bufferA
            Read(bufferA,*) apfData(functionCounter)%functionParameters(1)
          End If
          If(fCount.eq.3)Then
            Read(fileRow,*) bufferA, bufferB, bufferC
            Read(bufferA,*) apfData(functionCounter)%functionParameters(1)
            Read(bufferB,*) apfData(functionCounter)%functionParametersLB(1)
            Read(bufferC,*) apfData(functionCounter)%functionParametersUB(1)
          End If
          print *,"basic density",fCount,apfData(functionCounter)%functionParameters(1)
          i = i + 1
          fileRow = analyticFile(i)
          fCount = CountRowFields(fileRow)
          If(fCount.eq.1)Then
            Read(fileRow,*) bufferA
            Read(bufferA,*) apfData(functionCounter)%functionParameters(2)
          End If
          If(fCount.eq.3)Then
            Read(fileRow,*) bufferA, bufferB, bufferC
            Read(bufferA,*) apfData(functionCounter)%functionParameters(2)
            Read(bufferB,*) apfData(functionCounter)%functionParametersLB(2)
            Read(bufferC,*) apfData(functionCounter)%functionParametersUB(2)
          End If
        End If
        !---------------------------------------------------------------------
        If(bufferB.eq."BASICEMBED")Then  ! [3] Basic Embed, 3 parameters
        !---------------------------------------------------------------------
          apfData(functionCounter)%functionForm = 3
          apfData(functionCounter)%functionParameterCount = 3
          i = i + 1
          fileRow = analyticFile(i)
          fCount = CountRowFields(fileRow)
          If(fCount.eq.1)Then
            Read(fileRow,*) bufferA
            Read(bufferA,*) apfData(functionCounter)%functionParameters(1)
          End If
          If(fCount.eq.3)Then
            Read(fileRow,*) bufferA, bufferB, bufferC
            Read(bufferA,*) apfData(functionCounter)%functionParameters(1)
            Read(bufferB,*) apfData(functionCounter)%functionParametersLB(1)
            Read(bufferC,*) apfData(functionCounter)%functionParametersUB(1)
          End If
          i = i + 1
          fileRow = analyticFile(i)
          fCount = CountRowFields(fileRow)
          If(fCount.eq.1)Then
            Read(fileRow,*) bufferA
            Read(bufferA,*) apfData(functionCounter)%functionParameters(2)
          End If
          If(fCount.eq.3)Then
            Read(fileRow,*) bufferA, bufferB, bufferC
            Read(bufferA,*) apfData(functionCounter)%functionParameters(2)
            Read(bufferB,*) apfData(functionCounter)%functionParametersLB(2)
            Read(bufferC,*) apfData(functionCounter)%functionParametersUB(2)
          End If
          i = i + 1
          fileRow = analyticFile(i)
          fCount = CountRowFields(fileRow)
          If(fCount.eq.1)Then
            Read(fileRow,*) bufferA
            Read(bufferA,*) apfData(functionCounter)%functionParameters(3)
          End If
          If(fCount.eq.3)Then
            Read(fileRow,*) bufferA, bufferB, bufferC
            Read(bufferA,*) apfData(functionCounter)%functionParameters(3)
            Read(bufferB,*) apfData(functionCounter)%functionParametersLB(3)
            Read(bufferC,*) apfData(functionCounter)%functionParametersUB(3)
          End If
        End If
      End If
! End loop
    End Do
! Store function count
    apfData(1)%functionCount = functionCounter

    Call outputEAMAnParams(apfData)

    print *,functionCounter


  End Subroutine readAnalyticFile

! --------------------------------------------------------------------------------------------------!

  Subroutine makePotentialFile(generatedPotFile)
! Read in parameters
    Implicit None   ! Force declaration of all variables
! Private
    Character(len=128) :: generatedPotFile
    Integer(kind=StandardInteger) :: functionID
    Integer(kind=StandardInteger) :: i

    open(unit=999,file=trim(generatedPotFile))
    write(999,"(A4)") "#EAM"

    Do functionID=1,apfData(1)%functionCount
! Make points
      Call makePoints(functionID)
      write(999,"(A4,A1,A4,A1,A4)") apfData(functionID)%functionType," ",&
      apfData(functionID)%elementA," ",apfData(functionID)%elementB
      Do i=1,apfData(functionID)%dataPointCount
        write(999,"(E16.8,A1,E16.8)") apfData(functionID)%dataPoints(i,1)," ",&
        apfData(functionID)%dataPoints(i,2)
      End Do
    End Do
! Close file
    close(999)

  End Subroutine makePotentialFile

! --------------------------------------------------------------------------------------------------!

  Subroutine makePoints(functionID)
! Read in parameters
    Implicit None   ! Force declaration of all variables
! Private
    Integer(kind=StandardInteger) :: functionID
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: x, y, zblS, zblM, zblE, xB, xBd
    Real(kind=DoubleReal) :: xIncrement
    Real(kind=DoubleReal), Dimension(1:3) :: yArray
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients

    x = apfData(functionID)%fMin
    xIncrement = (apfData(functionID)%fMax-apfData(functionID)%fMin)/&
    (1.0D0*(apfData(functionID)%dataPointCount-1.0D0))



! Spline coefficients from ZBL to pair function
    zblE = apfData(functionID)%pairZBL
    If(apfData(functionID)%pairZBL.gt.0.0D0)Then
      zblS = 0.0D0
      zblM = 0.8D0*zblE
      yArray = ZblFull(zblM,apfData(functionID)%pairZBL_A,apfData(functionID)%pairZBL_B)
      xB = LJ(apfData(functionID)%functionParameters, zblE)
      xBd = LJ_Deriv(apfData(functionID)%functionParameters, zblE)
      coefficients = SplineExpThird(zblM,yArray(1),yArray(2),zblE,xB,xBd)
    End If
    Do i=1,apfData(functionID)%dataPointCount
!--------------------------
! calculate y
!--------------------------
! Make points
      If(zblE.gt.0.0D0.and.x.le.zblE)Then
        If(x.le.zblM)Then
          y = Zbl(x,apfData(functionID)%pairZBL_A,apfData(functionID)%pairZBL_B)
        End If
        If(x.gt.zblM.and.x.le.zblE)Then
          y = CalcPolynomialExp(coefficients, x)
        End If

      Else
        If(apfData(functionID)%functionForm.eq.1)Then      ! LJ
          y = LJ(apfData(functionID)%functionParameters, x)
        End If
        If(apfData(functionID)%functionForm.eq.2)Then      ! Basic Density
          y = BasicDensity(apfData(functionID)%functionParameters, x)
        End If
        If(apfData(functionID)%functionForm.eq.3)Then      ! Basic Embe
          y = BasicEmbed(apfData(functionID)%functionParameters, x)
        End If
        If(apfData(functionID)%functionForm.eq.4)Then      ! Morse
          y = LJ(apfData(functionID)%functionParameters, x)
        End If
      End If
! Store x-y
      apfData(functionID)%dataPoints(i,1)=x
      apfData(functionID)%dataPoints(i,2)=y
! Increment x
      x = x + xIncrement
    End Do
    !apfData(functionID)%functionForm,
  End Subroutine makePoints

! --------------------------------------------------------------------------------------------------!

  Subroutine updateAnalytic()
! Read in parameters
    Implicit None   ! Force declaration of all variables
! Private
    Integer(kind=StandardInteger) :: functionID
    Do functionID=1,eamFunctionCount
      Call makePoints(functionID)
    End Do
  End Subroutine updateAnalytic

! --------------------------------------------------------------------------------------------------!

  Subroutine varyAnalytic(apfData_Fin, varyAmount)
! Read in parameters
    Implicit None   ! Force declaration of all variables
! Private
    Type(analyticFunctions), Dimension(1:100) :: apfData_Fin
    Real(kind=DoubleReal) :: varyAmount, randVariation, randFloat
    Integer(kind=StandardInteger) :: functionID, paramID
! Load input apfData
    apfData = apfData_Fin
! Loop through functions
    Do functionID=1,apfData(1)%functionCount
! Vary function parameters
      Do paramID=1,apfData(functionID)%functionParameterCount
        randVariation = RandomDist('H','N',0.0D0,1.0D0)
        randVariation = varyAmount * randVariation * apfData(functionID)%functionParameterCount
        randFloat = RandomFloat(0.0D0,1.0D0)
        If(randFloat.le.0.5D0)Then
          apfData(functionID)%functionParameterCount = apfData(functionID)%functionParameterCount - randVariation
        Else
          apfData(functionID)%functionParameterCount = apfData(functionID)%functionParameterCount + randVariation
        End If
      End Do
! Make eam function points
      Call makePoints(functionID)
    End Do
!
  End Subroutine varyAnalytic

! --------------------------------------------------------------------------------------------------!

  Subroutine varyAnalyticRandom(apfData_Fin)
! Vary parameters randomly within the bounds given by the user
    Implicit None   ! Force declaration of all variables
! Private
    Type(analyticFunctions), Dimension(1:100) :: apfData_Fin
    Real(kind=DoubleReal) :: randFloat
    Real(kind=DoubleReal) :: lower, upper
    Integer(kind=StandardInteger) :: functionID, paramID
! Load input apfData
    apfData = apfData_Fin
! Loop through functions
    Do functionID=1,apfData(1)%functionCount
! Vary function parameters
      Do paramID=1,apfData(functionID)%functionParameterCount
        lower = apfData(functionID)%functionParametersLB(paramID)
        upper = apfData(functionID)%functionParametersUB(paramID)
        If(lower.eq.0.0D0.and.upper.eq.0.0D0)Then
          ! leave settings
        Else
          randFloat = RandomFloat(0.0D0,1.0D0)
          apfData(functionID)%functionParameters(paramID) = lower+randFloat*(upper-lower)
        End If
      End Do
! Make eam function points
      Call makePoints(functionID)
    End Do
!
  End Subroutine varyAnalyticRandom



! --------------------------------------------------------------------------------------------------!
!                                                                                                   !
! MODULE FUNCTIONS                                                                                  !
!                                                                                                   !
!                                                                                                   !
! --------------------------------------------------------------------------------------------------!

! -----------------------------------------------------------------------------------
  Function LJ(parameters, r) result (vr)
! Lennard Jones Potential
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: vr
! Private
    vr = 4.0D0*parameters(1)*((parameters(2)/r)**12-(parameters(2)/r)**6)
  End Function LJ
!-----------
  Function LJ_Deriv(parameters, r) result (dvr)
! Lennard Jones Potential
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: dvr
! Private
    dvr = 4.0D0*parameters(1)*&
          (12.0D0*((parameters(2)/r)**11)*((r-parameters(2))/r**2)-&
          6.0D0*((parameters(2)/r)**5)*((r-parameters(2))/r**2))
  End Function LJ_Deriv
! -----------------------------------------------------------------------------------
  Function BasicDensity(parameters, r) result (p)
! Lennard Jones Potential
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: p
! Private
    p = parameters(1)*parameters(2)*r**2*exp(-1.0D0*parameters(1)*r**2)
  End Function BasicDensity
! -----------------------------------------------------------------------------------
  Function BasicEmbed(parameters, p) result (F)
! Lennard Jones Potential
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: p
! Out
    Real(kind=DoubleReal) :: F
! Private
    F = parameters(1)*p**0.5+parameters(2)*p**2+parameters(3)*p**4
  End Function BasicEmbed
! -----------------------------------------------------------------------------------
  Function Morse(parameters, r) result (vr)
! Lennard Jones Potential
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: vr
! Private
    vr = parameters(1)*&
      ((1.0D0-exp(-1.0D0*parameters(2)*(r-parameters(3))))**2.0D0-1.0D0)
  End Function Morse
!-----------
  Function Morse_Deriv(parameters, r) result (dvr)
! Lennard Jones Potential
    Implicit None   ! Force declaration of all variables
! In
    Real(kind=DoubleReal), Dimension(:) :: parameters
    Real(kind=DoubleReal) :: r
! Out
    Real(kind=DoubleReal) :: dvr
! Private
    dvr = parameters(1)*&
      ((1.0D0-exp(-1.0D0*parameters(2)*(r-parameters(3))))**2.0D0-1.0D0)
  End Function Morse_Deriv




! --------------------------------------------------------------------------------------------------!


End Module analytic
