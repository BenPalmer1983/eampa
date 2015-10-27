Module readinput
! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Read user input file
! ----------------------------------------
! Updated: 25th Aug 2014
! ----------------------------------------
! Setup Modules
  Use kinds
  Use msubs
  Use constants
  Use maths
  Use general
  Use units
  Use initialise
  Use loadData
  Use globals
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: readUserInput
  Contains
! ---------------------------------------------------------------------------------------------------
  Subroutine readUserInput()
! force declaration of all variables
    Implicit None
! Private variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger) :: ios
    Character(len=255) :: fileRow
    Character(len=32) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
    Real(kind=DoubleReal) :: timeStartRI, timeEndRI
    Real(kind=DoubleReal) :: dftEnergyTemp
    Integer(kind=StandardInteger) :: dftAtomCountTemp
! Print out
    If(mpiProcessID.eq.0)Then
      print *,"READ INPUT:"
    End If
! Start Time
    Call cpu_time(timeStartRI)
! Read in command line arguments
    Call get_command_argument(1,inputFilePath)
! Write to output.dat file
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(F8.4,A2,A24,A60)") ProgramTime(),"  ",&
      "Reading user input file ",inputFilePath
    End If
! Step 1 - read to memory and strip input file of comment lines etc
! Read input file into memory
    j = 0
    Open(UNIT=1,FILE=trim(inputFilePath))
    Do i=1,maxFileRows
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
      fileRow = Trim(Adjustl(fileRow))
      fileRow = RemoveComments(fileRow)
      fileRow = RemoveQuotes(fileRow)
      If(fileRow(1:1).eq."!".or.fileRow(1:1).eq." ")Then
! skip
      Else
        j = j + 1
        If(fileRow(1:1).eq."#")Then
          userInputData(j) = Trim(Adjustl(StrToUpper(fileRow)))
        Else
          userInputData(j) = Trim(Adjustl(fileRow))
        End If
! if(mpiProcessID.eq.0)Then
!  print *,j,trim(userInputData(j))
! End If
      End If
    End Do
    Close(1)
!
! Step 2 - Check for verbose
!
! Read file  - check for verbose
    j = 1
    Do i=1,maxFileRows
      fileRow = userInputData(j)
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
      If(fileRow(1:6).eq."#PRINT")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        printToTerminal = UserTrueFalse(fileRow)
      End If
      j = j + 1
    End Do
! Verbose
    If(TerminalPrint())Then
      Print *,"Reading user input file ",trim(inputFilePath)
    End If
!
! Step 3 - Read In User Settings
!
! ------------------
! Set some default values
    eampaRunType = "EVAL"
! Read in user variables
    j = 1
    Do i=1,maxFileRows
      fileRow = userInputData(j)
      If(fileRow(1:1).eq." ")Then
        EXIT
      End If
! ----------------------------------
! Type of run
! ----------------------------------
      If(fileRow(1:8).eq."#RUNTYPE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        fileRow = StrToUpper(Trim(Adjustl(fileRow)))
! Only allow certain run types
        If(fileRow(1:4).eq."EVAL")Then
          eampaRunType = fileRow(1:4)  ! EVAL
        End If
        If(fileRow(1:4).eq."OPTI")Then
          eampaRunType = fileRow(1:4)  ! OPTI
        End If
      End If
! ----------------------------------
! EAM Potential
! ----------------------------------
      If(fileRow(1:10).eq."#POTENTIAL")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        eamFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:13).eq."#EAMFILETYPE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) eamFileType
      End If
      If(fileRow(1:12).eq."#EAMPREPFILE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        eamSaveFile = trim(adjustl(fileRow))
      End If
      If(fileRow(1:16).eq."#EAMINTERPPOINTS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) eamInterpPoints
      End If
      If(fileRow(1:8).eq."#ZBLCORE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA
        Read(bufferA,*) zblHardCore(1)
      End If
      If(fileRow(1:12).eq."#SPLINENODES")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Call strToIntArr(fileRow,splineNodeCount)
        splineTotalNodes = 0
      End If
      If(fileRow(1:13).eq."#EAMNODESFILE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        eamNodesFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:15).eq."#EAMFORCESPLINE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        eamForceSpline = UserTrueFalse(fileRow)
      End If
      If(fileRow(1:13).eq."#EAMMAKEALLOY")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Call strToStrArr(fileRow,eamMakeAlloy)
      End If
      If(fileRow(1:10).eq."#EAMCHARTS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        makeEAMCharts = UserTrueFalse(fileRow)
      End If
! ----------------------------------
! Config Details
! ----------------------------------
      If(fileRow(1:15).eq."#CONFIGURATIONS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        configFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:11).eq."#UNITVECTOR")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) globalConfigUnitVector(1,1)
        Read(bufferB,*) globalConfigUnitVector(1,2)
        Read(bufferC,*) globalConfigUnitVector(1,3)
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) globalConfigUnitVector(2,1)
        Read(bufferB,*) globalConfigUnitVector(2,2)
        Read(bufferC,*) globalConfigUnitVector(2,3)
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) globalConfigUnitVector(3,1)
        Read(bufferB,*) globalConfigUnitVector(3,2)
        Read(bufferC,*) globalConfigUnitVector(3,3)
      End If
      If(fileRow(1:15).eq."#SAVECONFIGFILE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        saveConfigFile = trim(adjustl(fileRow))
      End If
      If(fileRow(1:18).eq."#SAVEEXPCONFIGFILE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        saveExpConfigFile = trim(adjustl(StrToUpper(fileRow)))
      End If
! ----------------------------------
! DFT Settings
! ----------------------------------
      If(fileRow(1:15).eq."#OPTENSTART")Then
        Do k=1,300
          j = j + 1
          fileRow = userInputData(j)   !read next line
          fileRow = Trim(StrToUpper(fileRow))
          If(fileRow(1:1).eq."#")Then
            Exit
          Else
            Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
            Read(bufferB,*) dftEnergyTemp
            Read(bufferD,*) dftAtomCountTemp
            dftElement(k) = StrToUpper(bufferA(1:2))
            dftOptEnergy(k) = UnitConvert(dftEnergyTemp, bufferC, "EV")
            dftOptEnergy(k) = dftOptEnergy(k)/(1.0D0*dftAtomCountTemp)
            Read(bufferE,*) dftEnergyTemp
            dftCohEnergy(k) = UnitConvert(dftEnergyTemp, bufferF, "EV")
          End If
        End Do
      End If
! ----------------------------------
! Bulk Property Config Details
! ----------------------------------
      If(fileRow(1:17).eq."#BPCONFIGURATIONS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        bpConfigFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:12).eq."#BPPRINTDATA")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        bpPrintData = UserTrueFalse(fileRow)
      End If
! ----------------------------------
! Neighbour List Settings
! ----------------------------------
      If(fileRow(1:9).eq."#NLCUTOFF")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        fileRow = trim(adjustl(StrToUpper(fileRow)))
        Read(fileRow,*) nlCutoff
      End If
      If(fileRow(1:13).eq."#NLTESTCUTOFF")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        fileRow = trim(adjustl(StrToUpper(fileRow)))
        Read(fileRow,*) nlTestCutoff
      End If
! ----------------------------------
! Calculation options
! ----------------------------------
      If(fileRow(1:10).eq."#CALCEQVOL")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        calcEqVol = trim(adjustl(StrToUpper(fileRow)))
      End If
      If(fileRow(1:12).eq."#REFINEEQVOL")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        fileRow = trim(adjustl(fileRow))
        If(fileRow(1:1).eq."Y")Then
          refineEqVol = "YES"
        End If
      End If
      If(fileRow(1:11).eq."#SAVEFORCES")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        saveForcesToFile = UserTrueFalse(fileRow)
      End If
      If(fileRow(1:7).eq."#SAVENL")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        saveNLToFile = UserTrueFalse(fileRow)
      End If
      If(fileRow(1:9).eq."#EOSCHART")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        eosChart = UserTrueFalse(fileRow)
      End If
! ----------------------------------
! Optimise options
! ----------------------------------
      If(fileRow(1:12).eq."#OPTFORCEZBL")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        optForceZBL = UserTrueFalse(fileRow)
      End If
      If(fileRow(1:10).eq."#VARYNODES")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Call strToDPArr(fileRow,varyNodeOptions)
      End If
      If(fileRow(1:9).eq."#OPTLOOPS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) optLoops
      End If
      If(fileRow(1:9).eq."#OPTLOOPS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) optLoops
      End If
      If(fileRow(1:7).eq."#SAOPTS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
        Read(bufferA,*) saConfigIn%temp
        Read(bufferB,*) saConfigIn%tempEnd
        Read(bufferC,*) saConfigIn%varLoops
        Read(bufferD,*) saConfigIn%maxVar  
        Read(bufferE,*) saConfigIn%minVar 
      End If
      If(fileRow(1:15).eq."#VARYFIXEDNODES")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) varyFixedNodes
      End If
      If(fileRow(1:12).eq."#JUMBLENODES")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) jumbleNodesOpt
      End If
      If(fileRow(1:12).eq."#REDUCENODES")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) reduceNodes
      End If
      If(fileRow(1:12).eq."#EMBERESCALE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) embeRescale
      End If
      If(fileRow(1:13).eq."#FORCEEMBEFIT")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) forceEmbeFitOpt
      End If
! ----------------------------------
! RSS calculation options
! ----------------------------------
      If(fileRow(1:13).eq."#RSSWEIGHTING")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Call strToDPArr(fileRow,rssWeighting)
      End If
! ----------------------------------
! Testing Options
! ----------------------------------
      If(fileRow(1:7).eq."#EOSFIT")Then       !DELETE
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) testingFitChoice
      End If
      If(fileRow(1:10).eq."#EOSRSSFIT")Then   !DELETE
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) eosFitRSSOption
      End If
      j = j + 1
    End Do
! Synch MPI processes
    Call M_synchProcesses()
! End Time
    Call cpu_time(timeEndRI)
! Store Time
    Call storeTime(8,timeEndRI-timeStartRI)
! Output
    If(TerminalPrint())Then
      print *,"User input read: ",(timeEndRI-timeStartRI),"s"
    End If
  End Subroutine readUserInput

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!

  Function UserTrueFalse (inputText) RESULT (outputBool)
! Takes user input and figures out if they wanted true or false
    Character(*) :: inputText
    Character(len=6) :: userChoice
    Logical :: outputBool
! Init string
    userChoice = "      "
    userChoice = StrToUpper(Trim(Adjustl(inputText)))
    outputBool = .false.
    If(userChoice(1:1).eq."1")Then
      outputBool = .true.
    End If
    If(userChoice(1:1).eq."Y")Then
      outputBool = .true.
    End If
    If(userChoice(1:1).eq."+")Then
      outputBool = .true.
    End If
    If(userChoice(1:4).eq."TRUE")Then
      outputBool = .true.
    End If
    If(userChoice(1:6).eq.".TRUE.")Then
      outputBool = .true.
    End If
  End Function UserTrueFalse
End Module readinput
