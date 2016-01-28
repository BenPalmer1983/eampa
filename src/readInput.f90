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
    Integer(kind=StandardInteger) :: dftAtomCountTemp, fileRows
! Read in command line arguments
    Call get_command_argument(1,inputFilePath)    
! Step 1 - read to memory and strip input file of comment lines etc
! Read input file into memory    
    Call readFile(trim(inputFilePath), userInputData, fileRows)
! Step 2 - Check for verbose and # to caps
    Do i=1,fileRows
      fileRow = userInputData(i)
      If(fileRow(1:1).eq."#")Then
        userInputData(i) = StrToUpper(userInputData(i))
        fileRow = userInputData(i)
      End If
      If(fileRow(1:6).eq."#PRINT")Then
        fileRow = userInputData(i+1)   !read next line
        printToTerminal = UserTrueFalse(fileRow)
      End If  
    End Do
!
! Step 3 - Read In User Settings
!
! ------------------
! Set some default values
    eampaRunType = "EVAL"
! Read in user variables
    j = 1
    Do i=1,fileRows
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
      If(fileRow(1:8).eq."#POTTYPE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        fileRow = StrToUpper(fileRow)
        potentialType = 1 ! assume tabulated
        If(fileRow(1:1).eq."2".or.fileRow(1:8).eq."ANALYTIC")Then
          potentialType = 2
        End If      
      End If
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
      !If(fileRow(1:12).eq."#OPTFORCEZBL")Then
      !  j = j + 1
      !  fileRow = userInputData(j)   !read next line
      !  optForceZBL = UserTrueFalse(fileRow)
      !End If
      If(fileRow(1:15).eq."#ANOPTRANDLOOPS")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA
        Read(bufferA,*) anOptRandLoops
      End If
      
      
      If(fileRow(1:7).eq."#SAOPTS")Then
        Do k = 1,20
! Check next row
          fileRow = userInputData(j+1)
          If(fileRow(1:1).eq."#")Then
            fileRow = userInputData(j)
            Exit
          Else 
            j = j + 1
            Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
            bufferA = StrToUpper(bufferA)
            saConfigIn(k)%saType = 1
            If(bufferA(1:1).eq."2".or.bufferA(1:1).eq."A")Then
              saConfigIn(k)%saType = 2           
            End If        
            Read(bufferB,*) saConfigIn(k)%tempStart
            Read(bufferC,*) saConfigIn(k)%tempEnd
            Read(bufferD,*) saConfigIn(k)%varLoops
            Read(bufferE,*) saConfigIn(k)%maxVar  
            Read(bufferF,*) saConfigIn(k)%minVar 
            saConfigIn(1)%saCount = k 
          End If
        End Do  
      End If
      If(fileRow(1:11).eq."#OPTRUNTYPE")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        Read(fileRow,*) bufferA
        Read(bufferA,*) optRunType
      End If
      If(fileRow(1:11).eq."#OPTFROM")Then
        j = j + 1
        fileRow = userInputData(j)   !read next line
        optFrom = 1 ! Assume from potential
        Read(fileRow,*) bufferA
        bufferA = StrToUpper(bufferA)
        If(bufferA(1:1).eq."2")Then
          optFrom = 2
        End If       
        If(bufferA(1:7).eq."SCRATCH")Then
          optFrom = 2
        End If  
      End If
      
      
      
      !If(fileRow(1:15).eq."#VARYFIXEDNODES")Then
      !  j = j + 1
      !  fileRow = userInputData(j)   !read next line
      !  Read(fileRow,*) varyFixedNodes
      !End If
      !If(fileRow(1:12).eq."#JUMBLENODES")Then
      !  j = j + 1
      !  fileRow = userInputData(j)   !read next line
      !  Read(fileRow,*) jumbleNodesOpt
      !End If
      !If(fileRow(1:12).eq."#REDUCENODES")Then
      !  j = j + 1
      !  fileRow = userInputData(j)   !read next line
      !  Read(fileRow,*) reduceNodes
      !End If
      !If(fileRow(1:12).eq."#EMBERESCALE")Then
      !  j = j + 1
      !  fileRow = userInputData(j)   !read next line
      !  Read(fileRow,*) embeRescale
      !End If
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
