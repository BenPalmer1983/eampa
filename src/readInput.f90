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
    Integer(kind=StandardInteger) :: i, j
    Integer(kind=StandardInteger) :: ios
    Character(len=8) :: tempName
    Character(len=255) :: fileRow
    Character(len=32) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
    Real(kind=DoubleReal) :: timeStartRI, timeEndRI
    Real(kind=DoubleReal) :: dftEnergyTemp
    Integer(kind=StandardInteger) :: dftAtomCountTemp
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
! Step 1 - strip input file of comment lines etc
!
! Read input file, and make new stripped temp input file
    If(mpiProcessID.eq.0)Then
      Call tempFileName(tempName)
      inputFilePathT = trim(tempDirectory)//"/"//tempName//".temp.in"
    End If
    Call M_distChar(inputFilePathT)
    Call fileToClean(inputFilePathT)
! Master process only
    If(mpiProcessID.eq.0)Then
      Open(UNIT=1,FILE=trim(inputFilePath))
      Open(UNIT=2,FILE=trim(inputFilePathT))
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
          If(fileRow(1:1).eq."#")Then
            write(2,"(A)") Trim(Adjustl(StrToUpper(fileRow)))
          Else
            write(2,"(A)") Trim(Adjustl(fileRow))
          End If
        End If
      End Do
      Close(2)
      Close(1)
    End If
! Synch processes
    Call M_synchProcesses()
!
! Step 2 - Check for verbose
!
! Read file  - check for verbose
    Open(UNIT=1,FILE=trim(inputFilePathT))
    Do i=1,maxFileRows
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
      If(fileRow(1:6).eq."#PRINT")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        If(StrToUpper(fileRow(1:1)).eq."Y")Then
          printToTerminal = 1
        End If
      End If
    End Do
    Close(1)
! Verbose
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"Reading user input file ",trim(inputFilePathT)
    End If
!
! Step 3 - Read In User Settings
!
! Read in user variables
    Open(UNIT=1,FILE=trim(inputFilePathT))
    Do i=1,maxFileRows
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
! ----------------------------------
! Type of run
! ----------------------------------
      If(fileRow(1:8).eq."#RUNTYPE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        fileRow = StrToUpper(Trim(Adjustl(fileRow)))
        eampaRunType = fileRow(1:4)
        If(fileRow(1:4).eq."ENER")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionCalcEnergies = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: ENERGY"
          End If
        End If
        If(fileRow(1:4).eq."EVAL")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionEval = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: EVALUATION"
          End If
        End If
        If(fileRow(1:4).eq."EVAF")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionEvalFull = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: EVALUATION [FULL]"
          End If
        End If
        If(fileRow(1:4).eq."EVAT")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionEvalFull = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: EVALUATION [EAM Testing - BM, Alat, EC, EoS]"
          End If
        End If
        If(fileRow(1:3).eq."PWB")Then
          optionRunPWBatch = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: PWSCF BATCH INPUT FILES"
          End If
        End If
        If(fileRow(1:4).eq."OPTI")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionCalcEnergies = 0
          optionOptimise = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: OPTIMISE POTENTIAL"
          End If
        End If
        If(fileRow(1:4).eq."OPTF")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionCalcEnergies = 0
          optionOptimise = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: OPTIMISE POTENTIAL [FULL]"
          End If
        End If
        If(fileRow(1:4).eq."OPTT")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionCalcEnergies = 0
          optionOptimise = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: OPTIMISE POTENTIAL [EAM Testing - BM, Alat, EC, EoS]"
          End If
        End If
        If(fileRow(1:4).eq."OPTE")Then
          optionReadEAM = 1
          optionReadConf = 1
          optionNeighbourList = 1
          optionCalcEnergies = 0
          optionOptimise = 1
          optionOutput = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: OPTIMISE POTENTIAL [EXTENSIVE]"
          End If
        End If
        If(fileRow(1:4).eq."EAMP")Then
          optionReadEAM = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: EAM Potential Prepare"
          End If
        End If
        If(fileRow(1:4).eq."TEST")Then
          optionReadEAM = 1
          optionMakeConf = 1
          optionTestEAM = 1
          If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
            print *,"Run Type: Test EAM Potential"
          End If
        End If
      End If
! ----------------------------------
! MPI Options
! ----------------------------------
      If(fileRow(1:10).eq."#MPIENERGY")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) mpiEnergy
      End If
! ----------------------------------
! EAM Potential
! ----------------------------------
      If(fileRow(1:10).eq."#POTENTIAL")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        eamFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:12).eq."#EAMPREPFILE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        eamSaveFile = trim(adjustl(fileRow))
      End If
      If(fileRow(1:16).eq."#EAMINTERPPOINTS")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) eamInterpPoints
      End If
      If(fileRow(1:8).eq."#ZBLCORE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
        Read(bufferA,*) zblHardCore(1)
        Read(bufferB,*) zblHardCore(2)
        Read(bufferC,*) zblHardCore(3)
        Read(bufferD,*) zblHardCore(4)
        Read(bufferE,*) zblHardCore(5)
        Read(bufferF,*) zblHardCore(6)
      End If
      If(fileRow(1:12).eq."#SPLINENODES")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Call strToIntArr(fileRow,splineNodeCount)
        splineTotalNodes = 0
      End If
      If(fileRow(1:13).eq."#EAMNODESFILE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        eamNodesFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:15).eq."#EAMFORCESPLINE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) eamForceSpline
      End If
      If(fileRow(1:12).eq."#EAMFORCEZBL")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) eamForceZBL
      End If
      If(fileRow(1:13).eq."#EAMMAKEALLOY")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Call strToStrArr(fileRow,eamMakeAlloy)
      End If
      If(fileRow(1:13).eq."#EAMFILETYPE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) eamFileType
      End If
! ----------------------------------
! Config Details
! ----------------------------------
      If(fileRow(1:15).eq."#CONFIGURATIONS")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        configFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:11).eq."#UNITVECTOR")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) globalConfigUnitVector(1,1)
        Read(bufferB,*) globalConfigUnitVector(1,2)
        Read(bufferC,*) globalConfigUnitVector(1,3)
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) globalConfigUnitVector(2,1)
        Read(bufferB,*) globalConfigUnitVector(2,2)
        Read(bufferC,*) globalConfigUnitVector(2,3)
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) globalConfigUnitVector(3,1)
        Read(bufferB,*) globalConfigUnitVector(3,2)
        Read(bufferC,*) globalConfigUnitVector(3,3)
      End If
      If(fileRow(1:15).eq."#SAVECONFIGFILE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        saveConfigFile = trim(adjustl(fileRow))
      End If
      If(fileRow(1:18).eq."#SAVEEXPCONFIGFILE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        saveExpConfigFile = trim(adjustl(StrToUpper(fileRow)))
      End If
! ----------------------------------
! DFT Settings
! ----------------------------------
      If(fileRow(1:15).eq."#OPTENSTART")Then
        Do j=1,300
          Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
          fileRow = Trim(StrToUpper(fileRow))
          If(fileRow(1:1).eq."#")Then
            Exit
          Else
            Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
            Read(bufferB,*) dftEnergyTemp
            Read(bufferD,*) dftAtomCountTemp
            dftElement(j) = bufferA(1:2)
            dftOptEnergy(j) = UnitConvert(dftEnergyTemp, bufferC, "EV")
            dftOptEnergy(j) = dftOptEnergy(j)/(1.0D0*dftAtomCountTemp)
            Read(bufferE,*) dftEnergyTemp
            dftCohEnergy(j) = UnitConvert(dftEnergyTemp, bufferF, "EV")
          End If
        End Do
      End If
! ----------------------------------
! Neighbour List Settings
! ----------------------------------
      If(fileRow(1:9).eq."#NLCUTOFF")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        fileRow = trim(adjustl(StrToUpper(fileRow)))
        Read(fileRow,*) nlCutoff
      End If
      If(fileRow(1:13).eq."#NLTESTCUTOFF")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        fileRow = trim(adjustl(StrToUpper(fileRow)))
        Read(fileRow,*) nlTestCutoff
      End If
! ----------------------------------
! Calculation options
! ----------------------------------
      If(fileRow(1:10).eq."#CALCEQVOL")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        calcEqVol = trim(adjustl(StrToUpper(fileRow)))
      End If
      If(fileRow(1:12).eq."#REFINEEQVOL")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        fileRow = trim(adjustl(fileRow))
        If(fileRow(1:1).eq."Y")Then
          refineEqVol = "YES"
        End If
      End If
      If(fileRow(1:11).eq."#SAVEFORCES")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) saveForcesToFile
      End If
      If(fileRow(1:7).eq."#SAVENL")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) saveNLToFile
      End If
! ----------------------------------
! Optimise options
! ----------------------------------
      If(fileRow(1:10).eq."#VARYNODES")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Call strToDPArr(fileRow,varyNodeOptions)
      End If
      If(fileRow(1:9).eq."#OPTLOOPS")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) optLoops
      End If
      If(fileRow(1:9).eq."#OPTLOOPS")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) optLoops
      End If
      If(fileRow(1:7).eq."#SAOPTS")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
        Read(bufferA,*) saTemp
        Read(bufferB,*) saTempLoops
        Read(bufferC,*) saVarLoops
        Read(bufferD,*) saMaxVariation
      End If
      If(fileRow(1:15).eq."#VARYFIXEDNODES")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) varyFixedNodes
      End If
      If(fileRow(1:12).eq."#JUMBLENODES")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) jumbleNodesOpt
      End If
      If(fileRow(1:12).eq."#REDUCENODES")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) reduceNodes
      End If
      If(fileRow(1:12).eq."#EMBERESCALE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) embeRescale
      End If
      If(fileRow(1:13).eq."#FORCEEMBEFIT")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) forceEmbeFitOpt
      End If
! ----------------------------------
! RSS calculation options
! ----------------------------------
      If(fileRow(1:13).eq."#RSSWEIGHTING")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Call strToDPArr(fileRow,rssWeighting)
      End If
! ----------------------------------
! Testing Options
! ----------------------------------
      If(fileRow(1:7).eq."#EOSFIT")Then       !DELETE
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) testingFitChoice
      End If
      If(fileRow(1:10).eq."#EOSRSSFIT")Then   !DELETE
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) eosFitRSSOption
      End If
! ----------------------------------
! Testing Reference Data
! ----------------------------------
      If(fileRow(1:7).eq."#FCCREF")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Call strToDPArr(fileRow,fccReferenceValues)
      End If
! ----------------------------------
! PWscf Batch Files
! ----------------------------------
      If(fileRow(1:10).eq."#PWBTYPE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        fileRow = StrToUpper(fileRow)
        pwbRunType = trim(adjustl(fileRow))
      End If
      If(fileRow(1:10).eq."#PWBCONFIG")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        pwbConfigFilePath = trim(adjustl(fileRow))
      End If
      If(fileRow(1:12).eq."#PWBBATCHDIR")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        pwbBatchDir = trim(adjustl(fileRow))
      End If
      If(fileRow(1:18).eq."#PWBVARIANCESWITCH")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) pwbVarianceSwitch
      End If
      If(fileRow(1:18).eq."#PWBRANDOMVARIANCE")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) pwbVarianceType
        Read(bufferB,*) pwbVarianceMax
        Read(bufferC,*) pwbVarianceSigma
      End If
      If(fileRow(1:16).eq."#PWBINTERSTITIAL")Then
        Read(1,"(A255)",IOSTAT=ios) fileRow   !read next line
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferA,*) pwbInterstitialAtom
        Read(bufferB,*) pwbInterstitialDetails(1)
        Read(bufferC,*) pwbInterstitialDetails(2)
      End If
    End Do
    Close(1)
! Synch MPI processes
    Call M_synchProcesses()
! End Time
    Call cpu_time(timeEndRI)
! Store Time
    Call storeTime(8,timeEndRI-timeStartRI)
  End Subroutine readUserInput
End Module readinput
