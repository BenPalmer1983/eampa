Module globals

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Declare all global variables

! ----------------------------------------
! Updated: 18th June 2015
! ----------------------------------------

! Setup Modules
  Use libBP
  Use types
  Use msubs
  
! Force declaration of all variables
  Implicit None
!------------------------------------------------------------------------------ 
! Global Parameters
!
  Integer(kind=StandardInteger), Parameter :: atomsPerConfig = 1024
  Integer(kind=StandardInteger), Parameter :: nlLenPerConfig = 16000
  Integer(kind=StandardInteger), Parameter :: maxConfigs = 64    ! Maximum number of configurations able to load
  Integer(kind=StandardInteger), Parameter :: maxConfigsBP = 32    ! Maximum number of bulk-property testing configurations
  Integer(kind=StandardInteger), Parameter :: maxConfigsRelax = 4    ! Maximum number of bulk-property testing configurations
  Integer(kind=StandardInteger), Parameter :: bpCopies = 4         ! Size (nxnxn) of bulk property crystals
  Integer(kind=StandardInteger), Parameter :: aLatSamples = 7      ! number of energy-volume samples used to plot equation of state
  Integer(kind=StandardInteger), Parameter :: ecSamples = 4        ! number of samples used to calculate cubic elastic constants
  Real(kind=DoubleReal), Parameter :: fccAlatBP = 4.0D0
  Real(kind=DoubleReal), Parameter :: bccAlatBP = 2.6D0
  Real(kind=DoubleReal), Parameter :: bpCutoffNL = 6.5D0
  Real(kind=DoubleReal), Parameter :: bpCutoff = 6.5D0
  Character(len=4), Parameter, Dimension(1:10) :: eamFunctionTypes = &
  ["PAIR","DENS","EMBE","DDEN","SDEN","DEMB","SEMB","    ","    ","    "]
!     1      2      3      4      5      6      7      8      9     10  
!------------------------------------------------------------------------------ 

!------------------------------------------------------------------------------ 
! Program Configuration/Settings
!
  Type(progSettings) :: eampaSettings


!------------------------------------------------------------------------------ 
! My standard variables
!
  Character(len=64) :: compileLine
  Real(kind=DoubleReal) :: programStartTime, programEndTime
  Real(kind=DoubleReal), Dimension(1:100) :: cpuTime
  Character(Len=64), Dimension(1:100) :: cpuTimeLabels
  Character(len=255) :: currentWorkingDirectory
  Character(len=512), Dimension(1:100) :: fileCleanupList
  Character(len=255) :: outputDirectory
  Character(len=255) :: tempDirectory
! MPI Global Variables
  Integer(kind=StandardInteger) :: mpiProcessCount, mpiProcessID
! System Variables
  Real(kind=DoubleReal) :: largeArraySize
  Integer(kind=StandardInteger), Dimension(1:maxConfigs) :: processMap
  Integer(kind=StandardInteger), Dimension(1:maxConfigs) :: processMapBP
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:10) :: processMapE
! Other useful flags
  Logical :: quietOverride  

  
!------------------------------------------------------------------------------ 
! Initialise
! 
  Character(len=255) :: outputFile
  Character(len=255) :: outputFileEnergies
  Character(len=255) :: outputFileForces
  
  
!------------------------------------------------------------------------------ 
! Load isotope data
!   
  
  
  
!------------------------------------------------------------------------------ 
! Read Input
!   
! user file array
  Character(len=255), Dimension(1:1024) :: userInputData
! Verbose
  !Integer(kind=StandardInteger) :: printToTerminal
  Logical :: printToTerminal
  Character(len=4) :: eampaRunType
! Input File - User Input  
  Character(len=255) :: inputFilePath
! EAM Details - User Input
  Character(len=255) :: eamFilePath
  Character(len=255) :: eamNodesFilePath
  Character(len=64) :: eamSaveFile
  Integer(kind=StandardInteger) :: potentialType = 1 ! 1 tabulated  2 analytic
  Integer(kind=StandardInteger) :: eamInterpPoints
  Real(kind=DoubleReal), Dimension(1:6) :: zblHardCore = 1.0D0                            ! 1 Pair ZBL end, 2 Pair Spline End, 3 Dens Value, 4 De
  Integer(kind=StandardInteger), Dimension(1:50) :: splineNodeCount
  Integer(kind=StandardInteger) :: splineTotalNodes
  Character(len=2), Dimension(1:10) :: eamMakeAlloy
  Integer(kind=StandardInteger) :: eamFileType
  Logical :: makeEAMCharts
! Config Details - User Input
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: globalConfigUnitVector
  Character(len=255) :: configFilePath
  Character(len=255) :: configFilePathT
 ! Character(len=255) :: saveConfigFile
  Character(len=255) :: saveExpConfigFile
! BP Config Details - User Input
  Character(len=255) :: bpConfigFilePath
  Logical :: bpPrintData = .false.
! DFT Settings
  Character(len=2), Dimension(1:300) :: dftElement
  Real(kind=DoubleReal), Dimension(1:300) :: dftOptEnergy
  Real(kind=DoubleReal), Dimension(1:300) :: dftCohEnergy
! Neighbour List Settings
  Real(kind=DoubleReal) :: nlCutoff       ! Standard calculations
  Real(kind=DoubleReal) :: nlTestCutoff   ! Test calculations
! Calculation details
  Character(len=8) :: calcEqVol
  Character(len=3) :: refineEqVol
! Input File Logical
  Logical :: saveForcesToFile
  Logical :: saveNLToFile
  Logical :: eamForceSpline
  Logical :: eamForceZBL
  Logical :: eosChart
! Optimise options
  Integer(kind=StandardInteger) :: anOptRandLoops=0
  Real(kind=DoubleReal), Dimension(1:10) :: varyNodeOptions
  Integer(kind=StandardInteger) :: optRunType=0
  Integer(kind=StandardInteger) :: optFrom=1
!  
  Integer(kind=StandardInteger) :: optLoops
  Integer(kind=StandardInteger) :: reduceNodes
  Integer(kind=StandardInteger) :: embeRescale
! RSS calculation options
  Real(kind=DoubleReal), Dimension(1:20) :: rssWeighting
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configWeighting
! PW Batch Files - User Input
  Character(len=16) :: pwbRunType
  Character(len=255) :: pwbConfigFilePath                                          ! 255Bytes
  Character(len=255) :: pwbConfigFilePathT                                         ! 255Bytes
  Character(len=255) :: pwbBatchDir
  Integer(kind=StandardInteger) :: pwbVarianceSwitch
  Character(len=4) :: pwbVarianceType
  Real(kind=DoubleReal) :: pwbVarianceMax
  Real(kind=DoubleReal) :: pwbVarianceSigma
  Integer(kind=StandardInteger) :: pwbInterstitialAtom
  Character(len=16), Dimension(1:3) :: pwbInterstitialDetails  
  Type(bulkProperties), Dimension(1:maxConfigsBP) :: refBulkProperties
    
!------------------------------------------------------------------------------ 
! Read EAM File
!   
! EAM file array  
  Character(len=255), Dimension(1:65536) :: eamInputData
  Character(len=2), Dimension(1:300) :: elements                                   ! 0.6KB
  Integer(kind=StandardInteger) :: elementsCount
  Integer(kind=StandardInteger), Dimension(1:300) :: elementsCharge                ! 1.2KB
! Analytic potential functions  
  Type(analyticFunctions), Dimension(1:100) :: apfData
  Type(analyticFunctions), Dimension(1:100) :: apfDataIn
  Type(analyticFunctions), Dimension(1:100) :: apfDataOpt
! Read EAM File
  Integer(kind=StandardInteger) :: eamFunctionCount, eamPairCount, eamDensCount, eamEmbeCount
  Integer(kind=StandardInteger) :: eamDdenCount, eamSdenCount, eamDembCount, eamSembCount
  Integer(kind=StandardInteger) :: eamType                                         ! 4bit        1=EAM, 2=2BMEAM
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: eamKey                     ! 1.2KB       1 atomA, 2 atomB, 3 function/al type, 4 func start, 5 func length, 6 func end
  Real(kind=DoubleReal), Dimension(1:100000,1:4) :: eamData                        ! 3.2MB       1 x, 2 y(x), 3 y'(x), 4 y''(x)
! Spline arrays
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: splineNodesKey, splineNodesKeyTemp
  Real(kind=DoubleReal), Dimension(1:10000,1:6) :: splineNodesData, splineNodesDataTemp
  Real(kind=DoubleReal), Dimension(1:10000,1:2) :: splineNodesResponse
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: eamKeyInput
  Real(kind=DoubleReal), Dimension(1:100000,1:4) :: eamDataInput
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: eamKeyOpt
  Real(kind=DoubleReal), Dimension(1:100000,1:4) :: eamDataOpt
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: splineNodesKeyOpt, splineNodesKeyBest
  Real(kind=DoubleReal), Dimension(1:10000,1:6) :: splineNodesDataOpt, splineNodesDataBest
    
!------------------------------------------------------------------------------ 
! Read Config File
!   
! Arrays to read config file and dft files into
  Character(len=255), Dimension(1:65536) :: configInputData
  Character(len=255), Dimension(1:65536) :: configInputDataDFT
  Character(len=255), Dimension(1:65536) :: configInputDataTemp
  Character(len=255), Dimension(1:65536) :: configInputDataDFTTemp
! atom labels to replace
  Character(len=16), Dimension(1:maxConfigs,1:2) :: configLabelReplace
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:9) :: crystalUnitCell
! Counter  
  Integer(kind=StandardInteger) :: configCount
! Input coords and forces
  Integer(kind=StandardInteger) :: coordCount                                      ! 4bit
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:3) :: configurationCoordsKey   ! 13KB        1 start, 2 length, 3 end
  Integer(kind=StandardInteger), Dimension(1:50000,1:1) :: configurationCoordsI   !             1 atomID
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationCoordsR           !             1 x, 2 y, 3 z
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationForcesR     
! Generated/expanded coords and forces
  Integer(kind=StandardInteger) :: coordCountG                                     ! 4bit
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:3) :: configurationCoordsKeyG  ! 13KB        1 start, 2 length, 3 end
  Integer(kind=StandardInteger), Dimension(1:100000,1:1) :: configurationCoordsIG  !          1 atomID
  Real(kind=DoubleReal), Dimension(1:100000,1:6) :: configurationCoordsRG          !       1 x, 2 y, 3 z
! Volumes, total atoms 
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configVolume
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configVolumeOpt      !             1 fx, 2 fy, 3 fz
  Integer(kind=StandardInteger) :: configsAtomTotal
    
!------------------------------------------------------------------------------ 
! Read BP Config File
!  
  Character(len=255), Dimension(1:65536) :: bpConfigInputData
  Type(bpIn), Dimension(1:maxConfigsBP) :: bpInArr 

 
!------------------------------------------------------------------------------ 
! Bulk Properties Config
!   
! Counter  
  Integer(kind=StandardInteger) :: configCountBP
  Integer(kind=StandardInteger), Dimension(1:maxConfigsBP,1:3) :: configurationCoordsKeyBP  !     1 start, 2 length, 3 end
  Integer(kind=StandardInteger), Dimension(1:8192,1:1) :: configurationCoordsIBP  !     1 atomID
  Real(kind=DoubleReal), Dimension(1:8192,1:6) :: configurationCoordsRBP          !     1 x, 2 y, 3 z
  Integer(kind=StandardInteger) :: configsAtomTotalBP
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP) :: aLatBP
  Integer(kind=StandardInteger), Dimension(1:maxConfigsBP) :: bpUnitCellCount
  Type(bulkProperties), Dimension(1:maxConfigsBP) :: calcBulkProperties
  
 
!------------------------------------------------------------------------------ 
! Relax Config
!   
! Counter  
  Integer(kind=StandardInteger) :: configCountRelax
  Integer(kind=StandardInteger), Dimension(1:maxConfigsRelax,1:3) :: configurationCoordsKeyRelax  !     1 start, 2 length, 3 end
  Integer(kind=StandardInteger), Dimension(1:8192,1:1) :: configurationCoordsIRelax  !     1 atomID
  Real(kind=DoubleReal), Dimension(1:8192,1:6) :: configurationCoordsRRelax          !     1 x, 2 y, 3 z
  Integer(kind=StandardInteger) :: configsAtomTotalRelax
  Real(kind=DoubleReal), Dimension(1:maxConfigsRelax) :: aLatRelax
  Integer(kind=StandardInteger), Dimension(1:maxConfigsRelax) :: relaxUnitCellCount

!------------------------------------------------------------------------------ 
! Neighbour List up to 1024 configs, nl per config approx 10,000
! Current settings limit to approx 80, increase before compiling if neccessary
! 
  Integer(kind=StandardInteger), Dimension(1:nlLenPerConfig*maxConfigs) :: nlUniqueKeys               ! 2.0MB
  Integer(kind=StandardInteger) :: neighbourListCount
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:3) :: neighbourListKey
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:6) :: neighbourListKeyR
  Integer(kind=StandardInteger), Dimension(1:nlLenPerConfig*maxConfigs,1:6) :: neighbourListI         ! 12.0MB
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigs) :: neighbourListR                     ! 4.0MB
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigs,1:12) :: neighbourListCoords             ! 48.0MB
  Integer(kind=StandardInteger), Dimension(1:2000) :: atomSeparationSpread
! Temporary NL arrays
 ! Integer(kind=StandardInteger), Dimension(1:800000,1:6) :: neighbourListIT         ! 16.0MB
 ! Real(kind=DoubleReal), Dimension(1:800000) :: neighbourListRT                     ! 5.0MB
 ! Real(kind=DoubleReal), Dimension(1:800000,12) :: neighbourListCoordsT             ! 60.0MB


!------------------------------------------------------------------------------ 
! Neighbour List for Bulk Property Configs
! Maximum usually 32
! 
  Integer(kind=StandardInteger), Dimension(1:nlLenPerConfig*maxConfigsBP) :: nlUniqueKeysBP              ! 2.0MB
  Integer(kind=StandardInteger) :: neighbourListCountBP
  Integer(kind=StandardInteger), Dimension(1:maxConfigsBP,1:3) :: neighbourListKeyBP
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:6) :: neighbourListKeyRBP
  Integer(kind=StandardInteger), Dimension(1:nlLenPerConfig*maxConfigsBP,1:6) :: neighbourListIBP      ! 12.0MB
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsBP) :: neighbourListRBP               ! 4.0MB
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsBP,1:12) :: neighbourListCoordsBP             ! 48.0MB
  Integer(kind=StandardInteger), Dimension(1:2000) :: atomSeparationSpreadBP
! Temporary NL arrays
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsBP) :: neighbourListRBP_T  
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsBP,1:6) :: neighbourListCoordsBP_T  



!------------------------------------------------------------------------------ 
! Neighbour List for geometric relaxation
! Maximum usually 4 (probably 1 config at a time)
! 
  Integer(kind=StandardInteger), Dimension(1:nlLenPerConfig*maxConfigsRelax) :: nlUniqueKeysRelax    
  Integer(kind=StandardInteger) :: neighbourListCountRelax
  Integer(kind=StandardInteger), Dimension(1:maxConfigsRelax,1:3) :: neighbourListKeyRelax
  Real(kind=DoubleReal), Dimension(1:maxConfigsRelax,1:6) :: neighbourListKeyRRelax                            ! 1 rcutoff, 2 rmin, 3 rmax, 6 alat sim cube 
  Integer(kind=StandardInteger), Dimension(1:nlLenPerConfig*maxConfigsRelax,1:11) :: neighbourListIRelax      ! 12.0MB
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsRelax) :: neighbourListRRelax               ! 4.0MB
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsRelax,1:12) :: neighbourListCoordsRelax             ! 48.0MB
  Integer(kind=StandardInteger), Dimension(1:2000) :: atomSeparationSpreadRelax
! Temporary NL arrays
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsRelax) :: neighbourListRRelax_T  
  Real(kind=DoubleReal), Dimension(1:nlLenPerConfig*maxConfigsRelax,1:6) :: neighbourListCoordsRelax_T    
  
!------------------------------------------------------------------------------ 
! Precalc
!   
!


!------------------------------------------------------------------------------ 
! Evaluate EAM
!   
!
  Type(rssConfig), Dimension(1:maxConfigs) :: rssConfigsArr
  Type(rssConfig) :: rssConfigsArrTotal


!------------------------------------------------------------------------------ 
! CalcEAM
!   
!    
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configCalcEnergies  
  Real(kind=DoubleReal), Dimension(1:atomsPerConfig*maxConfigs,1:3) :: configCalcForces  
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:9) :: configCalcStresses
  
  
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:20) :: configCalc  
    
  
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configCalcEV
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configCalcEE
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configCalcEL
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configCalcBM
  
  Real(kind=DoubleReal) :: maxDensity


!------------------------------------------------------------------------------ 
! Evaluate BP
!   
!
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP) :: undistortedCellEnergies = 0.0D0
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:aLatSamples) :: alatEnergies = 0.0D0 
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:aLatSamples) :: configVolBP = 0.0D0
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:ecSamples) :: ecEnergiesBP = 0.0D0 
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:ecSamples) :: ecStrainBP = 0.0D0 
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:2*ecSamples+1) :: ecEnergiesBP_T = 0.0D0 
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP,1:2*ecSamples+1) :: ecStrainBP_T = 0.0D0 
  Type(rssBP), Dimension(1:maxConfigsBP) :: rssBPArr
  Type(rssBP) :: rssBPArrTotal
    

!------------------------------------------------------------------------------ 
! BpCalcEAM
!   
!  
  Real(kind=DoubleReal), Dimension(1:maxConfigsBP) :: configCalcEnergiesBP  
  
      

!------------------------------------------------------------------------------ 
! relaxCalcEAM
!   
!  
  Real(kind=DoubleReal), Dimension(1:maxConfigsRelax) :: configCalcEnergiesRelax  
  Real(kind=DoubleReal), Dimension(1:atomsPerConfig*maxConfigsRelax,1:3) :: configCalcForcesRelax 
  Real(kind=DoubleReal), Dimension(1:atomsPerConfig*maxConfigsRelax,1:2) :: configAtomEnergyRelax
  
  
  
!------------------------------------------------------------------------------ 
! Optimise
!   
!   
  Type(saConfig), Dimension(1:20) :: saConfigIn
  Integer(kind=StandardInteger) :: optLogCounter
  Logical :: optForceZBL
! LMA calc and ref values
  Real(kind=DoubleReal), Dimension(1:100000,1:2) :: calcRef
  Integer(kind=StandardInteger) :: countCalcRef
  Integer(kind=StandardInteger) :: crCount
  Integer(kind=StandardInteger) :: optEmbeddingFit
  Integer(kind=StandardInteger) :: optDensityFit

! -----------------------
! Default variables
  
! Declare variables - debug options

! Declare variables - run options



! ----------------------------------------------
! Read Input
! ----------------------------------------------

! ----------------------------------------------
! Read EAM
! ----------------------------------------------

! ----------------------------------------------
! Read Config
! ----------------------------------------------


! -----------------------
! Read EAM File + Read Configuration File    < 20MB

! Read Configuration File
  Integer(kind=StandardInteger) :: configCountT, configCountRI                      ! 4bit
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:20) :: configurationsI         ! 41KB        1 xcopy, 2 ycopy, 3 zcopy, 4 forces,
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:30) :: configurationsR                 ! 164KB       1 lp, 2-10 xx-zz, 11 rc, 12 vc   21-29 configUnitVector


! Configuration Reference/Calculated Values
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:20) :: configRef                       !             1 Energy PA, 2 EqVol
                    !             1 Energy, 2 EqVol    (maybe)
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configRefForces                !       1 fx, 2 fy, 3 fz

  Real(kind=DoubleReal), Dimension(1:100000,1:2) :: configAtomEnergy
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:9) :: configRefStresses
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configRefEnergies

  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configRefEV                          ! Equilibrium volume


  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configRefBM

  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:10) :: configRSS                       ! 1 energy, 2 forces, 3 stresses
  Real(kind=DoubleReal), Dimension(1:20) :: testConfigRSS                          ! 1FccALat,2FccEMin,3FccBM,4FccEoS,5FccC11,6FccC12,7FccC44,8BccALat,9BccEMin,10BccBM,11BccEoS,12BccC11,13BccC12,14BccC44
  Real(kind=DoubleReal) :: totalRSS
  Real(kind=DoubleReal) :: optimumRSS, startRSS, configTotalRSS, bestRSS
! Optimisation
  Real(kind=DoubleReal) :: nodeVariationAmount
  Real(kind=DoubleReal) :: saTemp, saMaxVariation
  Integer(kind=StandardInteger) :: saTempLoops, saVarLoops
  Integer(kind=StandardInteger) :: varyFixedNodes
  Integer(kind=StandardInteger) :: jumbleNodesOpt
  Integer(kind=StandardInteger) :: forceEmbeFitOpt
! DFT Config
  Character(len=8), Dimension(1:10,1:2) :: dftReplaceLabel

! ----------------------------------------------
! Neighbour List
! ----------------------------------------------


! ----------------------------------------------
! Calculations
! ----------------------------------------------
  Real(kind=DoubleReal) :: forceStressSwitch
  Real(kind=DoubleReal), Dimension(1:50000) :: calculationDensity    ! Electron density, D-band
  Real(kind=DoubleReal), Dimension(1:50000) :: calculationDensityS   ! S-Band
  Real(kind=DoubleReal), Dimension(1:100000) :: pairForce

! ----------------------------------------------
! Results
! ----------------------------------------------
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: calcConfigEnergies
  Real(kind=DoubleReal), Dimension(1:maxConfigs*atomsPerConfig,1:3) :: calcConfigForces
  
  
! ----------------------------------------------
! MD
! ----------------------------------------------
  Real(kind=DoubleReal), Dimension(1:maxConfigs*atomsPerConfig,1:3) :: mdAcceleration
  Real(kind=DoubleReal), Dimension(1:maxConfigs*atomsPerConfig,1:3) :: mdVelocity
  Real(kind=DoubleReal), Dimension(1:maxConfigs*atomsPerConfig,1:3) :: mdCoordChange
  

! ----------------------------------------------
! EAM Testing
! ----------------------------------------------

  Real(kind=DoubleReal), Dimension(1:12) :: fccCalcValues, bccCalcValues
  Real(kind=DoubleReal) :: fccALat, fccEMin, fccVolMin, fccBM, fccBMP
  Real(kind=DoubleReal) :: bccALat, bccEMin, bccVolMin, bccBM, bccBMP
  Real(kind=DoubleReal), Dimension(1:10) :: fccEC, bccEC
  Real(kind=DoubleReal) :: fccALatMurn, fccEMinMurn, fccVolMinMurn, fccBMMurn, fccBMPMurn
  Real(kind=DoubleReal) :: bccALatMurn, bccEMinMurn, bccVolMinMurn, bccBMMurn, bccBMPMurn
  Real(kind=DoubleReal) :: fccALatBirchMurn, fccEMinBirchMurn, fccVolMinBirchMurn, fccBMBirchMurn, fccBMPBirchMurn
  Real(kind=DoubleReal) :: bccALatBirchMurn, bccEMinBirchMurn, bccVolMinBirchMurn, bccBMBirchMurn, bccBMPBirchMurn
  Real(kind=DoubleReal), Dimension(1:10) :: fccECMurn, bccECMurn
  Integer(kind=StandardInteger) :: printTestingData, outputTestingData
  Real(kind=DoubleReal) :: testingRSS, eosFitRSS
  Integer(kind=StandardInteger) :: testingFitChoice  ! 1 Murn, 2 BirchMurn
  Integer(kind=StandardInteger) :: eosFitRSSOption  ! Include RSS of EAM model to ref EoS
  Real(kind=DoubleReal) :: testingALatRSS, testingEMinRSS
  Real(kind=DoubleReal) :: testingBMRSS, testingECRSS

! ----------------------------------------------
! Input Config Neighbour List
! ----------------------------------------------
  Integer(kind=StandardInteger) :: neighbourListCountInput
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:3) :: neighbourListKeyInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:1) :: neighbourListKeyRInput
  Integer(kind=StandardInteger), Dimension(1:800000,1:6) :: neighbourListIInput
  Real(kind=DoubleReal), Dimension(1:800000) :: neighbourListRInput
  Real(kind=DoubleReal), Dimension(1:800000,12) :: neighbourListCoordsInput
  
  Integer(kind=StandardInteger) :: configCountInput
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:20) :: configurationsIInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:30) :: configurationsRInput
  Integer(kind=StandardInteger) :: coordCountInput
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:3) :: configurationCoordsKeyInput
  Integer(kind=StandardInteger), Dimension(1:50000,1:1) :: configurationCoordsIInput
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationCoordsRInput
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationForcesRInput
  Integer(kind=StandardInteger) :: coordCountGInput
  Integer(kind=StandardInteger), Dimension(1:maxConfigs,1:3) :: configurationCoordsKeyGInput
  Integer(kind=StandardInteger), Dimension(1:100000,1:1) :: configurationCoordsIGInput
  Real(kind=DoubleReal), Dimension(1:100000,1:6) :: configurationCoordsRGInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configVolumeInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:20) :: configRefInput
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configRefForcesInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs,1:9) :: configRefStressesInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configRefEnergiesInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configRefEVInput
  Real(kind=DoubleReal), Dimension(1:maxConfigs) :: configRefBMInput

! Times
  Real(kind=DoubleReal) :: timeStart, timeEnd, timeDuration
  Real(kind=DoubleReal) :: globInitTime, eamLoadTime, nlTime, nlTimeBP, configLoadTime
  Real(kind=DoubleReal) :: efsCalcTime, efsCalcTimeBP, evalTimeBP

  Private

!------------------------------------------------------------------------------ 
! Global Parameters
!

! Subroutine
  Public :: initGlobals
  Public :: storeTime
! Initialise subroutine variables
  Public :: compileLine
  Public :: programStartTime, programEndTime
  Public :: cpuTime, cpuTimeLabels
  Public :: outputFile
  Public :: outputFileEnergies
  Public :: outputFileForces
  Public :: currentWorkingDirectory
  Public :: fileCleanupList
  Public :: outputDirectory
  Public :: tempDirectory
! MPI Variables
  Public :: mpiProcessCount
  Public :: mpiProcessID
  Public :: quietOverride
  
! System Variables
  Public :: largeArraySize
  Public :: processMap, processMapBP, processMapE
  Public :: maxConfigs, maxConfigsBP, maxConfigsRelax
  Public :: bpCopies, fccAlatBP, bccAlatBP, bpCutoffNL, bpCutoff
  Public :: aLatSamples, ecSamples
! Default variables
  Public :: eamFunctionTypes
! Set defaults - debug options
  

!------------------------------------------------------------------------------ 
! Program Configuration/Settings
!
  Public :: eampaSettings  

  
  Public :: eamForceSpline
  Public :: eamForceZBL
! Config Details - User Input
  Public :: globalConfigUnitVector
  Public :: configFilePath, configFilePathT
  Public :: saveExpConfigFile
! BP Config Details - User Input
  Public :: bpConfigFilePath, bpPrintData
! DFT Settings
  Public :: dftElement
  Public :: dftOptEnergy
  Public :: dftCohEnergy
! Neighbour List Settings
  Public :: nlCutoff
  Public :: nlTestCutoff
! Calculation details
  Public :: calcEqVol
  Public :: refineEqVol
  Public :: saveForcesToFile, saveNLToFile
  Public :: eosChart
! Optimise options
  Public :: anOptRandLoops
  Public :: varyNodeOptions
  Public :: optRunType
  Public :: optFrom
!  
  Public :: optLoops
  Public :: reduceNodes
  Public :: embeRescale
  Public :: forceEmbeFitOpt
! RSS calculation options
  Public :: rssWeighting
  Public :: configWeighting
! PW Batch Files - User Input
  Public :: pwbRunType
  Public :: pwbConfigFilePath, pwbConfigFilePathT
  Public :: pwbBatchDir
  Public :: pwbVarianceSwitch
  Public :: pwbVarianceType
  Public :: pwbVarianceMax
  Public :: pwbVarianceSigma
  Public :: pwbInterstitialAtom
  Public :: pwbInterstitialDetails



! ----------------------------------------------
! Read Config
! ----------------------------------------------
  Public :: configInputData
  Public :: configInputDataDFT
  Public :: configInputDataTemp
  Public :: configInputDataDFTTemp
  Public :: configLabelReplace
  Public :: crystalUnitCell
  Public :: configsAtomTotal

! -----------------------
! Read EAM File + Read Configuration File
  Public :: elements
  Public :: elementsCount
  Public :: elementsCharge
! Analytic potential functions  
  Public :: apfData, apfDataIn, apfDataOpt
! Read EAM File
  Public :: eamFunctionCount, eamPairCount, eamDensCount, eamEmbeCount
  Public :: eamDdenCount, eamSdenCount, eamDembCount, eamSembCount
  Public :: eamType
  Public :: eamKey
  Public :: eamData
  Public :: splineNodesData, splineNodesDataTemp
  Public :: splineNodesKey, splineNodesKeyTemp
  Public :: splineNodesResponse
  Public :: eamKeyInput
  Public :: eamDataInput
  Public :: eamKeyOpt
  Public :: eamDataOpt
  Public :: splineNodesKeyOpt, splineNodesKeyBest
  Public :: splineNodesDataOpt, splineNodesDataBest
! Read Configuration File
  Public :: configCount, configCountT, configCountRI
  Public :: configurationsI, configurationsR
  Public :: coordCount
  Public :: configurationCoordsKey, configurationCoordsI
  Public :: configurationCoordsR, configurationForcesR
  Public :: coordCountG
  Public :: configurationCoordsKeyG, configurationCoordsIG
  Public :: configurationCoordsRG
  Public :: configVolume
  Public :: configVolumeOpt


!------------------------------------------------------------------------------ 
! Read Input
!   
! user file array
  Public :: userInputData
! Verbose
  Public :: printToTerminal, eampaRunType
! Input File - User Input  
  Public :: inputFilePath
! EAM Details - User Input
  Public :: eamFilePath, eamNodesFilePath, eamSaveFile, potentialType
  Public :: eamInterpPoints, zblHardCore, splineNodeCount
  Public :: splineTotalNodes, eamMakeAlloy, eamFileType
  Public :: makeEAMCharts
  Public :: refBulkProperties

  
!------------------------------------------------------------------------------ 
! Read EAM File
!   
! EAM file array  
  Public :: eamInputData
    
!------------------------------------------------------------------------------ 
! Read BP Config File
!  
  Public :: bpInArr, bpConfigInputData   
  
!------------------------------------------------------------------------------ 
! Bulk Properties Config
!     
  Public :: configCountBP, configurationCoordsKeyBP, configurationCoordsIBP
  Public :: configurationCoordsRBP, configsAtomTotalBP, aLatBP, bpUnitCellCount
  Public :: calcBulkProperties
  
!------------------------------------------------------------------------------ 
! Relax Config
!     
  Public :: configCountRelax, configurationCoordsKeyRelax
  Public :: configurationCoordsIRelax, configurationCoordsRRelax
  Public :: configsAtomTotalRelax, aLatRelax
  Public :: relaxUnitCellCount
  
!------------------------------------------------------------------------------ 
! Neighbour List
!   
  Public :: nlUniqueKeys, neighbourListCount
  Public :: neighbourListKey, neighbourListKeyR
  Public :: neighbourListI, neighbourListR
  Public :: neighbourListCoords, atomSeparationSpread
    
!------------------------------------------------------------------------------ 
! Neighbour List for Bulk Property Configs
! 
  Public :: nlUniqueKeysBP, neighbourListCountBP
  Public :: neighbourListKeyBP, neighbourListKeyRBP
  Public :: neighbourListIBP, neighbourListRBP
  Public :: neighbourListCoordsBP, atomSeparationSpreadBP  
! Temporary NL arrays
  Public :: neighbourListRBP_T, neighbourListCoordsBP_T    
      
!------------------------------------------------------------------------------ 
! Neighbour List for Relax Configs
! 
  Public :: nlUniqueKeysRelax, neighbourListCountRelax
  Public :: neighbourListKeyRelax, neighbourListKeyRRelax
  Public :: neighbourListIRelax, neighbourListRRelax
  Public :: neighbourListCoordsRelax, atomSeparationSpreadRelax  
  Public :: neighbourListRRelax_T, neighbourListCoordsRelax_T
  
!------------------------------------------------------------------------------ 
! Precalc
!   
!



!------------------------------------------------------------------------------ 
! Evaluate EAM
!   
!
  Public :: rssConfigsArr
  Public :: rssConfigsArrTotal


!------------------------------------------------------------------------------ 
! CalcEAM
!   
!      
 Public ::  maxDensity
  
  
!------------------------------------------------------------------------------ 
! Evaluate BP
!   
!
  Public :: undistortedCellEnergies
  Public :: alatEnergies, configVolBP
  Public :: ecEnergiesBP, ecStrainBP
  Public :: ecEnergiesBP_T, ecStrainBP_T
  Public :: rssBPArr 
  Public :: rssBPArrTotal 

!------------------------------------------------------------------------------ 
! BpCalcEAM
!   
!    
  Public :: configCalcEnergiesBP  

!------------------------------------------------------------------------------ 
! relaxCalcEAM
!   
!    
  Public :: configCalcEnergiesRelax  
  Public :: configCalcForcesRelax  
  Public :: configAtomEnergyRelax  
  
!------------------------------------------------------------------------------ 
! Optimise
!   
!   
  Public :: saConfigIn 
  Public :: optLogCounter
  Public :: optForceZBL
! LMA calc and ref values
  Public :: calcRef
  Public :: countCalcRef
  Public :: crCount
  Public :: optEmbeddingFit, optDensityFit
  
! Configuration Reference/Calculated Values
  Public :: configRef
  Public :: configCalc
  Public :: configRefForces
  Public :: configCalcForces
  Public :: configAtomEnergy
  Public :: configRefStresses
  Public :: configCalcStresses
  Public :: configRefEnergies
  Public :: configCalcEnergies
  Public :: configRefEV
  Public :: configCalcEV
  Public :: configCalcEE
  Public :: configCalcEL
  Public :: configRefBM
  Public :: configCalcBM
  Public :: configRSS, configTotalRSS, testConfigRSS
  Public :: totalRSS, optimumRSS, startRSS, bestRSS
  
  
  
  
! Optimisation
  Public :: nodeVariationAmount
  Public :: saTemp
  Public :: saTempLoops, saVarLoops, saMaxVariation
  Public :: varyFixedNodes
  Public :: jumbleNodesOpt
! Neighbour List

! Temporary NL arrays
  !Public :: neighbourListIT
  !Public :: neighbourListRT
  !Public :: neighbourListCoordsT
! Calculations
  Public :: forceStressSwitch
  Public :: calculationDensity
  Public :: calculationDensityS
  Public :: pairForce
  
! ----------------------------------------------
! Results
! ----------------------------------------------

  Public :: calcConfigEnergies, calcConfigForces
  
    
! ----------------------------------------------
! MD
! ----------------------------------------------
  
  Public :: mdAcceleration, mdVelocity, mdCoordChange
  
  
! EAM Testing
  Public :: fccCalcValues, bccCalcValues
  Public :: fccALat, fccEMin, fccVolMin, fccBM, fccBMP
  Public :: bccALat, bccEMin, bccVolMin, bccBM, bccBMP
  Public :: fccEC, bccEC
  Public :: fccALatMurn, fccEMinMurn, fccVolMinMurn, fccBMMurn, fccBMPMurn
  Public :: bccALatMurn, bccEMinMurn, bccVolMinMurn, bccBMMurn, bccBMPMurn
  Public :: fccALatBirchMurn, fccEMinBirchMurn, fccVolMinBirchMurn, fccBMBirchMurn, fccBMPBirchMurn
  Public :: bccALatBirchMurn, bccEMinBirchMurn, bccVolMinBirchMurn, bccBMBirchMurn, bccBMPBirchMurn
  Public :: fccECMurn, bccECMurn
  Public :: printTestingData, outputTestingData
  Public :: testingRSS
  Public :: testingFitChoice
  Public :: eosFitRSS
  Public :: eosFitRSSOption
  Public :: testingALatRSS, testingEMinRSS
  Public :: testingBMRSS, testingECRSS
! DFT Config
  Public :: dftReplaceLabel

! ----------------------------------------------

! Input Config Neighbour List
  Public :: neighbourListCountInput
  Public :: neighbourListKeyInput
  Public :: neighbourListKeyRInput
  Public :: neighbourListIInput
  Public :: neighbourListRInput
  Public :: neighbourListCoordsInput
  Public :: configCountInput
  Public :: configurationsIInput
  Public :: configurationsRInput
  Public :: coordCountInput
  Public :: configurationCoordsKeyInput
  Public :: configurationCoordsIInput
  Public :: configurationCoordsRInput
  Public :: configurationForcesRInput
  Public :: coordCountGInput
  Public :: configurationCoordsKeyGInput
  Public :: configurationCoordsIGInput
  Public :: configurationCoordsRGInput
  Public :: configVolumeInput
  Public :: configRefInput
  Public :: configRefForcesInput
  Public :: configRefStressesInput
  Public :: configRefEnergiesInput
  Public :: configRefEVInput
  Public :: configRefBMInput

! Times
  Public :: timeStart, timeEnd, timeDuration
  Public :: globInitTime, eamLoadTime, nlTime, nlTimeBP, configLoadTime
  Public :: efsCalcTime, efsCalcTimeBP, evalTimeBP

  Contains

! Init global variables
  Subroutine initGlobals()
    Implicit None
    Real(kind=DoubleReal) :: globalsTimeStart, globalsTimeEnd
! Global Init Start time
    Call cpu_time(globalsTimeStart)
! Initialise Subroutine Variable
    compileLine = "01:00:08  22/12/2015"
    PROGRAMEndTime = 0.0D0
      quietOverride = .false.
      timeStart = 0.0D0
      timeEnd = 0.0D0
      timeDuration = 0.0D0
      nlTime = 0.0D0
      nlTimeBP = 0.0D0
      configLoadTime = 0.0D0
      efsCalcTime = 0.0D0
      efsCalcTimeBP = 0.0D0
      evalTimeBP = 0.0D0
      cpuTime = 0.0D0
      cpuTimeLabels(1) = "Globals Init"
      cpuTimeLabels(2) = "Evaluation Calculations"
      cpuTimeLabels(3) = "Equilibrium Volume/Energy"
      cpuTimeLabels(4) = "Read EAM"
      cpuTimeLabels(5) = "Optimise Potential Functions"
      cpuTimeLabels(6) = "E-F-S Calculations"
      cpuTimeLabels(7) = "Read Configs"
      cpuTimeLabels(8) = "Read User Input"
      cpuTimeLabels(9) = "Make Configs"
      cpuTimeLabels(10) = "Make Neighbour List"
      cpuTimeLabels(100) = "Program Time"
      currentWorkingDirectory = BlankString(currentWorkingDirectory)
      outputFile = BlankString(outputFile)
      outputFileEnergies = BlankString(outputFileEnergies)
      outputFileForces = BlankString(outputFileForces)
      fileCleanupList = BlankStringArray(fileCleanupList)
      outputDirectory = BlankString(outputDirectory)
      tempDirectory = BlankString(tempDirectory)
! MPI Global Variables
      mpiProcessCount = 0
      mpiProcessID = 0
! System Variables
      largeArraySize = 0.0D0
      processMap = -1
      processMapBP = -1
      processMapE = -1
! Debug options
      printToTerminal = .false.
! Run options
      eampaRunType = BlankString(eampaRunType)
      
! Input File - User Input
      inputFilePath = BlankString(inputFilePath)
! ----------------------------------------------
! Read Input
! ----------------------------------------------
      userInputData = BlankStringArray(userInputData)
      makeEAMCharts = .true.
      eosChart = .false.
! ----------------------------------------------
! EAM Input
! ----------------------------------------------
      eamInputData = BlankStringArray(eamInputData)
      elements = "ZZ"
      eamType = 0
! EAM Details - User Input
      eamFilePath = BlankString(eamFilePath)
      eamSaveFile = BlankString(eamSaveFile)
      eamNodesFilePath = BlankString(eamNodesFilePath)
      eamInterpPoints = 4
      zblHardCore = 0.0D0
      splineNodeCount = 0
      splineTotalNodes = 0
      eamForceSpline = .false.
      eamForceZBL = .false.
      eamMakeAlloy = BlankStringArray(eamMakeAlloy)
      eamFileType = 1
! ----------------------------------------------
! Read Config
! ----------------------------------------------
      configInputData = BlankStringArray(configInputData)
      configInputDataDFT = BlankStringArray(configInputDataDFT)
      configInputDataTemp = BlankStringArray(configInputDataTemp)
      configInputDataDFTTemp = BlankStringArray(configInputDataDFTTemp)
      configLabelReplace = BlankStringArray(configLabelReplace)
      crystalUnitCell = 0.0D0
! Config Details - User Input
      globalConfigUnitVector = 0.0D0
      configFilePath = BlankString(configFilePath)
      configFilePathT = BlankString(configFilePathT)
      saveExpConfigFile = BlankString(saveExpConfigFile)
! DFT Settings
      dftElement = BlankStringArray(dftElement)
      dftOptEnergy = 0.0D0
      dftCohEnergy = 0.0D0
! Neighbour List Settings
      nlCutoff = -1.0D0
! Calculation details
      calcEqVol = BlankString(calcEqVol)
      refineEqVol = "NO "
      saveNLToFile = .false.
      saveForcesToFile = .true.
    
!------------------------------------------------------------------------------ 
! Read BP Config File
!  
      bpConfigInputData = BlankStringArray(bpConfigInputData)         
      
!------------------------------------------------------------------------------ 
! Optimise
!   
!       
! Optimise options
      varyNodeOptions = 0.0D0
      optLoops = 1
      saTemp = 0.0D0
      saTempLoops = 0
      saVarLoops = 0
      reduceNodes = 0
! Calculations
      forceStressSwitch = 1        ! Calculate stress and force by default
      calculationDensity = 0.0D0
      pairForce = 0.0D0
      optLogCounter = 0
      optForceZBL = .true.
      optEmbeddingFit = 0
      optDensityFit = 0
! Global Init End Time
      Call cpu_time(globalsTimeEnd)
! Store time duration
      globInitTime = globalsTimeEnd-globalsTimeStart
    End Subroutine initGlobals

! Init global variables
    Subroutine storeTime(i,time)
      Implicit None
! Private variables
      Integer(kind=StandardInteger) :: i
      Real(kind=DoubleReal) :: time
      cpuTime(i) = cpuTime(i) + time
    End Subroutine storeTime

  End Module globals
