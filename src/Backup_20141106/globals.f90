Module globals

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Declare all global variables 

!----------------------------------------
! Updated: 12th Aug 2014
!----------------------------------------

! Setup Modules
  Use kinds
  Use msubs
  Use general
! Force declaration of all variables
  Implicit None
! Initialise Subroutine Variable
  Real(kind=DoubleReal) :: programStartTime, programEndTime
  Real(kind=DoubleReal) :: timeStart, timeEnd, timeDuration
  Real(kind=DoubleReal) :: globalsTimeStart, globalsTimeEnd
  Real(kind=DoubleReal), Dimension(1:100) :: cpuTime
  Character(Len=64), Dimension(1:100) :: cpuTimeLabels  
  Character(len=255) :: currentWorkingDirectory
  Character(len=255) :: outputFile
  Character(len=255) :: outputFileEnergies
  Character(len=255) :: outputFileForces
  Character(len=512), Dimension(1:100) :: fileCleanupList
  Character(len=255) :: outputDirectory
  Character(len=255) :: tempDirectory
! MPI Global Variables
  Integer(kind=StandardInteger) :: mpiProcessCount, mpiProcessID    
! System Variables  
  Real(kind=DoubleReal) :: largeArraySize
  Integer(kind=StandardInteger), Dimension(1:1024,1:10) :: processMap
! Parameters
  Integer(kind=StandardInteger), Parameter :: maxConfigs = 1024
  
!-----------------------  
! Default variables
  Character(len=4), Dimension(1:10) :: eamFunctionTypes
! Declare variables - debug options
  Integer(kind=StandardInteger) :: printToTerminal
! Declare variables - run options
  Character(len=4) :: eampaRunType
  Integer(kind=StandardInteger) :: optionReadConf
  Integer(kind=StandardInteger) :: optionMakeConf
  Integer(kind=StandardInteger) :: optionReadEAM
  Integer(kind=StandardInteger) :: optionRunPrep
  Integer(kind=StandardInteger) :: optionNeighbourList
  Integer(kind=StandardInteger) :: optionCalcEnergies
  Integer(kind=StandardInteger) :: optionEval
  Integer(kind=StandardInteger) :: optionEvalFull
  Integer(kind=StandardInteger) :: optionOptimise
  Integer(kind=StandardInteger) :: optionTestEAM
  Integer(kind=StandardInteger) :: optionRunPWBatch
  Integer(kind=StandardInteger) :: optionOutput
! Input File - User Input
  Character(len=255) :: inputFilePath
  Character(len=255) :: inputFilePathT
! MPI Options
  Integer(kind=StandardInteger) :: mpiEnergy
! EAM Details - User Input
  Character(len=255) :: eamFilePath
  Character(len=255) :: eamFilePathT
  Character(len=255) :: eamNodesFilePath
  Character(len=64) :: eamSaveFile
  Integer(kind=StandardInteger) :: eamInterpPoints
  Real(kind=DoubleReal), Dimension(1:6) :: zblHardCore                             ! 1 Pair ZBL end, 2 Pair Spline End, 3 Dens Value, 4 De
  Integer(kind=StandardInteger), Dimension(1:50) :: splineNodeCount
  Integer(kind=StandardInteger) :: splineTotalNodes
  Integer(kind=StandardInteger) :: eamForceSpline
  Integer(kind=StandardInteger) :: eamForceZBL
  Character(len=2), Dimension(1:10) :: eamMakeAlloy
  Integer(kind=StandardInteger) :: eamFileType
! Config Details - User Input  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: globalConfigUnitVector
  Character(len=255) :: configFilePath
  Character(len=255) :: configFilePathT
  Character(len=255) :: saveConfigFile
  Character(len=255) :: saveExpConfigFile
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
  Integer(kind=StandardInteger) :: saveForcesToFile
  Integer(kind=StandardInteger) :: saveNLToFile
! Optimise options  
  Real(kind=DoubleReal), Dimension(1:10) :: varyNodeOptions
  Integer(kind=StandardInteger) :: optLoops
  Integer(kind=StandardInteger) :: reduceNodes  
  Integer(kind=StandardInteger) :: embeRescale
! RSS calculation options  
  Real(kind=DoubleReal), Dimension(1:20) :: rssWeighting    
  Real(kind=DoubleReal), Dimension(1:1024) :: configWeighting
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
  
  
!-----------------------
! Read EAM File + Read Configuration File    < 20MB
  Character(len=2), Dimension(1:300) :: elements                                   ! 0.6KB    
  Integer(kind=StandardInteger) :: elementsCount                                                     
  Integer(kind=StandardInteger), Dimension(1:300) :: elementsCharge                ! 1.2KB
! Read EAM File 
  Integer(kind=StandardInteger) :: eamFunctionCount, eamPairCount, eamDensCount, eamEmbeCount
  Integer(kind=StandardInteger) :: eamDdenCount, eamSdenCount, eamDembCount, eamSembCount
  Integer(kind=StandardInteger) :: eamType                                         ! 4bit        1=EAM, 2=2BMEAM
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: eamKey                     ! 1.2KB       1 atomA, 2 atomB, 3 function/al type, 4 func start, 5 func length, 6 func end
  Real(kind=DoubleReal), Dimension(1:100000,1:4) :: eamData                        ! 3.2MB       1 x, 2 y(x), 3 y'(x), 4 y''(x) 
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: splineNodesKey
  Real(kind=DoubleReal), Dimension(1:10000,1:6) :: splineNodesData
  Real(kind=DoubleReal), Dimension(1:10000,1:2) :: splineNodesResponse
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: eamKeyInput
  Real(kind=DoubleReal), Dimension(1:100000,1:4) :: eamDataInput
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: eamKeyOpt
  Real(kind=DoubleReal), Dimension(1:100000,1:4) :: eamDataOpt
  Integer(kind=StandardInteger), Dimension(1:50,1:6) :: splineNodesKeyOpt
  Real(kind=DoubleReal), Dimension(1:10000,1:6) :: splineNodesDataOpt
  !Integer(kind=StandardInteger), Dimension(1:50,1:6) :: splineNodesKeyIn
  !Real(kind=DoubleReal), Dimension(1:10000,1:4) :: splineNodesDataIn
! Read Configuration File  
  Integer(kind=StandardInteger) :: configCount, configCountT, configCountRI                      ! 4bit
  Integer(kind=StandardInteger), Dimension(1:1024,1:20) :: configurationsI         ! 41KB        1 xcopy, 2 ycopy, 3 zcopy, 4 forces,   
  Real(kind=DoubleReal), Dimension(1:1024,1:30) :: configurationsR                 ! 164KB       1 lp, 2-10 xx-zz, 11 rc, 12 BM   21-29 configUnitVector
! Input coords and forces
  Integer(kind=StandardInteger) :: coordCount                                      ! 4bit
  Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: configurationCoordsKey   ! 13KB        1 start, 2 length, 3 end
  Integer(kind=StandardInteger), Dimension(1:50000,1:1) :: configurationCoordsI   !             1 atomID
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationCoordsR           !             1 x, 2 y, 3 z
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationForcesR           !             1 fx, 2 fy, 3 fz 
! Generated/expanded coords and forces 
  Integer(kind=StandardInteger) :: coordCountG                                     ! 4bit
  Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: configurationCoordsKeyG  ! 13KB        1 start, 2 length, 3 end
  Integer(kind=StandardInteger), Dimension(1:100000,1:1) :: configurationCoordsIG  !          1 atomID
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configurationCoordsRG          !       1 x, 2 y, 3 z
  Real(kind=DoubleReal), Dimension(1:1024) :: configVolume
  Real(kind=DoubleReal), Dimension(1:1024) :: configVolumeOpt
! Configuration Reference/Calculated Values
  Real(kind=DoubleReal), Dimension(1:1024,1:20) :: configRef                       !             1 Energy PA, 2 EqVol 
  Real(kind=DoubleReal), Dimension(1:1024,1:20) :: configCalc                      !             1 Energy, 2 EqVol    (maybe)
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configRefForces                !       1 fx, 2 fy, 3 fz
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configCalcForces
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configRefStresses
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configCalcStresses
  Real(kind=DoubleReal), Dimension(1:1024) :: configRefEnergies
  Real(kind=DoubleReal), Dimension(1:1024) :: configCalcEnergies
  Real(kind=DoubleReal), Dimension(1:1024) :: configRefEV                          ! Equilibrium volume
  Real(kind=DoubleReal), Dimension(1:1024) :: configCalcEV
  Real(kind=DoubleReal), Dimension(1:1024) :: configCalcEE
  Real(kind=DoubleReal), Dimension(1:1024) :: configCalcEL
  Real(kind=DoubleReal), Dimension(1:1024) :: configRefBM
  Real(kind=DoubleReal), Dimension(1:1024) :: configCalcBM
  Real(kind=DoubleReal), Dimension(1:1024,1:10) :: configRSS                       ! 1 energy, 2 forces, 3 stresses
  Real(kind=DoubleReal) :: totalRSS, optimumRSS, startRSS, configTotalRSS
! Optimisation  
  Real(kind=DoubleReal) :: nodeVariationAmount
  Real(kind=DoubleReal) :: saTemp, saMaxVariation
  Integer(kind=StandardInteger) :: saTempLoops, saVarLoops
  Integer(kind=StandardInteger) :: varyFixedNodes
  Integer(kind=StandardInteger) :: jumbleNodesOpt
  Integer(kind=StandardInteger) :: forceEmbeFitOpt
! DFT Config
  Character(len=8), Dimension(1:10,1:2) :: dftReplaceLabel 
  

!----------------------------------------------
! Neighbour List  
!----------------------------------------------
  Integer(kind=StandardInteger), Dimension(1:800000) :: nlUniqueKeys               ! 2.0MB
  Integer(kind=StandardInteger) :: neighbourListCount
  Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: neighbourListKey
  Real(kind=DoubleReal), Dimension(1:1024,1:1) :: neighbourListKeyR
  Integer(kind=StandardInteger), Dimension(1:800000,1:6) :: neighbourListI         ! 12.0MB
  Real(kind=DoubleReal), Dimension(1:800000) :: neighbourListR                     ! 4.0MB
  Real(kind=DoubleReal), Dimension(1:800000,12) :: neighbourListCoords             ! 48.0MB
  Integer(kind=StandardInteger), Dimension(1:2000) :: atomSeparationSpread
! Temporary NL arrays    
  Integer(kind=StandardInteger), Dimension(1:800000,1:6) :: neighbourListIT         ! 16.0MB
  Real(kind=DoubleReal), Dimension(1:800000) :: neighbourListRT                     ! 5.0MB
  Real(kind=DoubleReal), Dimension(1:800000,12) :: neighbourListCoordsT             ! 60.0MB
  
  
!----------------------------------------------
! Calculations  
!----------------------------------------------
  Real(kind=DoubleReal), Dimension(1:50000) :: calculationDensity
  Real(kind=DoubleReal), Dimension(1:100000) :: pairForce
  
!----------------------------------------------
! Results  
!---------------------------------------------- 
  Real(kind=DoubleReal), Dimension(1:1024) :: calcConfigEnergies
  Real(kind=DoubleReal), Dimension(1:200000,1:3) :: calcConfigForces
  
  
  
!----------------------------------------------
! EAM Testing  
!----------------------------------------------   
  Real(kind=DoubleReal), Dimension(1:12) :: fccReferenceValues, bccReferenceValues
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
  
!----------------------------------------------
! PWscf Batch File Globals  
!----------------------------------------------
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: pwbUnitVector, pwbUnitVectorWorking      
  Real(kind=DoubleReal) :: pwbLatticeParameter
  Integer(kind=StandardInteger) :: pwbXCopy, pwbYCopy, pwbZCopy  
  Real(kind=DoubleReal) :: pwbXd, pwbYd, pwbZd  
  Character(len=8), Dimension(1:1024)   :: pwbAtomLabelsInput
  Real(kind=DoubleReal), Dimension(1:1024,1:3)   :: pwbAtomCoordsInput  
  Character(len=8), Dimension(1:4096)   :: pwbAtomLabels, pwbAtomLabelsWorking                  ! 2x33KB
  Real(kind=DoubleReal), Dimension(1:4096,1:3)   :: pwbAtomCoords, pwbAtomCoordsWorking        ! 2x99KB
  Character(len=8), Dimension(1:128)   :: pwbAtomicSpeciesL 
  Character(len=64), Dimension(1:128)   :: pwbAtomicSpeciesPP    
  Real(kind=DoubleReal), Dimension(1:128)   :: pwbAtomicSpeciesDP
! PWscf options
  Character(len=64) :: pwbRestartMode, pwbCalculation, pwbOutDir, pwbPseudoDir, pwbPrefix,&
                       pwbDiskIO, pwbOccupations, pwbSmearing, &
                       pwbDiagonalization, pwbMixingMode, pwbIonDynamics, pwbCellDynamics, &
                       pwbKpoints             
  Character(len=6) ::  pwbTprnfor, pwbTstress
  Real(kind=DoubleReal) :: pwbEtotConvThr, pwbForcConvThr, pwbDegauss, pwbMixingBeta, &
                           pwbConvThr, pwbPress, pwbCellFactor
  Integer(kind=StandardInteger) :: pwbNstep, pwbIbrav, pwbNat, pwbNtyp, pwbEcutwfc, pwbEcutrho
  Integer(kind=StandardInteger) :: pwbNbnd, pwbFixedAtoms
! PWscf working variables  
  Integer(kind=StandardInteger) :: pwbNatWorking, pwbNtypWorking
  
!----------------------------------------------
! Input Config Neighbour List  
!----------------------------------------------
  Integer(kind=StandardInteger) :: neighbourListCountInput
  Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: neighbourListKeyInput
  Real(kind=DoubleReal), Dimension(1:1024,1:1) :: neighbourListKeyRInput
  Integer(kind=StandardInteger), Dimension(1:800000,1:6) :: neighbourListIInput
  Real(kind=DoubleReal), Dimension(1:800000) :: neighbourListRInput
  Real(kind=DoubleReal), Dimension(1:800000,12) :: neighbourListCoordsInput
  Integer(kind=StandardInteger) :: configCountInput                              
  Integer(kind=StandardInteger), Dimension(1:1024,1:20) :: configurationsIInput   
  Real(kind=DoubleReal), Dimension(1:1024,1:30) :: configurationsRInput      
  Integer(kind=StandardInteger) :: coordCountInput                                      
  Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: configurationCoordsKeyInput   
  Integer(kind=StandardInteger), Dimension(1:50000,1:1) :: configurationCoordsIInput  
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationCoordsRInput         
  Real(kind=DoubleReal), Dimension(1:50000,1:3) :: configurationForcesRInput         
  Integer(kind=StandardInteger) :: coordCountGInput                                    
  Integer(kind=StandardInteger), Dimension(1:1024,1:3) :: configurationCoordsKeyGInput 
  Integer(kind=StandardInteger), Dimension(1:100000,1:1) :: configurationCoordsIGInput 
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configurationCoordsRGInput        
  Real(kind=DoubleReal), Dimension(1:1024) :: configVolumeInput
  Real(kind=DoubleReal), Dimension(1:1024,1:20) :: configRefInput          
  Real(kind=DoubleReal), Dimension(1:100000,1:3) :: configRefForcesInput  
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configRefStressesInput
  Real(kind=DoubleReal), Dimension(1:1024) :: configRefEnergiesInput
  Real(kind=DoubleReal), Dimension(1:1024) :: configRefEVInput      
  Real(kind=DoubleReal), Dimension(1:1024) :: configRefBMInput
  
  
  Private    
!-----------------------
! Subroutine  
  Public :: initGlobals
  Public :: storeTime
! Initialise subroutine variables
  Public :: programStartTime, programEndTime
  Public :: timeStart, timeEnd, timeDuration  
  Public :: cpuTime, cpuTimeLabels
  Public :: globalsTimeStart, globalsTimeEnd
  Public :: outputFile          
  Public :: outputFileEnergies     
  Public :: outputFileForces      
  Public :: currentWorkingDirectory  
  Public :: fileCleanupList
  Public :: outputDirectory
  Public :: tempDirectory
!MPI Variables  
  Public :: mpiProcessCount
  Public :: mpiProcessID  
! System Variables  
  Public :: largeArraySize
  Public :: processMap
  Public :: maxConfigs
! Default variables
  Public :: eamFunctionTypes 
! Set defaults - debug options  
  Public :: printToTerminal
! Public Variables - run options
  Public :: eampaRunType
  Public :: optionReadEAM
  Public :: optionReadConf
  Public :: optionMakeConf
  Public :: optionRunPrep
  Public :: optionNeighbourList
  Public :: optionCalcEnergies
  Public :: optionEval
  Public :: optionEvalFull
  Public :: optionTestEAM
  Public :: optionOptimise
  Public :: optionRunPWBatch
  Public :: optionOutput
! Input File - User Input
  Public :: inputFilePath, inputFilePathT
! MPI Options
  Public :: mpiEnergy
! EAM Details - User Input
  Public :: eamFilePath, eamFilePathT
  Public :: eamNodesFilePath
  Public :: eamSaveFile
  Public :: eamInterpPoints
  Public :: zblHardCore
  Public :: splineNodeCount
  Public :: splineTotalNodes
  Public :: eamForceSpline
  Public :: eamForceZBL
  Public :: eamMakeAlloy
  Public :: eamFileType
! Config Details - User Input   
  Public :: globalConfigUnitVector
  Public :: configFilePath, configFilePathT
  Public :: saveConfigFile, saveExpConfigFile
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
! Optimise options    
  Public :: varyNodeOptions
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
!-----------------------  
! Read EAM File + Read Configuration File
  Public :: elements
  Public :: elementsCount
  Public :: elementsCharge
! Read EAM File  
  Public :: eamFunctionCount, eamPairCount, eamDensCount, eamEmbeCount
  Public :: eamDdenCount, eamSdenCount, eamDembCount, eamSembCount
  Public :: eamType
  Public :: eamKey
  Public :: eamData
  Public :: splineNodesData
  Public :: splineNodesKey
  Public :: splineNodesResponse
  Public :: eamKeyInput
  Public :: eamDataInput
  Public :: eamKeyOpt
  Public :: eamDataOpt
  Public :: splineNodesKeyOpt
  Public :: splineNodesDataOpt
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
! Configuration Reference/Calculated Values
  Public :: configRef                     
  Public :: configCalc                      
  Public :: configRefForces               
  Public :: configCalcForces
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
  Public :: configRSS, configTotalRSS
  Public :: totalRSS, optimumRSS, startRSS
! Optimisation  
  Public :: nodeVariationAmount  
  Public :: saTemp
  Public :: saTempLoops, saVarLoops, saMaxVariation
  Public :: varyFixedNodes
  Public :: jumbleNodesOpt
! Neighbour List
  Public :: nlUniqueKeys  
  Public :: neighbourListCount  
  Public :: neighbourListKey
  Public :: neighbourListKeyR
  Public :: neighbourListI        
  Public :: neighbourListR       
  Public :: neighbourListCoords   
  Public :: atomSeparationSpread
! Temporary NL arrays    
  Public :: neighbourListIT        
  Public :: neighbourListRT       
  Public :: neighbourListCoordsT 
! Calculations  
  Public :: calculationDensity
  Public :: pairForce
! Results  
  Public :: calcConfigEnergies
  Public :: calcConfigForces
! EAM Testing  
  Public :: fccReferenceValues, bccReferenceValues
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
  
  
!----------------------------------------------
! PWscf Batch File Globals 
  Public :: pwbUnitVector, pwbUnitVectorWorking      
  Public :: pwbLatticeParameter
  Public :: pwbXCopy, pwbYCopy, pwbZCopy  
  Public :: pwbXd, pwbYd, pwbZd  
  Public :: pwbAtomLabelsInput
  Public :: pwbAtomCoordsInput  
  Public :: pwbAtomLabels, pwbAtomLabelsWorking
  Public :: pwbAtomCoords, pwbAtomCoordsWorking 
  Public :: pwbAtomicSpeciesL, pwbAtomicSpeciesPP    
  Public :: pwbAtomicSpeciesDP
!PWscf options
  Public :: pwbRestartMode, pwbCalculation, pwbOutDir, pwbPseudoDir, pwbPrefix,&
            pwbDiskIO, pwbOccupations, pwbSmearing, &
            pwbDiagonalization, pwbMixingMode, pwbIonDynamics, pwbCellDynamics, &
            pwbKpoints             
  Public :: pwbTprnfor, pwbTstress
  Public :: pwbEtotConvThr, pwbForcConvThr, pwbDegauss, pwbMixingBeta, &
            pwbConvThr, pwbPress, pwbCellFactor
  Public :: pwbNstep, pwbIbrav, pwbNat, pwbNtyp, pwbEcutwfc, pwbEcutrho
  Public :: pwbNbnd, pwbFixedAtoms
! PWscf working variables  
  Public :: pwbNatWorking, pwbNtypWorking  
  
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
  
  
Contains
  
! Init global variables
  Subroutine initGlobals() 
    Implicit None  
! Global Init Start time
    Call cpu_time(globalsTimeStart)
! Initialise Subroutine Variable
    programEndTime = 0.0D0
    timeStart = 0.0D0
    timeEnd = 0.0D0
    timeDuration = 0.0D0
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
! Default variables    
    eamFunctionTypes = BlankStringArray(eamFunctionTypes)
    eamFunctionTypes(1) = "PAIR"
    eamFunctionTypes(2) = "DENS"
    eamFunctionTypes(3) = "EMBE"
    eamFunctionTypes(4) = "DDEN"
    eamFunctionTypes(5) = "SDEN"
    eamFunctionTypes(6) = "DEMB"
    eamFunctionTypes(7) = "SEMB"
! Debug options  
    printToTerminal = 0
! Run options
    eampaRunType = BlankString(eampaRunType)
    optionReadConf = 0
    optionMakeConf = 0
    optionReadEAM = 0
    optionRunPrep = 0  
    optionNeighbourList = 0
    optionCalcEnergies = 0
    optionEval = 0
    optionEvalFull = 0
    optionTestEAM = 0
    optionOptimise = 0
    optionRunPWBatch = 0  
    optionOutput = 0    
! Input File - User Input  
    inputFilePath = BlankString(inputFilePath)
    inputFilePathT = BlankString(inputFilePathT)
! MPI Options
    mpiEnergy = 0
! EAM Details - User Input  
    eamFilePath = BlankString(eamFilePath)
    eamFilePathT = BlankString(eamFilePathT)
    eamSaveFile = BlankString(eamSaveFile)
    eamNodesFilePath = BlankString(eamNodesFilePath)
    eamInterpPoints = 4
    zblHardCore = 0.0D0
    splineNodeCount = 0
    splineTotalNodes = 0
    eamForceSpline = 0
    eamForceZBL = 0
    eamMakeAlloy = BlankStringArray(eamMakeAlloy)
    eamFileType = 1
! Config Details - User Input   
    globalConfigUnitVector = 0.0D0 
    configFilePath = BlankString(configFilePath)
    configFilePathT = BlankString(configFilePathT)
    saveConfigFile = BlankString(saveConfigFile)
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
    saveForcesToFile = 0
    saveNLToFile = 0    
! Optimise options    
    varyNodeOptions = 0.0D0 
    optLoops = 1
    saTemp = 0.0D0
    saTempLoops = 0
    saVarLoops = 0
    reduceNodes = 0
! RSS calculation options  
    rssWeighting = 0.0D0
    configWeighting = 1.0D0
! PW Batch Files - User Input   
    pwbRunType = BlankString(pwbRunType)
    pwbConfigFilePath = BlankString(pwbConfigFilePath)
    pwbConfigFilePathT = BlankString(pwbConfigFilePathT)
    pwbBatchDir = BlankString(pwbBatchDir)
    pwbVarianceSwitch = 0
    pwbVarianceType = BlankString(pwbVarianceType)
    pwbVarianceMax = 0.0D0
    pwbVarianceSigma = 0.0D0
    pwbInterstitialAtom = 0
    pwbInterstitialDetails = BlankStringArray(pwbInterstitialDetails)
! Read EAM + Config    
    elements = "ZZ"
    elementsCount = 0
    elementsCharge = -1
! Read EAM File
    eamFunctionCount = 0
    eamPairCount = 0
    eamDensCount = 0
    eamEmbeCount = 0
    eamDdenCount = 0
    eamSdenCount = 0
    eamDembCount = 0
    eamSembCount = 0
    eamType = 1
    eamKey = -1
    eamData = 0.0D0    
    splineNodesKey = -1
    splineNodesData = 0.0D0 
    splineNodesResponse = 0.0D0 
    eamKeyInput = -1
    eamDataInput = 0.0D0
    eamKeyOpt = -1
    eamDataOpt = 0.0D0
    splineNodesKeyOpt = -1
    splineNodesDataOpt = 0.0D0
! Read Configuration File 
    configCount = 0
    configCountT = 0
    configCountRI = 0
    configurationsI = 0
    configurationsR = 0.0D0
    coordCount = 0
    configurationCoordsKey = 0
    configurationCoordsI = 0
    configurationCoordsR = 0.0D0
    configurationForcesR = -2.1D20
    coordCountG = 0
    configurationCoordsKeyG = 0
    configurationCoordsIG = 0
    configurationCoordsRG = 0.0D0
    configVolume = 0.0D0
    configVolumeOpt = -2.1D0
! Configuration Reference/Calculated Values
    configRef = -2.1D20                    
    configCalc = -2.1D20                      
    configRefForces = -2.1D20               
    configCalcForces = -2.1D20
    configRefStresses = -2.1D20
    configCalcStresses = -2.1D20  
    configRefEnergies = -2.1D20  
    configCalcEnergies = -2.1D20       
    configRefEV = -2.1D20    
    configCalcEV = -2.1D20      
    configCalcEE = -2.1D20     
    configCalcEL = -2.1D20  
    configRefBM = -2.1D20
    configCalcBM = -2.1D20
    configRSS = 0.0D0
    configTotalRSS = 0.0D0
    totalRSS = 0.0D0
    optimumRSS = 0.0D0
    startRSS = 0.0D0
! Optimisation  
    nodeVariationAmount = 0.0D0    
    saTemp = 100.0D0
    saTempLoops = 10
    saVarLoops = 100
    saMaxVariation = 0.0D0
    varyFixedNodes = 0
    jumbleNodesOpt = 0
    embeRescale = 0
    forceEmbeFitOpt = 0
! Neighbour List
    nlUniqueKeys = 0
    neighbourListCount = 0
    neighbourListKey = 0
    neighbourListKeyR = 0.0D0
    neighbourListI = 0
    neighbourListR = 0.0D0
    neighbourListCoords = 0.0D0
    atomSeparationSpread = 0
! Temporary NL arrays    
    neighbourListIT = 0        
    neighbourListRT = 0.0D0       
    neighbourListCoordsT = 0.0D0 
! Calculation
    calculationDensity = 0.0D0
    pairForce = 0.0D0
! Results  
    calcConfigEnergies = 0.0D0
    calcConfigForces = 0.0D0
! EAM Testing  
    fccReferenceValues = -2.1D20
    bccReferenceValues = -2.1D20
    fccCalcValues = -2.1D20
    bccCalcValues = -2.1D20
    fccALat = 0.0D0
    fccEMin = 0.0D0
    fccVolMin = 0.0D0
    fccBM = 0.0D0
    fccBMP = 0.0D0
    bccALat = 0.0D0
    bccEMin = 0.0D0
    bccVolMin = 0.0D0
    bccBM = 0.0D0
    bccBMP = 0.0D0
    fccEC = -2.1D20
    bccEC = -2.1D20
    fccALatMurn = 0.0D0
    fccEMinMurn = 0.0D0
    fccVolMinMurn = 0.0D0
    fccBMMurn = 0.0D0
    fccBMPMurn = 0.0D0
    bccALatMurn = 0.0D0
    bccEMinMurn = 0.0D0
    bccVolMinMurn = 0.0D0
    bccBMMurn = 0.0D0
    bccBMPMurn = 0.0D0
    fccALatBirchMurn = 0.0D0
    fccEMinBirchMurn = 0.0D0
    fccVolMinBirchMurn = 0.0D0
    fccBMBirchMurn = 0.0D0    
    fccBMPBirchMurn = 0.0D0
    bccALatBirchMurn = 0.0D0
    bccEMinBirchMurn = 0.0D0
    bccVolMinBirchMurn = 0.0D0
    bccBMBirchMurn = 0.0D0
    bccBMPBirchMurn = 0.0D0
    fccECMurn = -2.1D20
    bccECMurn = -2.1D20
    
    printTestingData = 0
    outputTestingData = 0
    testingRSS = 0.0D0
    testingFitChoice = 1
    eosFitRSS = 0.0D0
    eosFitRSSOption = 0
    testingALatRSS = 0.0D0
    testingEMinRSS = 0.0D0
    testingBMRSS = 0.0D0
    testingECRSS = 0.0D0
! DFT Config
    dftReplaceLabel = BlankString2DArray(dftReplaceLabel)
  
    
!----------------------------------------------
! PWscf Batch File Globals - default values 
!---------------------------------------------- 
    pwbUnitVector = 0.0D0    
    pwbAtomLabelsInput = "#BLANK##"
    pwbAtomCoordsInput = -2.1D20
    pwbAtomicSpeciesL = "#BLANK##"
    pwbAtomicSpeciesPP = "#BLANK##"
    pwbAtomicSpeciesDP = -2.1D20
    pwbAtomLabelsWorking = BlankStringArray(pwbAtomLabelsWorking)
    pwbAtomCoordsWorking = -2.1D20
!Default pwscf file values, text
    pwbRestartMode = "from_scratch"
    pwbCalculation = "scf"
    pwbOutDir = "/gpfs/bb/bxp912/scratch"
    pwbPseudoDir = "/gpfs/bb/bxp912/pseudopotentials"
    pwbPrefix = "pwbatchfile"
    pwbDiskIO = "low"
    pwbOccupations = "smearing"
    pwbSmearing = "mv"
    pwbDiagonalization = "david"
    pwbMixingMode = "TF"
    pwbIonDynamics = "bfgs"
    pwbCellDynamics = "bfgs"
    pwbKpoints = "2 2 2 1 1 1"
!Default pwscf file values, boolean
    pwbTprnfor = ".true."
    pwbTstress = ".true."
!Default pwscf file values, double precision
    pwbEtotConvThr = 1.0D-4
    pwbForcConvThr = 1.0D-3
    pwbDegauss = 0.07
    pwbMixingBeta = 0.7
    pwbConvThr = 1.0D-6
    pwbPress = 0.0
    pwbCellFactor = 2.0
!Default pwscf file values, integer
    pwbNstep = 50
    pwbIbrav = 14
    pwbNat = 0
    pwbNtyp = 0
    pwbNbnd = 0
    pwbEcutwfc = 110
    pwbEcutrho = 320
!Fixed atoms in coords
    pwbFixedAtoms = 0    
! Working variables    
    pwbNatWorking = 0
    pwbNtypWorking  = 0 
    
    
    
! System Variables

! Estimate memory of large arrays
    largeArraySize = size(fileCleanupList) * 1.0D0 * 512.0D0 + &     ! Character
                     size(processMap) * 1.0D0 * 4.0D0 + &            ! Int
                     size(eamFunctionTypes) * 1.0D0 * 4.0D0 + &
                     size(zblHardCore) * 1.0D0 * 8.0D0 + &           ! DP
                     size(splineNodeCount) * 1.0D0 * 4.0D0 + &       ! Int
                     size(globalConfigUnitVector) * 8.0D0 + &
                     size(varyNodeOptions) * 1.0D0 * 8.0D0 + &     
                     size(rssWeighting) * 1.0D0 * 4.0D0 + &                            
                     size(elements) * 1.0D0 + &
                     size(elementsCharge) * 4.0D0 + &
                     size(eamKey) * 4.0D0 + &
                     size(eamData) * 8.0D0 + &
                     size(splineNodesKey) * 4.0D0 + &
                     size(splineNodesData) * 8.0D0 + &
                     size(eamKeyInput) * 4.0D0 + &
                     size(eamDataInput) * 8.0D0 + &
                     size(eamKeyOpt) * 4.0D0 + &
                     size(eamDataOpt) * 8.0D0 + &
                     size(splineNodesKeyOpt) * 4.0D0 + &
                     size(splineNodesDataOpt) * 8.0D0 + &
                     size(configurationsI) * 4.0D0 + &
                     size(configurationsR) * 8.0D0 + &
                     size(configurationCoordsKey) * 4.0D0 + &
                     size(configurationCoordsI) * 4.0D0 + &
                     size(configurationCoordsR) * 8.0D0 + &
                     size(configurationForcesR) * 8.0D0 + &
                     size(configurationCoordsKeyG) * 4.0D0 + &
                     size(configurationCoordsIG) * 4.0D0 + &
                     size(configurationCoordsRG) * 8.0D0 + &
                     size(nlUniqueKeys) * 4.0D0 + &
                     size(neighbourListKey) * 4.0D0 + &
                     size(neighbourListI) * 4.0D0 + &
                     size(neighbourListR) * 8.0D0 + &
                     size(neighbourListCoords) * 8.0D0 + &
                     size(neighbourListIT) * 4.0D0 + &
                     size(neighbourListRT) * 8.0D0 + &
                     size(neighbourListCoordsT) * 8.0D0 + &
                     size(pwbUnitVector) * 8.0D0 + &
                     size(pwbUnitVectorWorking) * 8.0D0 + &
                     size(pwbAtomLabelsInput) * 1.0D0 * 8.0D0 + &
                     size(pwbAtomCoordsInput) * 8.0D0 + &
                     size(pwbAtomLabels) * 1.0D0 * 8.0D0 + &
                     size(pwbAtomLabelsWorking) * 1.0D0 * 8.0D0 + &
                     size(pwbAtomCoords) *8.0D0 + &
                     size(pwbAtomCoordsWorking) *8.0D0 + &
                     size(pwbAtomicSpeciesL) * 1.0D0 * 8.0D0 + &
                     size(pwbAtomicSpeciesPP) * 1.0D0 * 64.0D0 + &
                     size(pwbAtomicSpeciesDP) * 8.0D0 + &
                     size(calculationDensity) * 8.0D0 + &
                     size(configRef) * 8.0D0 + &
                     size(configCalc) * 8.0D0 + &
                     size(configRefForces) * 8.0D0 + &
                     size(configCalcForces) * 8.0D0 + &
                     size(configRefStresses) * 8.0D0 + &
                     size(configCalcStresses) * 8.0D0
                     
                     
    largeArraySize = largeArraySize / 1.0D6                 
                     
    
! Global Init End Time
    Call cpu_time(globalsTimeEnd)    
! Store time duration
    Call storeTime(1,globalsTimeEnd-globalsTimeStart)   
  
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