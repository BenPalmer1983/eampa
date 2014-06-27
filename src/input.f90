Module input

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use strings		!string functions
  Use maths
  Use units
  Use initialise


!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'    
!declare global variables  
  Character(len=255)  :: inputFileName
  Character(len=255)  :: potentialFilePath
  Character(len=255)  :: potentialFilePathTemp
  Character(len=255)  :: configurationsFilePath
  Character(len=2), Dimension( : ), Allocatable :: elements
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: elementsCharge
  Integer(kind=StandardInteger) :: elementCount
  Integer(kind=StandardInteger) :: eamType
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKey
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamData
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: unitVector, globalUnitVector
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: configHeaderI
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configHeaderR
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configCoordsI
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configCoordsR
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configForcesR
  Integer(kind=StandardInteger), Parameter :: headerWidth = 14 
  Integer(kind=StandardInteger) :: numberPotentials
  Integer(kind=StandardInteger) :: pairCount, densCount, dendCount, embeCount, embdCount
!user input variables  
  Character(len=1) :: saveFileCoords,saveFileNeighbourList,saveFilePot
  Character(len=1) :: saveFileForces, saveFileBM
  Character(len=255) :: eamPreparedFile
  Character(len=255) :: dlpolyFileDir
  Character(len=255) :: pwbOptConfFile
  Character(len=255) :: pwbOutputDir
  Integer(kind=StandardInteger) :: potType
  Integer(kind=StandardInteger) :: calcRunType
  Integer(kind=StandardInteger) :: optionReadConf, optionReadEAM
  Integer(kind=StandardInteger) :: optionRunProcesses, optionRunPrep, optionRunDLPrep
  Integer(kind=StandardInteger) :: optionRunPrepEAM, optionRunPWBatch  
  Integer(kind=StandardInteger) :: configCount
  Integer(kind=StandardInteger) :: eamInterpType
  Integer(kind=StandardInteger) :: calcForcesOnOff
  Integer(kind=StandardInteger) :: calcStressOnOff
  Integer(kind=StandardInteger) :: calcBMCalcOnOff  !0 only with bm supplied, 1 on for all configurations
  Integer(kind=StandardInteger) :: calcECCalcOnOff	!0 only with ec supplied, 1 on for all configurations
  Integer(kind=StandardInteger) :: calcEVCalcOnOff  !0 only with eq vol supplied, 1 on for all configurations
  Real(kind=SingleReal) :: dftInRadiusCutoff
  Character(len=3), Dimension( : , : ), Allocatable :: dftInReplaceSym
  Real(kind=SingleReal) :: dftInOptEnergy
  Real(kind=SingleReal) :: dftInCohEnergy
  Character(len=128), Dimension(1:500) :: pwscfFilesList
  Integer(kind=StandardInteger) :: saCycles, saTStart, saTEnd, saTIncs
  Real(kind=DoubleReal) :: saTemp, saSpreadFactor
!Isotope/Element Arrays
  Character(len=2), Dimension( : ), Allocatable :: isotopesChar
  Integer, Dimension( : , : ), Allocatable :: isotopesInt
  Real, Dimension( : , : ), Allocatable :: isotopesReal
  Character(len=2), Dimension( : ), Allocatable :: elementSymbol
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: elementAtomicMass
!Optimisation eam reduced points size
  Integer(kind=StandardInteger) :: eamReducedPointsV, eamReducedPointsP, eamReducedPointsF
  Integer(kind=StandardInteger) :: eamTrialPointsV, eamTrialPointsP, eamTrialPointsF
!Optimization weights
  Integer(kind=StandardInteger) :: eamEnergyOptWeight, eamStressOptWeight, eamForceOptWeight
  Integer(kind=StandardInteger), Dimension(1:10) :: eamOptWeights
  Integer(kind=StandardInteger) :: eamOptWeightType
!zbl hard core
  Real(kind=DoubleReal) :: eamZBLPairLower, eamZBLPairUpper, eamZBLDensZero, eamZBLDensCutoff
  Real(kind=DoubleReal) :: eamZBLEmbeZero, eamZBLEmbeCutoff
!output
  Integer(kind=StandardInteger) :: printToTerminal

  
!Privacy of functions/subroutines/variables
  Private
!Variables - all run types       
  Public :: calcRunType 
  Public :: optionReadConf, optionReadEAM, optionRunProcesses 
  Public :: optionRunPrep, optionRunDLPrep, optionRunPrepEAM, optionRunPWBatch
  Public :: inputFileName
!Variables - Run calculation + Pot Prep
  Public :: potentialFilePath, potentialFilePathTemp
  Public :: elements 
  Public :: elementsCharge
  Public :: eamZBLPairLower, eamZBLPairUpper, eamZBLDensZero, eamZBLDensCutoff,&
            eamZBLEmbeZero, eamZBLEmbeCutoff
!Variables - Run calculation only
  Public :: configurationsFilePath	             	
  Public :: unitVector, globalUnitVector   
  Public :: configHeaderI, configHeaderR    
  Public :: configCoordsI, configCoordsR, configForcesR  
  
!prep eam potential only
  Public :: eamPreparedFile
     
!prep dlpoly potential only
  Public :: dlpolyFileDir
       	
!pwscf batch only		
  Public :: pwbOptConfFile
  Public :: pwbOutputDir
  
!save other data files
  Public :: saveFileCoords,saveFileNeighbourList,saveFilePot
  Public :: saveFileForces,saveFileBM
		
!Misc		
  Public :: headerWidth             
  Public :: numberPotentials        
  Public :: eamType           	     
  Public :: potType
  Public :: pairCount           	  
  Public :: densCount           	  
  Public :: dendCount           	   
  Public :: embeCount           	   
  Public :: embdCount           	  
  Public :: elementCount           	  
  Public :: eamKey               	   
  Public :: eamData            
  Public :: configCount     
  Public :: eamInterpType           
  Public :: calcForcesOnOff      
  Public :: calcStressOnOff
  Public :: calcBMCalcOnOff, calcECCalcOnOff, calcEVCalcOnOff
  Public :: dftInRadiusCutoff   
  Public :: dftInReplaceSym         
  Public :: dftInOptEnergy
  Public :: dftInCohEnergy
  Public :: elementSymbol
  Public :: isotopesChar		
  Public :: isotopesInt  		
  Public :: isotopesReal		
  Public :: eamReducedPointsV, eamReducedPointsP, eamReducedPointsF !Variable
  Public :: eamTrialPointsV, eamTrialPointsP, eamTrialPointsF
  Public :: eamEnergyOptWeight, eamStressOptWeight, eamForceOptWeight
  Public :: eamOptWeights
  Public :: eamOptWeightType
!User input variables
  Public :: saCycles, saTStart, saTEnd, saTIncs
  Public :: saTemp, saSpreadFactor
  Public :: readEamPot
  Public :: reReadEamPot
  Public :: printToTerminal
!Subroutines
  Public :: runInput	
  
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runInput()	
	!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
	Call readUserInput()             !Read input from command line i.e. input file name/location
	Call readIsotopeData()           !Read in isotope data from data/isotopes.dat
	Call readInputFile()		     !Read input file
	If(optionReadEAM.eq.1)Then
	  Call readEamPotPrep()            !Do for all run options
	End If  
	Call synchMpiProcesses()
	If(optionReadEAM.eq.1)Then
	  Call readEamPot()
	End If  
	If(optionReadConf.eq.1)Then
	  If(mpiProcessID.eq.0)Then
	    Call copyConfigFile()
	  End If  
	  Call synchMpiProcesses()
	  If(mpiProcessID.eq.0)Then
	    Call readDFTFiles()
	  End If
	  Call synchMpiProcesses()
	  Call readConfiguration()
	End If	
  End Subroutine runInput  
!================================================================================================================================================  
!================================================================================================================================================
!read in user input data, input from the command line
  Subroutine readUserInput()
  !force declaration of all variables
	Implicit None  
!Read in command line arguments
    call get_command_argument(1,inputFileName)
  End Subroutine readUserInput
!================================================================================================================================================  
!================================================================================================================================================
!Read in any data files (isotope data etc)  
!read in isotope data  
  Subroutine readIsotopeData()
!force declaration of all variables
	Implicit None	
!declare variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, isotopeCounter, fileRows
	Integer(kind=StandardInteger) :: maxZ
	Character(len=25) :: buffera, bufferb, bufferc, bufferd, buffere, bufferf,tempA
	Character(len=255) :: bufferLong
	Character(len=255) :: dataFilesDir,dataFilesFile,tempB, fileRow
!set default variables	
	dataFilesDir = "data"
!open output file
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
    !write(999,"(A36,F8.4)") "Read Isotope Data                   ",ProgramTime()	 
!read isotope file path from input file
    Open(UNIT=1,FILE=inputFileName) 
    Do i=1,maxFileRows 
	  !Read in line
	  Read(1,"(A255)",IOSTAT=ios) bufferLong
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if 
	  if(bufferLong(1:8).eq."#datadir")then
	    Read(1,*,IOSTAT=ios) bufferLong
	    dataFilesDir = TrimSpaces(bufferLong)
	  endif	  
    enddo
	CLOSE(1) !close file
!make file name/path	
    dataFilesFile = Trim(dataFilesDir)//"/"//"isotopes.dat"
!count isotopes
	isotopeCounter = 0
	fileRows = 0
	Open(UNIT=1,FILE=dataFilesFile) 
    do i=1,maxFileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera, bufferb
	  !Read(1,*,IOSTAT=ios) fileRow
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
	  !Read(fileRow,*) buffera
	  !Read(fileRow,*) bufferb
	  !Read(i,*) tempA
	  !tempB = "isotopeData["//trim(tempA)//"] = "//trim(fileRow)	  
	  !print *,tempB
	  fileRows = fileRows + 1
	  if(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Number")then
		isotopeCounter = isotopeCounter + 1
	  endif	  
    enddo
	CLOSE(1) !close file
!Allocate Arrays
    Allocate(isotopesChar(1:isotopeCounter))
    Allocate(isotopesInt(1:isotopeCounter,1:10))
    Allocate(isotopesReal(1:isotopeCounter,1:10))
!Read in data
	isotopeCounter = 0
	Open(UNIT=1,FILE=dataFilesFile) 
	maxZ = 0
    do i=1,fileRows 
	  !Read in line
	  Read(1,*,IOSTAT=ios) buffera, bufferb
	  !break out
	  if (ios /= 0) then
	    EXIT 
	  end if
	  if(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Number")then
		isotopeCounter = isotopeCounter + 1
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		read(bufferd,*) isotopesInt(isotopeCounter,1) 		
		if(isotopesInt(isotopeCounter,1).gt.maxZ)then
		  maxZ = isotopesInt(isotopeCounter,1)
		endif
	  endif	  
	  if(buffera(1:6).eq."Atomic".and.bufferb(1:6).eq."Symbol")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		isotopesChar(isotopeCounter) = StrToUpper(bufferd(1:2))		
	  endif	 
	  if(buffera(1:4).eq."Mass".and.bufferb(1:6).eq."Number")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		read(bufferd,*) isotopesInt(isotopeCounter,2) 		
	  endif	  
	  if(buffera(1:8).eq."Relative".and.bufferb(1:6).eq."Atomic")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere
		buffere = NumericOnly(buffere)
		if(buffere(1:5).eq."     ")then
		  isotopesReal(isotopeCounter,1) = 1.0 * isotopesInt(isotopeCounter,2)
		else
		  read(buffere,*) isotopesReal(isotopeCounter,1) 		
		endif
	  endif	  
	  if(buffera(1:8).eq."Isotopic".and.bufferb(1:11).eq."Composition")then
		!Re-read file line
		BACKSPACE(1)
		Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		bufferd = NumericOnly(bufferd)
		if(bufferd(1:5).eq."     ")then
		  isotopesReal(isotopeCounter,2) = 0.0
		else
		  read(bufferd,*) isotopesReal(isotopeCounter,2) 		
		endif
	  endif	
    enddo
	CLOSE(1) !close file	
!make element array
    Allocate(elementSymbol(0:maxZ))	
	!Allocate(elementAtomicMass(0:maxZ))
!store Z and element symbol
    elementSymbol(0) = "NN" !Neutron
    Do i=1,size(isotopesInt,1)
	  elementSymbol(isotopesInt(i,1)) = isotopesChar(i)
	  !elementAtomicMass
	End Do
!close output file
    close(999) 
  End Subroutine readIsotopeData  
!================================================================================================================================================  
!================================================================================================================================================
!---------------------------------------------------------------------------------------------------
!
!    Input File
!
!---------------------------------------------------------------------------------------------------    
!read in input file
  Subroutine readInputFile()      
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, elementCounter, headerRow, fileRowCount
	Integer(kind=StandardInteger) :: replaceSymCount
	Character(len=3) :: bufferShortA, bufferShortB
	Character(len=32) :: buffera, bufferb, bufferc, bufferd, buffere, bufferf
	Character(len=64) :: fileRow
	Character(len=4) :: lastRow
	Character(len=255) :: bufferLongA
	Integer(kind=StandardInteger) :: bufferIA, pwscfFileCount
	Integer(kind=StandardInteger) :: forcePrepEAM
	Real(kind=DoubleReal) :: bufferDA	
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
    If(mpiProcessID.eq.0)Then
	  write(999,"(F8.4,A2,A24,A60)") ProgramTime(),"  ",&
	  "Reading user input file ",inputFileName
	End If	
!Allocate arrays	
	Allocate(unitVector(1:3,1:3))
	Allocate(dftInReplaceSym(1:20,1:2))	
!----------------------------------------------------------	
!Start set defaults, overridden by user input
!----------------------------------------------------------	
! Print out to terminal
    printToTerminal = 0	
!Force ZBL Core Settings	
	eamZBLPairLower = -2.1D20
	eamZBLPairUpper = -2.1D20
	eamZBLDensZero = -2.1D20
	eamZBLDensCutoff = -2.1D20
	eamZBLEmbeZero = -2.1D20
	eamZBLEmbeCutoff = -2.1D20
!Default unit vectors
	unitVector(1,1) = 1.0D0			!x1
	unitVector(2,1) = 0.0D0			!x2
	unitVector(3,1) = 0.0D0			!x3
	unitVector(1,2) = 0.0D0			!y1
	unitVector(2,2) = 1.0D0			!y2
	unitVector(3,2) = 0.0D0			!y3
	unitVector(1,3) = 0.0D0			!z1
	unitVector(2,3) = 0.0D0			!z2
	unitVector(3,3) = 1.0D0			!z3
!Run type and run type options
	calcRunType = 1
	optionReadConf = 1
	optionReadEAM = 1
	optionRunProcesses = 1
	optionRunPrep = 1
	optionRunDLPrep = 1
	optionRunPrepEAM = 1
	optionRunPWBatch = 1
!Potential type
	potType = 1
	eamInterpType = 1
	calcForcesOnOff = 1
	calcBMCalcOnOff = 0
	calcECCalcOnOff = 0
	calcEVCalcOnOff = 0
	dftInRadiusCutoff = 6.5D0
	Do i=1,size(dftInReplaceSym,1)
	  Do j=1,size(dftInReplaceSym,2)
	    dftInReplaceSym(i,j) = "   "
	  End Do
	End Do
	eamReducedPointsV = 50
	eamReducedPointsP = 50
	eamReducedPointsF = 25
	eamTrialPointsV = 1001
	eamTrialPointsP = 1001
	eamTrialPointsF = 1001
!Save files
    saveFileForces = "N"
    saveFileForces = "N"
    saveFileForces = "N"
	saveFileBM = "N"
!Optimisation weighting
	eamEnergyOptWeight = 100
	eamStressOptWeight = 10
	eamForceOptWeight = 1
	eamOptWeights =0
	eamOptWeightType = 1
	saCycles = 50
	saTStart = 3 
	saTEnd = 1 
	saTIncs = 3
	pwscfFilesList = "__empty__"
	eamPreparedFile = "prepared_eam.pot"
	dlpolyFileDir = "output/dlpoly"
!pwbatch
    pwbOutputDir = "pwb"	
	
!End set defaults
!Count rows in file
	fileRowCount = 0
  	Open(UNIT=1,FILE=inputFileName) 
	If(mpiProcessID.eq.0)Then
  	  Open(UNIT=2,FILE="temp/tmpinput.in") 
	End If
    Do i=1,maxFileRows 
      fileRowCount = fileRowCount + 1
	  Read(1,"(A64)",IOSTAT=ios) fileRow
	  If (ios /= 0) Then
	    EXIT 
	  End If
	  If(mpiProcessID.eq.0)Then
	    fileRow = TrimSpaces(fileRow)
	    If(fileRow(1:1).eq."!".or.fileRow(1:1).eq." ")Then	  
	    Else
	      If(fileRow(1:1).eq."#")Then
		    write(2,"(A)") StrToUpper(fileRow)
		  Else
	        write(2,"(A)") fileRow
		  End If  
		End If 
	  End If
	End Do
	If(mpiProcessID.eq.0)Then
	  Close(2)
	End If
	Close(1)
!Synch mpi processes
	Call MPI_synchProcesses()
!Set other variables
    forcePrepEAM = 0	
!open & read in file	
  	Open(UNIT=1,FILE="temp/tmpinput.in") 
	pwscfFileCount = 0
    do i=1,fileRowCount 
!Read in line
	  Read(1,"(A64)",IOSTAT=ios) fileRow
!-----------------------------
! Print to terminal
!-----------------------------
	  If(StrToUpper(fileRow(1:6)).eq."#PRINT")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		If(StrToUpper(fileRow(1:1)).eq."Y")Then
		  printToTerminal = 1
	    End If
		If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		  print *,ProgramTime(),"Reading user input file"
		End If
	  End If
!-----------------------------
! Run type
!-----------------------------	  
	  If(StrToUpper(fileRow(1:8)).eq."#RUNTYPE")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		fileRow = StrToUpper(fileRow)
		If(StrToUpper(fileRow(1:4)).eq."ENER")then		!ENERGY
		  calcRunType = 1
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."BMOD")then		!BULKMODULUS
		  calcRunType = 3
		  calcBMCalcOnOff = 1
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."ECON")then		!ELASTICCONSTANTS
		  calcRunType = 4
		  calcBMCalcOnOff = 1
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."EVAL")then		!EVALUATE
		  calcRunType = 5
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."OPTI")then		!OPTIMISE
		  calcRunType = 6
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."OPTF")then		!OPTIMISE Full Eval
		  calcRunType = 7
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."EVAF")then		!EVALUATE
		  calcRunType = 8
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 0
		  optionRunProcesses = 1
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
!Potential only, no calculations
		if(StrToUpper(fileRow(1:4)).eq."PRP1")then		!PREP EAM File - read in, and output a formatted eam file
		  calcRunType = 15
!Set module/subroutine options
		  optionReadConf = 0
		  optionReadEAM = 1
		  optionRunPrep = 0
		  optionRunDLPrep = 0
		  optionRunProcesses = 0
		  optionRunPrepEAM = 1
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."PRP2")then		!PREP EAM File - read in, and output a formatted eam file
		  calcRunType = 16
!Set module/subroutine options
		  optionReadConf = 0
		  optionReadEAM = 1
		  optionRunPrep = 0
		  optionRunDLPrep = 0
		  optionRunProcesses = 0
		  optionRunPrepEAM = 1
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."PRP3")then		!PREP EAM File - read in, and output a formatted eam file
		  calcRunType = 17
!Set module/subroutine options
		  optionReadConf = 0
		  optionReadEAM = 1
		  optionRunPrep = 0
		  optionRunDLPrep = 0
		  optionRunProcesses = 0
		  optionRunPrepEAM = 1
		  optionRunPWBatch = 0
	    End If
		if(StrToUpper(fileRow(1:4)).eq."PRP4")then		!PREP EAM File, DLPOLY eam and conf file
		  calcRunType = 18
!Set module/subroutine options
		  optionReadConf = 1
		  optionReadEAM = 1
		  optionRunPrep = 1
		  optionRunDLPrep = 1
		  optionRunProcesses = 0
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 0
	    End If
		If(StrToUpper(fileRow(1:4)).eq."PWB1")then		!PWscf Batch Files 1 - input files
		  calcRunType = 31
!Set module/subroutine options
		  optionReadConf = 0
		  optionReadEAM = 0
		  optionRunPrep = 0
		  optionRunDLPrep = 0
		  optionRunProcesses = 0
		  optionRunPrepEAM = 0
		  optionRunPWBatch = 1
		End If
		if(StrToUpper(fileRow(1:4)).eq."TEST")then		!Test
		  calcRunType = 99
	    End If		
	    If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
		  print *,ProgramTime(),"Run type ",calcRunType," ",fileRow(1:4)
		End If
	  End If	  
!-----------------------------
! Potential details
!-----------------------------	  
	  if(StrToUpper(fileRow(1:10)).eq."#POTENTIAL")then
!read next line
	    Read(1,*,IOSTAT=ios) potentialFilePath
	  endif
	  if(StrToUpper(fileRow(1:8)).eq."#POTTYPE")then
!read next line
	    Read(1,*,IOSTAT=ios) bufferb
		If(StrToUpper(bufferb(1:2)).eq."NO")Then
		  potType = 1
		End If
		If(StrToUpper(bufferb(1:2)).eq."LA")Then
		  potType = 2
		End If	
		If(StrToUpper(bufferb(1:2)).eq."DL")Then
		  potType = 3
		End If		
	  endif
	  If(StrToUpper(fileRow(1:8)).eq."#ZBLCORE")then
!read next line
        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere, bufferf
		Read(buffera,*) eamZBLPairLower
		Read(bufferb,*) eamZBLPairUpper
		Read(bufferc,*) eamZBLDensZero
		Read(bufferd,*) eamZBLDensCutoff
		Read(buffere,*) eamZBLEmbeZero
		Read(bufferf,*) eamZBLEmbeCutoff
	  End If	  
	  
	  
	  
	  if(StrToUpper(fileRow(1:15)).eq."#CONFIGURATIONS")then
!read next line
	    Read(1,*,IOSTAT=ios) configurationsFilePath
	  endif
	  if(StrToUpper(fileRow(1:11)).eq."#SAVECOORDS")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFileCoords = StrToUpper(buffera(1:1))
	  endif
	  if(StrToUpper(fileRow(1:18)).eq."#SAVENEIGHBOURLIST")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFileNeighbourList = StrToUpper(buffera(1:1))
	  endif
	  if(StrToUpper(fileRow(1:8)).eq."#SAVEPOT")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFilePot = StrToUpper(buffera(1:1))
	  endif
	  if(StrToUpper(fileRow(1:11)).eq."#SAVEFORCES")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFileForces = StrToUpper(buffera(1:1))
	  endif
	  if(StrToUpper(fileRow(1:11)).eq."#SAVEBM")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFileBM = StrToUpper(buffera(1:1))
	  endif	  
	  
	  

	  
!Unit vector components
	  if(StrToUpper(fileRow(1:11)).eq."#UNITVECTOR")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		Read(buffera,*) unitVector(1,1)			!x1
		Read(bufferb,*) unitVector(2,1)			!x2
		Read(bufferc,*) unitVector(3,1)			!x3
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		Read(buffera,*) unitVector(1,2)			!y1
		Read(bufferb,*) unitVector(2,2)			!y2 
		Read(bufferc,*) unitVector(3,2)			!y3 
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		Read(buffera,*) unitVector(1,3)			!z1
		Read(bufferb,*) unitVector(2,3)			!z2 
		Read(bufferc,*) unitVector(3,3)			!z3 
	  endif
	  
!EAM Interpolation Types
	  if(StrToUpper(fileRow(1:10)).eq."#EAMINTERP")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)
		if(buffera(1:4).eq."NONE")then
		  eamInterpType = 0
		endif
		if(buffera(1:6).eq."LINEAR")then
		  eamInterpType = 1
		endif
		if(buffera(1:10).eq."THREEPOINT")then
		  eamInterpType = 2
		endif
		if(buffera(1:9).eq."FOURPOINT")then
		  eamInterpType = 3
		endif
		if(buffera(1:9).eq."FIVEPOINT")then
		  eamInterpType = 4
		endif
	  endif
	  
!#eamreducedpoints	  
	  if(StrToUpper(fileRow(1:17)).eq."#EAMREDUCEDPOINTS")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		Read(buffera,*) eamReducedPointsV
		Read(bufferb,*) eamReducedPointsP
		Read(bufferc,*) eamReducedPointsF
	  endif
	  
!#eamtrialpoints	  
	  if(StrToUpper(fileRow(1:15)).eq."#EAMTRIALPOINTS")then
!read next line
	    !Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		!Read(buffera,*) eamTrialPointsV
		!Read(bufferb,*) eamTrialPointsP
		!Read(bufferc,*) eamTrialPointsF
	  endif	 

!#weighting	  
	  if(StrToUpper(fileRow(1:10)).eq."#WEIGHTING")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		Read(buffera,*) eamOptWeights(1) !Energy
		Read(bufferb,*) eamOptWeights(2) !Stress
		Read(bufferc,*) eamOptWeights(3) !Force
		Read(bufferd,*) eamOptWeights(4) !Eq Vol
	  endif	 
!#weighting	type  
	  If(StrToUpper(fileRow(1:11)).eq."#WEIGHTTYPE")then	  
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		If(StrToUpper(buffera(1:3)).eq."ABS")Then	
		  eamOptWeightType = 1		!Absolute value, difference between ref and calculated
		End If
		If(StrToUpper(buffera(1:3)).eq."REL")Then
		  eamOptWeightType = 2		!Relative value, fraction difference between ref and calculated
		End If
	  End If

!#saCycles	  simulated annealing cycles	  
	  If(StrToUpper(fileRow(1:9)).eq."#SACYCLES")Then
	  	!read next line
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd
		Read(buffera,*) saCycles
		Read(bufferb,*) saTStart
		Read(bufferc,*) saTEnd
		Read(bufferd,*) saTIncs		
      End If
	  
!#saCycles	  simulated annealing cycles	  
	  If(StrToUpper(fileRow(1:7)).eq."#SATEMP")Then
	  	!read next line
	    Read(1,*,IOSTAT=ios) buffera, bufferb
		Read(buffera,*) saTemp
		Read(bufferb,*) saSpreadFactor
      End If	  
	  

		
!Calc Forces
	  If(StrToUpper(fileRow(1:11)).eq."#CALCFORCES")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		fileRow = StrToUpper(fileRow)	  
		If(fileRow(1:1).eq."Y".or.fileRow(1:1).eq."1")Then
	      calcForcesOnOff = 1
		Else
		  calcForcesOnOff = 0
		End If
	  End If 
	  
!Calc Stress
	  If(StrToUpper(fileRow(1:11)).eq."#CALCSTRESS")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)	  
		If(buffera(1:1).eq."Y".or.buffera(1:1).eq."1")Then
	      calcStressOnOff = 1
		Else
		  calcStressOnOff = 0
		End If
	  End If 
	  
!Calc All BM
	  If(StrToUpper(fileRow(1:10)).eq."#CALCALLBM")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)	  
		If(buffera(1:1).eq."Y".or.buffera(1:1).eq."1")Then
	      calcBMCalcOnOff = 1
		Else
		  calcBMCalcOnOff = 0
		End If
	  End If 
	  
!Calc All EC
	  If(StrToUpper(fileRow(1:10)).eq."#CALCALLEC")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)	  
		If(buffera(1:1).eq."Y".or.buffera(1:1).eq."1")Then
	      calcECCalcOnOff = 1
		Else
		  calcECCalcOnOff = 0
		End If
	  End If 	

	  
!Calc All Equilibrium Volumes
	  If(StrToUpper(fileRow(1:10)).eq."#CALCALLEV")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)	  
		If(buffera(1:1).eq."Y".or.buffera(1:1).eq."1")Then
	      calcEVCalcOnOff = 1
		Else
		  calcEVCalcOnOff = 0
		End If
	  End If 	  
	  
	  
	  	  
!Dft Input Settings
	  If(StrToUpper(fileRow(1:9)).eq."#DFTIN-RC")then
	    Read(1,*,IOSTAT=ios) buffera
	    Read(buffera,*) dftInRadiusCutoff
	  End If
	  If(StrToUpper(fileRow(1:17)).eq."#DFTIN-REPLACESYM")then
	    Read(1,*,IOSTAT=ios) buffera
		Read(buffera,*) replaceSymCount
		Do j=1,replaceSymCount
	      Read(1,*,IOSTAT=ios) bufferShortA, bufferShortB
          dftInReplaceSym(j,1) = TrimSpaces(bufferShortA)
          dftInReplaceSym(j,2) = TrimSpaces(bufferShortB)
		End Do  
	  End If
	  If(StrToUpper(fileRow(1:16)).eq."#DFTIN-OPTENERGY")then
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
	    Read(buffera,*) bufferDA
	    Read(bufferc,*) bufferIA
		dftInOptEnergy = UnitConvert(bufferDA,TrimSpaces(bufferb),"EV")
		dftInOptEnergy = dftInOptEnergy / (1.0D0 * bufferIA)
	  End If
	  If(StrToUpper(fileRow(1:16)).eq."#DFTIN-COHENERGY")then
	    Read(1,*,IOSTAT=ios) buffera, bufferb
	    Read(buffera,*) bufferDA
		dftInCohEnergy = UnitConvert(bufferDA,TrimSpaces(bufferb),"EV")
	  End If
	  If(StrToUpper(fileRow(1:15)).eq."#DFT-PWSCFFILES")then
	    Do j=1,1000
	      Read(1,*,IOSTAT=ios) buffera	!read next line
	      If(buffera(1:3).eq."   ".or.buffera(1:1).eq."#".or.ios/=0)Then
		    Backspace(1)
		    Exit
		  Else	
			pwscfFileCount = pwscfFileCount + 1
			pwscfFilesList(pwscfFileCount) = TrimSpaces(buffera)
		  End If		
	    End Do
	  End If	  
!-----------------------------
! EAM Prep Only Inputs
!-----------------------------	  
	  If(StrToUpper(fileRow(1:12)).eq."#EAMPREPFILE")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera		
		eamPreparedFile = trim(buffera)
		forcePrepEAM = 1
	  End If
!-----------------------------
! Make PWscf Batch Files Only Options
!-----------------------------	  
	  If(StrToUpper(fileRow(1:11)).eq."#PWBOPTCONF")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera		
		pwbOptConfFile = trim(buffera)
	  End If	
	  If(StrToUpper(fileRow(1:9)).eq."#BATCHDIR")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera		
		pwbOutputDir = trim(buffera)
	  End If	  
    End Do
!-----------------------------
! Force some options, depending on input
!-----------------------------
	If(forcePrepEAM.eq.1)Then
	  optionRunPrepEAM = 1
	End If
!Fill in any other variables	
	globalUnitVector = unitVector
!close file	
	CLOSE(1) 	  
!close output file
    close(999)	   
  End Subroutine readInputFile
!================================================================================================================================================  
!================================================================================================================================================
!---------------------------------------------------------------------------------------------------
!
!    Reading Potential File
!
!---------------------------------------------------------------------------------------------------     
  Subroutine readEamPotPrep()  
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, n  
	Integer(kind=StandardInteger) :: rowCount
	Character(len=255) :: fileRow
	Character(len=128) :: bufferA, bufferB, bufferC, bufferD, bufferE, bufferF
	Integer(kind=StandardInteger) :: bufferIntA, bufferIntB
	Real(kind=DoubleReal) :: bufferDoubleA, bufferDoubleB
	Integer(kind=StandardInteger) :: elementTypeCount, elementZ
	Integer(kind=StandardInteger) :: nrho, nr, potentialLinesRho, potentialLinesR
	Real(kind=DoubleReal) :: drho, dr, cutoff
	Real(kind=DoubleReal) :: radius, rho, tempDouble
	Character(len=2), Dimension( : ), Allocatable ::  functionElements	
!---------------------------------------------
! Set potentialFilePathTemp on all processes
!---------------------------------------------
	If(potType.eq.1)Then
	  !Normal type, use file as provided
	  potentialFilePathTemp = potentialFilePath
    End If
	If(potType.eq.2)Then
	  !LAMMPS type
	  potentialFilePathTemp = trim(potentialFilePath)//".temp"
	End If
	If(potType.eq.3)Then
	  !DLPOLY type
	  potentialFilePathTemp = trim(potentialFilePath)//".temp"
	End If
	If(mpiProcessID.eq.0)Then
!--------------------  
! NORMAL
!--------------------
    If(potType.eq.1)Then
	  !Normal type, use file as provided
    End If
!--------------------  
! LAMMPS
!--------------------	
	If(potType.eq.2)Then
	  !LAMMPS type
!open output potential file
	  Open(UNIT=2,FILE=trim(potentialFilePathTemp)) 
!open LAMMPS input potential file
	  Open(UNIT=1,FILE=trim(potentialFilePath)) 
	  rowCount = 0
      Do i=1,maxFileRows 
!Read in line
	    Read(1,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
	    If (ios /= 0) then
	      EXIT 
	    End If		
	    If(fileRow(1:1).ne."#")Then		  
		  rowCount = rowCount + 1
		  If(rowCount.eq.1)Then
		    Read(fileRow,*) elementTypeCount	
			Allocate(functionElements(1:elementTypeCount))
		  End If
		  If(rowCount.eq.2)Then
		    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE
		    Read(bufferA,*) nrho	!number of points at which p(r) is evaluated
		    Read(bufferB,*) drho	
		    Read(bufferC,*) nr	    !number of points at which V(r) and F(p) are evaluated
		    Read(bufferD,*) dr	
		    Read(bufferE,*) cutoff		    
			potentialLinesR = Ceil(1.0D0*(nrho/5))
			potentialLinesRho = Ceil((1.0D0*nrho)/5.0D0)
		  End If
		  If(rowCount.eq.3)Then
!Read in embedding function and density function for each element
		    Do j=1,elementTypeCount
			  If(j.gt.1)Then
			    !Read in next element row
                Read(1,"(A255)",IOSTAT=ios) fileRow	
			  End If
!Read element type	       	 
			  Read(fileRow,*) elementZ
!Embedding function
			  write(2,"(A8,A2)") "EMBE    ",elementSymbol(elementZ)
			  functionElements(j) = elementSymbol(elementZ)
			  rho = 0.0D0
!Embedding function
			  Do k=1,potentialLinesRho
			    Read(1,"(A255)",IOSTAT=ios) fileRow
				If(k.lt.potentialLinesRho.or.mod(nrho,5).eq.0)Then
			      Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE	 
				  Read(bufferA,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  rho = rho + drho 
				  Read(bufferB,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  rho = rho + drho
				  Read(bufferC,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  rho = rho + drho
				  Read(bufferD,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  rho = rho + drho
				  Read(bufferE,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  rho = rho + drho
				Else  
				  If(mod(nrho,5).eq.1)Then
				    Read(fileRow,*) bufferA	 
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  End If 
				  If(mod(nrho,5).eq.2)Then
				    Read(fileRow,*) bufferA, bufferB
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				    rho = rho + drho 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  End If
				  If(mod(nrho,5).eq.3)Then
				    Read(fileRow,*) bufferA, bufferB, bufferC
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				    rho = rho + drho 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				    rho = rho + drho 
				    Read(bufferC,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  End If
				  If(mod(nrho,5).eq.3)Then
				    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				    rho = rho + drho 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				    rho = rho + drho 
				    Read(bufferC,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				    rho = rho + drho 
				    Read(bufferD,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") rho,"  ",tempDouble
				  End If				  
				End If
			  End Do
!Density function
              write(2,"(A8,A2)") "DENS    ",elementSymbol(elementZ)
			  radius = 0.0D0  
!Density function
			  Do k=1,potentialLinesR+1
			    Read(1,"(A255)",IOSTAT=ios) fileRow
				If(k.lt.(potentialLinesR+1).or.(k.eq.(potentialLinesR+1).and.mod(nr,5).eq.0))Then
			      Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE	 
				  Read(bufferA,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  radius = radius + dr 
				  Read(bufferB,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  radius = radius + dr
				  Read(bufferC,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  radius = radius + dr
				  Read(bufferD,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  radius = radius + dr
				  Read(bufferE,*) tempDouble 
				  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  radius = radius + dr
				Else  
				  If(mod(nr,5).eq.1)Then
				    Read(fileRow,*) bufferA	 
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  End If 
				  If(mod(nr,5).eq.2)Then
				    Read(fileRow,*) bufferA, bufferB
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				    radius = radius + dr 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  End If
				  If(mod(nr,5).eq.3)Then
				    Read(fileRow,*) bufferA, bufferB, bufferC
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				    radius = radius + dr 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				    radius = radius + dr 
				    Read(bufferC,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  End If
				  If(mod(nr,5).eq.3)Then
				    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
				    Read(bufferA,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				    radius = radius + dr 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				    radius = radius + dr 
				    Read(bufferC,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				    radius = radius + dr 
				    Read(bufferD,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",tempDouble
				  End If				  
				End If
			  End Do
			End Do !Element count			
!Pair potentials
            Do j=1,elementTypeCount
			  Do n=j,elementTypeCount
			    write(2,"(A8,A2,A4,A2)") "PAIR    ",&
				functionElements(j),"    ",functionElements(n)
				radius = 0.0D0
!Pair function
			    Do k=1,potentialLinesR+1
			      Read(1,"(A255)",IOSTAT=ios) fileRow
				  If(k.lt.(potentialLinesR+1).or.&
				  (k.eq.(potentialLinesR+1).and.mod(nr,5).eq.0))Then
			        Read(fileRow,*) bufferA, bufferB, bufferC, bufferD, bufferE	 
				    Read(bufferA,*) tempDouble 
					If(radius.eq.0.0D0)Then
					  write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",0.0D0
					Else
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
					End If
				    radius = radius + dr 
				    Read(bufferB,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    radius = radius + dr
				    Read(bufferC,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    radius = radius + dr
				    Read(bufferD,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    radius = radius + dr
				    Read(bufferE,*) tempDouble 
				    write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    radius = radius + dr
				  Else  
				    If(mod(nr,5).eq.1)Then
				      Read(fileRow,*) bufferA	 
				      Read(bufferA,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    End If 
				    If(mod(nr,5).eq.2)Then
				      Read(fileRow,*) bufferA, bufferB
				      Read(bufferA,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				      radius = radius + dr 
				      Read(bufferB,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    End If
				    If(mod(nr,5).eq.3)Then
				      Read(fileRow,*) bufferA, bufferB, bufferC
				      Read(bufferA,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				      radius = radius + dr 
				      Read(bufferB,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				      radius = radius + dr 
				      Read(bufferC,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    End If
				    If(mod(nr,5).eq.3)Then
				      Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
				      Read(bufferA,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				      radius = radius + dr 
				      Read(bufferB,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				      radius = radius + dr 
				      Read(bufferC,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				      radius = radius + dr 
				      Read(bufferD,*) tempDouble 
				      write(2,"(E24.16E3,A2,E24.16E3)") radius,"  ",(tempDouble/(1.0D0*radius))
				    End If				  
				  End If
			    End Do
			  End Do
            End Do
		  End If
		End If		
	  End Do
    End If
!--------------------  
! DLPOLY
!--------------------	
	If(potType.eq.3)Then
	  !DLPOLY type
!open output potential file
	  Open(UNIT=2,FILE=trim(potentialFilePathTemp)) 
!open DLPOLY input potential file
	  Open(UNIT=1,FILE=trim(potentialFilePath)) 
	  rowCount = 0
      Do i=1,maxFileRows 
!Read in line
	    Read(1,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
	    If (ios /= 0) then
	      EXIT 
	    End If		
	    If(fileRow(1:1).ne."#")Then		  
		  rowCount = rowCount + 1
		End If
      End Do
    End If	  
	End If  !If mpiprocessid=0  
    Close(1)
	Close(2)
!-----------------------------
! Deallocate arrays
!-----------------------------
	If(Allocated(functionElements))Then
	  Deallocate(functionElements)
	End If	
  End Subroutine readEamPotPrep
!================================================================================================================================================  
!================================================================================================================================================    
!read eam.pot subroutine  
  Subroutine readEamPot()	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, functionResult, fileRows
	Integer(kind=StandardInteger) :: potentialCounter, elementCounter, dataCounter
	Integer(kind=StandardInteger) :: atomA, atomB, potType, atomMax, atomMin
	Integer(kind=StandardInteger) :: potKey
	Integer(kind=StandardInteger) :: potDataStart, potDataLength
	Character(len=2), Dimension( : ), Allocatable :: elementsTemp 
	Character(len=4) :: potTypeText
	Character(len=32) :: buffera, bufferb, bufferc, bufferd
	Character(len=255) :: bufferLongA
	Character(len=255) :: fileRowData
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!open output potfile
    if(saveFilePot.eq."Y".and.mpiProcessID.eq.0)then
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.pot"
	  open(unit=21,file=trim(outputFile))
	endif
!write to output file
    If(mpiProcessID.eq.0)Then
	  write(999,"(F8.4,A2,A27,A60)") ProgramTime(),"  ",&
	  "Reading EAM potential from ",potentialFilePathTemp
	  !"Reading EAM potential from ",potentialFilePath
	endif	
!Deallocate already allocated arrays
    If(Allocated(elements))Then
	  Deallocate(elements)
	End If
    If(Allocated(elementsTemp))Then
	  Deallocate(elementsTemp)
	End If
    If(Allocated(eamKey))Then
	  Deallocate(eamKey)
	End If
    If(Allocated(eamData))Then
	  Deallocate(eamData)
	End If	
!Set EAM type
    eamType = 1		!1 = default type/standard   2 = 2BMEAM	
!allocate elements array
    Allocate(elements(1:300))
    Allocate(elementsTemp(1:300))
!"blank" elements array
	do i=1,size(elements)
	  elements(i) = "ZZ"
	enddo	
!Count file rows
	fileRows = 0
  	Open(UNIT=1,FILE=potentialFilePathTemp) 
    do i=1,maxFileRows 
!count file rows
	  fileRows = fileRows + 1
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
	enddo
!close file	
	CLOSE(1)  		
!Read eam pot file and store unique elements
	potentialCounter = 0
  	Open(UNIT=1,FILE=potentialFilePathTemp) 
    do i=1,fileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
!skip blank or commented line      
	  if(buffera(1:1).eq." ".or.buffera(1:1).eq."#".or.buffera(1:1).eq."!")then
!skip
	  else 
!check if PAIR DENS or EMBE
        if(StrToUpper(buffera(1:4)).eq."PAIR")then
	      BACKSPACE(1)
		  Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		  potentialCounter = potentialCounter + 1
		  functionResult = AddUniqueElement(bufferb(1:2))
		  functionResult = AddUniqueElement(bufferc(1:2))
	    elseif(&
	    StrToUpper(buffera(1:4)).eq."DENS".or.StrToUpper(buffera(1:4)).eq."EMBE"&
	    .or.StrToUpper(buffera(1:4)).eq."EMBS"&
	    .or.StrToUpper(buffera(1:4)).eq."DEND".or.StrToUpper(buffera(1:4)).eq."EMBD"&
	    )then	  
!adjust potential type 
	      if(&
	        StrToUpper(buffera(1:4)).eq."DEND".or.StrToUpper(buffera(1:4)).eq."EMBD"&
		    .or.StrToUpper(buffera(1:4)).eq."EMBS"&
	      )then	
		    eamType = 2	!2BMEAM
		  endif
!re-read line, count potential and add element if unique
	      BACKSPACE(1)
		  Read(1,*,IOSTAT=ios) buffera, bufferb
		  potentialCounter = potentialCounter + 1
		  functionResult = AddUniqueElement(bufferb(1:2))
	    endif
	  endif
	enddo
!close file	
	CLOSE(1)   	  
!adjust elements array
    elementCounter = 0
	do i=1,size(elements)
	  if(elements(i).eq."ZZ")then
	    exit
	  else 	  
	    elementCounter = elementCounter + 1
		elementsTemp(i) = elements(i)
	  endif
	enddo
	Deallocate(elements)
	Allocate(elements(1:elementCounter))
	do i=1,elementCounter
	  elements(i) = elementsTemp(i)
	enddo
	Deallocate(elementsTemp)
	elementCount = elementCounter	
!store/output elements	
    If(mpiProcessID.eq.0)Then
      write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	  Do i=1,size(elements)
	    write(999,"(A8,I4,A2,A2)") "Element ",i,": ",elements(i)
	  End Do
	End If			
!count data points
	dataCounter = 0	
Open(UNIT=1,FILE=potentialFilePathTemp) 
    do i=1,fileRows 
!Read in line
	  Read(1,'(A64)',IOSTAT=ios) fileRowData
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  endif
!skip blank or commented line   
	  if(fileRowData(1:8).ne."        ".and.fileRowData(1:1).ne."#".and.fileRowData(1:1).ne."!")then
!check if PAIR DENS or EMBE - store pot data start/length	
		read(fileRowData,*) potTypeText 
		potTypeText = trim(potTypeText)
	    if(StrToUpper(potTypeText(1:4)).eq."PAIR".or.&
	    StrToUpper(potTypeText(1:4)).eq."DENS".or.&
	    StrToUpper(potTypeText(1:4)).eq."SDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."DDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."EMBE".or.&
	    StrToUpper(potTypeText(1:4)).eq."SEMB".or.&
	    StrToUpper(potTypeText(1:4)).eq."DEMB")then
!potential header line          
		else		
	      dataCounter = dataCounter + 1
		endif
      endif
	enddo
	CLOSE(1) 
!potential counts
    if(eamType.eq.1)then
	  numberPotentials = (elementCounter * ( elementCounter + 5)) / 2
	  pairCount = (elementCounter * ( elementCounter + 1)) / 2
	  densCount = elementCounter
	  embeCount = elementCounter
	elseif(eamType.eq.2)then
	  if(elementCounter.eq.1)then
	    numberPotentials = 1 + (elementCounter * ( elementCounter + 7)) / 2
		pairCount = (elementCounter * ( elementCounter + 1)) / 2
	    densCount = 1
		dendCount = elementCounter
	    embeCount = elementCounter
	    embdCount = elementCounter
	  else
	    numberPotentials = (elementCounter * ( elementCounter + 3))
		pairCount = (elementCounter * ( elementCounter + 1)) / 2
	    densCount = (elementCounter * ( elementCounter - 1)) / 2
		dendCount = elementCounter
	    embeCount = elementCounter
	    embdCount = elementCounter
	  endif
	endif	
!store/output potential count etc
    If(mpiProcessID.eq.0)Then
    write(999,"(A6,A16,F8.4)") "      ","Program time:   ",ProgramTime()
	write(999,"(A6,A28,I4)") "      ","Element count:              ",elementCounter
	write(999,"(A6,A28,I4)") "      ","EAM Type:                   ",eamType
	write(999,"(A6,A28,I4)") "      ","Potentials count:           ",numberPotentials
	write(999,"(A6,A28,I4)") "      ","Potentials expected:        ",potentialCounter
    if(eamType.eq.1)then
	  write(999,"(A6,A28,I4)") "      ","Pairs functions:            ",pairCount
	  write(999,"(A6,A28,I4)") "      ","Density functions:          ",densCount
	  write(999,"(A6,A28,I4)") "      ","Embedding functions:        ",embeCount	
	elseif(eamType.eq.2)then
	  write(999,"(A6,A28,I4)") "      ","Pairs functions:            ",pairCount
	  write(999,"(A6,A28,I4)") "      ","S-band Density functions:   ",densCount
	  write(999,"(A6,A28,I4)") "      ","D-band Density functions:   ",dendCount
	  write(999,"(A6,A28,I4)") "      ","S-band Embedding functions: ",embeCount	
	  write(999,"(A6,A28,I4)") "      ","S-band Embedding functions: ",embdCount	
	endif	
	write(999,"(A6,A28,I8)") "      ","Data points:                ",dataCounter
	write(999,"(A6,A28,A1)") "      ","Save coords file:           ",saveFileCoords
	write(999,"(A6,A28,A1)") "      ","Save neighbour list file:   ",saveFileNeighbourList
	write(999,"(A6,A28,A1)") "      ","Save potential file:        ",saveFilePot
	write(999,"(A6,A28,I8)") "      ","Calculation run type:       ",calcRunType
	write(999,"(A6,A28)") "      ","Unit vector:                "
	write(999,"(A6,F8.4,F8.4,F8.4)") "      ",unitVector(1,1),unitVector(2,1),unitVector(3,1)
	write(999,"(A6,F8.4,F8.4,F8.4)") "      ",unitVector(1,2),unitVector(2,2),unitVector(3,2)
	write(999,"(A6,F8.4,F8.4,F8.4)") "      ",unitVector(1,3),unitVector(2,3),unitVector(3,3)
	End If		
!Allocate eam keys/data
    Allocate(eamKey(1:numberPotentials,1:5))	
    Allocate(eamData(1:dataCounter,1:3))
!Initialise data array
    Do i=1,size(eamData,1)
      Do j=1,size(eamData,2)
	    eamData(i,j) = 0.0D0
	  End Do
	End Do
	potentialCounter = 0
	eamType = 1
	dataCounter = 0
	potDataStart = 1
	potDataLength = 0
	Open(UNIT=1,FILE=potentialFilePathTemp) 
    do i=1,fileRows 
!Read in line
	  Read(1,'(A64)',IOSTAT=ios) fileRowData
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  endif
!skip blank or commented line   
	  If(fileRowData(1:8).ne."        ".and.fileRowData(1:1)&
	    .ne."#".and.fileRowData(1:1).ne."!")Then
!check if PAIR DENS or EMBE - store pot data start/length	
		read(fileRowData,*) potTypeText 
		potTypeText = trim(potTypeText)
	    if(StrToUpper(potTypeText(1:4)).eq."PAIR".or.&
	    StrToUpper(potTypeText(1:4)).eq."DENS".or.&
	    StrToUpper(potTypeText(1:4)).eq."SDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."DDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."EMBE".or.&
	    StrToUpper(potTypeText(1:4)).eq."SEMB".or.&
	    StrToUpper(potTypeText(1:4)).eq."DEMB")then
!read potential header line
          potentialCounter = potentialCounter + 1
	      if(StrToUpper(potTypeText(1:4)).eq."PAIR".or.&
	      StrToUpper(potTypeText(1:4)).eq."SDEN")then
		    Read(fileRowData,*) buffera, bufferb, bufferc
		  else
		    Read(fileRowData,*) buffera, bufferb
		  endif		  
!save start/length details of last potential         
          if(potentialCounter.gt.1)then
	        eamKey(potKey,4) = potDataStart
	        eamKey(potKey,5) = potDataLength
		  endif
!get potential key
          potKey = 0
		  potType = 0 
		  if(StrToUpper(potTypeText(1:4)).eq."PAIR")then	
		    potType = 1	    
		    atomA = QueryUniqueElement(bufferb(1:2))
		    atomB = QueryUniqueElement(bufferc(1:2))
            atomMax = max(atomA,atomB)-1
            atomMin = min(atomA,atomB)-1
		    potKey = 1+atomMin+(atomMax*(atomMax+1))/2
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."DENS")then	
			potType = 2	    
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."SDEN")then	
			potType = 2
            eamType = 2		  
		    atomA = QueryUniqueElement(bufferb(1:2))
		    atomB = QueryUniqueElement(bufferc(1:2))
			if(elementCounter.eq.1)then
              potKey = pairCount + 1
            else
              atomMax = max(atomA,atomB)-1
              atomMin = min(atomA,atomB)-1
              potKey = pairCount+1+atomMin+(atomMax*(atomMax-1))/2
            endif
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."DDEN")then	
			potType = 4
            eamType = 2		  
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
			potKey = pairCount + densCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."EMBE")then	
			potType = 3	    
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + densCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."SEMB")then	
			potType = 3	    
            eamType = 2
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + densCount + dendCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."DEMB")then	
			potType = 5	    
            eamType = 2
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + densCount + dendCount + embdCount + atomA
		  endif
!store potential details
          eamKey(potKey,1) = atomA
          eamKey(potKey,2) = atomB
          eamKey(potKey,3) = potType
!reset start/length details
		  potDataStart = potDataStart + potDataLength
		  potDataLength = 0
		else		
	      dataCounter = dataCounter + 1
		  potDataLength = potDataLength + 1
		  read(fileRowData,*) buffera, bufferb
		  read(buffera,*) eamData(dataCounter,1)
		  read(bufferb,*) eamData(dataCounter,2)		  	  
		endif
      endif
	enddo
!save last start/length entry
	eamKey(potKey,4) = potDataStart
	eamKey(potKey,5) = potDataLength
	CLOSE(1) 	
!Fill in elementsCharge array
    If(Allocated(elementsCharge))Then
	  Deallocate(elementsCharge)
	End If
    Allocate(elementsCharge(1:size(elements)))
    Do i=1,size(elements)
	  Do j=0,size(elementSymbol,1)-1
	    If(elements(i).eq.elementSymbol(j))Then
		  elementsCharge(i) = j
	    End If
	  End Do
	End Do
!output/store potentials	
    If(mpiProcessID.eq.0)Then
      write(999,"(A6,A32)") "      ","Summary of potential functions "
      do i=1,size(eamKey,1)
	    If(eamKey(i,2).gt.0)Then
	      write(999,"(A6,I8,A4,A4,I8,I8,I8,I8,I8)") "      ",i,&
		  elements(eamKey(i,1)),elements(eamKey(i,2)),&
		  eamKey(i,1),eamKey(i,2),&
	      eamKey(i,3),eamKey(i,4),eamKey(i,5)
		Else
		  write(999,"(A6,I8,A4,A4,I8,I8,I8,I8,I8)") "      ",i,&
		  elements(eamKey(i,1)),"    ",&
		  eamKey(i,1),eamKey(i,2),&
	      eamKey(i,3),eamKey(i,4),eamKey(i,5)
		End If
      enddo
	End If
!close output file
    close(999)	
!-----------------------------
! Deallocate arrays
!-----------------------------
	If(Allocated(elementsTemp))Then
	  Deallocate(elementsTemp)
	End If			
  End Subroutine readEamPot
!================================================================================================================================================  
!================================================================================================================================================	
!re-read eam potential file, assume element list is the same 
  Subroutine reReadEamPot(inputEamFile)	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, functionResult, fileRows
	Integer(kind=StandardInteger) :: potentialCounter, elementCounter, dataCounter
	Integer(kind=StandardInteger) :: atomA, atomB, potType, atomMax, atomMin
	Integer(kind=StandardInteger) :: potKey
	Integer(kind=StandardInteger) :: potDataStart, potDataLength
	Character(len=2), Dimension( : ), Allocatable :: elementsTemp 
	Character(len=4) :: potTypeText
	Character(len=32) :: buffera, bufferb, bufferc, bufferd
	Character(len=255) :: bufferLongA
	Character(len=255) :: fileRowData
	Character(*) :: inputEamFile
!Deallocate already allocated arrays
    If(Allocated(elements))Then
	  Deallocate(elements)
	End If
    If(Allocated(elementsTemp))Then
	  Deallocate(elementsTemp)
	End If
    If(Allocated(eamKey))Then
	  Deallocate(eamKey)
	End If
    If(Allocated(eamData))Then
	  Deallocate(eamData)
	End If	
!Set EAM type
    eamType = 1		!1 = default type/standard   2 = 2BMEAM	
!allocate elements array
    Allocate(elements(1:300))
    Allocate(elementsTemp(1:300))
!"blank" elements array
	do i=1,size(elements)
	  elements(i) = "ZZ"
	enddo	
!Count file rows
	fileRows = 0
  	Open(UNIT=1,FILE=trim(inputEamFile)) 
    do i=1,maxFileRows 
!count file rows
	  fileRows = fileRows + 1
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
	enddo
!close file	
	CLOSE(1)  		
!Read eam pot file and store unique elements
	potentialCounter = 0
  	Open(UNIT=1,FILE=trim(inputEamFile)) 
    do i=1,fileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
!skip blank or commented line      
	  if(buffera(1:1).eq." ".or.buffera(1:1).eq."#".or.buffera(1:1).eq."!")then
!skip
	  else 
!check if PAIR DENS or EMBE
        if(StrToUpper(buffera(1:4)).eq."PAIR")then
	      BACKSPACE(1)
		  Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		  potentialCounter = potentialCounter + 1
		  functionResult = AddUniqueElement(bufferb(1:2))
		  functionResult = AddUniqueElement(bufferc(1:2))
	    elseif(&
	    StrToUpper(buffera(1:4)).eq."DENS".or.StrToUpper(buffera(1:4)).eq."EMBE"&
	    .or.StrToUpper(buffera(1:4)).eq."EMBS"&
	    .or.StrToUpper(buffera(1:4)).eq."DEND".or.StrToUpper(buffera(1:4)).eq."EMBD"&
	    )then	  
!adjust potential type 
	      if(&
	        StrToUpper(buffera(1:4)).eq."DEND".or.StrToUpper(buffera(1:4)).eq."EMBD"&
		    .or.StrToUpper(buffera(1:4)).eq."EMBS"&
	      )then	
		    eamType = 2	!2BMEAM
		  endif
!re-read line, count potential and add element if unique
	      BACKSPACE(1)
		  Read(1,*,IOSTAT=ios) buffera, bufferb
		  potentialCounter = potentialCounter + 1
		  functionResult = AddUniqueElement(bufferb(1:2))
	    endif
	  endif
	enddo
!close file	
	CLOSE(1)   	  
!adjust elements array
    elementCounter = 0
	do i=1,size(elements)
	  if(elements(i).eq."ZZ")then
	    exit
	  else 	  
	    elementCounter = elementCounter + 1
		elementsTemp(i) = elements(i)
	  endif
	enddo
	Deallocate(elements)
	Allocate(elements(1:elementCounter))
	do i=1,elementCounter
	  elements(i) = elementsTemp(i)
	enddo
	Deallocate(elementsTemp)
	elementCount = elementCounter	
!store/output elements	
    If(mpiProcessID.eq.0)Then
      write(999,"(A16,F8.4)") "Program time:   ",ProgramTime()
	  Do i=1,size(elements)
	    write(999,"(A8,I4,A2,A2)") "Element ",i,": ",elements(i)
	  End Do
	End If			
!count data points
	dataCounter = 0	
Open(UNIT=1,FILE=trim(inputEamFile)) 
    do i=1,fileRows 
!Read in line
	  Read(1,'(A64)',IOSTAT=ios) fileRowData
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  endif
!skip blank or commented line   
	  if(fileRowData(1:8).ne."        ".and.fileRowData(1:1).ne."#".and.fileRowData(1:1).ne."!")then
!check if PAIR DENS or EMBE - store pot data start/length	
		read(fileRowData,*) potTypeText 
		potTypeText = trim(potTypeText)
	    if(StrToUpper(potTypeText(1:4)).eq."PAIR".or.&
	    StrToUpper(potTypeText(1:4)).eq."DENS".or.&
	    StrToUpper(potTypeText(1:4)).eq."SDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."DDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."EMBE".or.&
	    StrToUpper(potTypeText(1:4)).eq."SEMB".or.&
	    StrToUpper(potTypeText(1:4)).eq."DEMB")then
!potential header line          
		else		
	      dataCounter = dataCounter + 1
		endif
      endif
	enddo
	CLOSE(1) 
!potential counts
    if(eamType.eq.1)then
	  numberPotentials = (elementCounter * ( elementCounter + 5)) / 2
	  pairCount = (elementCounter * ( elementCounter + 1)) / 2
	  densCount = elementCounter
	  embeCount = elementCounter
	elseif(eamType.eq.2)then
	  if(elementCounter.eq.1)then
	    numberPotentials = 1 + (elementCounter * ( elementCounter + 7)) / 2
		pairCount = (elementCounter * ( elementCounter + 1)) / 2
	    densCount = 1
		dendCount = elementCounter
	    embeCount = elementCounter
	    embdCount = elementCounter
	  else
	    numberPotentials = (elementCounter * ( elementCounter + 3))
		pairCount = (elementCounter * ( elementCounter + 1)) / 2
	    densCount = (elementCounter * ( elementCounter - 1)) / 2
		dendCount = elementCounter
	    embeCount = elementCounter
	    embdCount = elementCounter
	  endif
	endif	
!store/output potential count etc
!Allocate eam keys/data
    Allocate(eamKey(1:numberPotentials,1:5))	
    Allocate(eamData(1:dataCounter,1:3))
!Initialise data array
    Do i=1,size(eamData,1)
      Do j=1,size(eamData,2)
	    eamData(i,j) = 0.0D0
	  End Do
	End Do
	potentialCounter = 0
	eamType = 1
	dataCounter = 0
	potDataStart = 1
	potDataLength = 0
	Open(UNIT=1,FILE=trim(inputEamFile)) 
    do i=1,fileRows 
!Read in line
	  Read(1,'(A64)',IOSTAT=ios) fileRowData
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  endif
!skip blank or commented line   
	  If(fileRowData(1:8).ne."        ".and.fileRowData(1:1)&
	    .ne."#".and.fileRowData(1:1).ne."!")Then
!check if PAIR DENS or EMBE - store pot data start/length	
		read(fileRowData,*) potTypeText 
		potTypeText = trim(potTypeText)
	    if(StrToUpper(potTypeText(1:4)).eq."PAIR".or.&
	    StrToUpper(potTypeText(1:4)).eq."DENS".or.&
	    StrToUpper(potTypeText(1:4)).eq."SDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."DDEN".or.&
	    StrToUpper(potTypeText(1:4)).eq."EMBE".or.&
	    StrToUpper(potTypeText(1:4)).eq."SEMB".or.&
	    StrToUpper(potTypeText(1:4)).eq."DEMB")then
!read potential header line
          potentialCounter = potentialCounter + 1
	      if(StrToUpper(potTypeText(1:4)).eq."PAIR".or.&
	      StrToUpper(potTypeText(1:4)).eq."SDEN")then
		    Read(fileRowData,*) buffera, bufferb, bufferc
		  else
		    Read(fileRowData,*) buffera, bufferb
		  endif		  
!save start/length details of last potential         
          if(potentialCounter.gt.1)then
	        eamKey(potKey,4) = potDataStart
	        eamKey(potKey,5) = potDataLength
		  endif
!get potential key
          potKey = 0
		  potType = 0 
		  if(StrToUpper(potTypeText(1:4)).eq."PAIR")then	
		    potType = 1	    
		    atomA = QueryUniqueElement(bufferb(1:2))
		    atomB = QueryUniqueElement(bufferc(1:2))
            atomMax = max(atomA,atomB)-1
            atomMin = min(atomA,atomB)-1
		    potKey = 1+atomMin+(atomMax*(atomMax+1))/2
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."DENS")then	
			potType = 2	    
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."SDEN")then	
			potType = 2
            eamType = 2		  
		    atomA = QueryUniqueElement(bufferb(1:2))
		    atomB = QueryUniqueElement(bufferc(1:2))
			if(elementCounter.eq.1)then
              potKey = pairCount + 1
            else
              atomMax = max(atomA,atomB)-1
              atomMin = min(atomA,atomB)-1
              potKey = pairCount+1+atomMin+(atomMax*(atomMax-1))/2
            endif
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."DDEN")then	
			potType = 4
            eamType = 2		  
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
			potKey = pairCount + densCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."EMBE")then	
			potType = 3	    
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + densCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."SEMB")then	
			potType = 3	    
            eamType = 2
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + densCount + dendCount + atomA
		  endif
		  if(StrToUpper(potTypeText(1:4)).eq."DEMB")then	
			potType = 5	    
            eamType = 2
		    atomA = QueryUniqueElement(bufferb(1:2))
			atomB = 0
		    potKey = pairCount + densCount + dendCount + embdCount + atomA
		  endif
!store potential details
          eamKey(potKey,1) = atomA
          eamKey(potKey,2) = atomB
          eamKey(potKey,3) = potType
!reset start/length details
		  potDataStart = potDataStart + potDataLength
		  potDataLength = 0
		else		
	      dataCounter = dataCounter + 1
		  potDataLength = potDataLength + 1
		  read(fileRowData,*) buffera, bufferb
		  read(buffera,*) eamData(dataCounter,1)
		  read(bufferb,*) eamData(dataCounter,2)		  	  
		endif
      endif
	enddo
!save last start/length entry
	eamKey(potKey,4) = potDataStart
	eamKey(potKey,5) = potDataLength
	CLOSE(1) 	
!-----------------------------
! Deallocate arrays
!-----------------------------
	If(Allocated(elementsTemp))Then
	  Deallocate(elementsTemp)
	End If		
  End Subroutine reReadEamPot
!================================================================================================================================================  
!================================================================================================================================================
!---------------------------------------------------------------------------------------------------
!
!    Read Configuration File
!
!---------------------------------------------------------------------------------------------------  
  Subroutine copyConfigFile()
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, fileRows
    Character(len=255)  :: configurationsFilePathTemp
    Character(len=255)  :: fileRow
!temp config file
    configurationsFilePathTemp = trim(configurationsFilePath)//".temp"
!make temp config file
    Open(unit=301,file=trim(configurationsFilePathTemp))
	Open(UNIT=1,FILE=configurationsFilePath)   	
    Do i=1,maxFileRows 
	  Read(1,"(A255)",IOSTAT=ios) fileRow
	  If (ios /= 0) Then !break out if end of file
	    EXIT 
	  End If
	  write(301,"(A255)") fileRow
	End Do 
!close input file
	close(1)
!close output file
	close(301)		
  End Subroutine copyConfigFile
!================================================================================================================================================  
!================================================================================================================================================	
!Read in configuration file
  Subroutine readDFTFiles()
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, fileRows
	Character(len=255) :: filePath
!Read PWSCF files
    If(size(pwscfFilesList,1).gt.0)Then
      Do i=1,size(pwscfFilesList,1)
	    filePath = pwscfFilesList(i)
	    If(filePath(1:9).eq."__empty__")Then
	      Exit
		Else
	      filePath = TrimSpaces(pwscfFilesList(i))
	      Call readPWSCFFile(filePath)			
		End If
	  End Do
	End If		
  End Subroutine readDFTFiles
!================================================================================================================================================  
!================================================================================================================================================
!Read in configuration file
  Subroutine readPWSCFFile(filePath)
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, fileRows
	Character(len=255) :: filePath
	Character(len=255) :: fileLineBuffer
	Character(len=255) :: tempBuffer, bufferA, bufferB, bufferC
	Character(len=255) :: bufferD, bufferE, bufferF
	Character(len=3) :: atomTypeBuffer
    Integer(kind=StandardInteger) :: readType, lastScf
	Integer(kind=StandardInteger) :: numberOfAtoms
	Real(kind=DoubleReal) :: aLat
	Real(kind=DoubleReal) :: configTotalEnergy
	Real(kind=DoubleReal) :: configEnergyPerAtom
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: stress
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: crystalAxes	!Unit vector
	Character(len=3), Dimension( : ), Allocatable :: atomType
	Real(kind=SingleReal), Dimension( : , : ), Allocatable :: atomCoords
	Real(kind=SingleReal), Dimension( : , : ), Allocatable :: atomForcess
	Character(len=255) :: configurationsFilePathTemp
!set temp config file name
    configurationsFilePathTemp = trim(configurationsFilePath)//".temp"
!set variables
    readType = 0    	
!Check file type
  	Open(UNIT=1,FILE=filePath) 
    Do i=1,maxFileRows 
!Read in line
	  Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
	  If(fileLineBuffer(1:54).eq.&
	  "     A final scf calculation at the relaxed structure.")Then
	    readType = 1
	  End If
	  If(fileLineBuffer(1:33).eq.&
	  "     number of atoms/cell      = ")Then
	    tempBuffer = fileLineBuffer(34:100)
		read(tempBuffer,*) numberOfAtoms
	  End If	  
!break out if end of file
	  If (ios /= 0) then
	    EXIT 
	  End If
	End Do
!close file	
	CLOSE(1) 	
!Allocate arrays
    Allocate(atomType(1:numberOfAtoms))
    Allocate(atomCoords(1:numberOfAtoms,1:3))
    Allocate(atomForcess(1:numberOfAtoms,1:3))	
!Read in data
	If(readType.eq.1)Then
	  lastScf = 0
	Else
	  lastScf = 1
	End If
  	Open(UNIT=1,FILE=filePath) 
    Do i=1,maxFileRows 
!Read in line
	  Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
	  If(fileLineBuffer(1:23).eq.&
	  "Begin final coordinates")Then
	    lastScf = 1
	  End If
      If(lastScf.eq.1)Then
!Total energy
		If(fileLineBuffer(1:17).eq."!    total energy")Then
		  tempBuffer = fileLineBuffer(33:100)
		  read(tempBuffer,*) bufferA, bufferB
		  read(bufferA,*) configTotalEnergy
		  configTotalEnergy = UnitConvert(configTotalEnergy,"RY","EV")
		  configEnergyPerAtom = (configTotalEnergy/(1.0D0*numberOfAtoms))-&
		    (dftInOptEnergy-dftInCohEnergy)  
		End If
!Stresses		
		If(fileLineBuffer(1:38).eq."          total   stress  (Ry/bohr**3)")Then
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  read(fileLineBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) stress(1,1)
		  read(bufferB,*) stress(1,2)
		  read(bufferC,*) stress(1,3)
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  read(fileLineBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) stress(2,1)
		  read(bufferB,*) stress(2,2)
		  read(bufferC,*) stress(2,3)
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  read(fileLineBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) stress(3,1)
		  read(bufferB,*) stress(3,2)
		  read(bufferC,*) stress(3,3)
		  Do j=1,3
		    Do k=1,3
			  stress(j,k) = UnitConvert(stress(j,k),"RYBOH3","GPA")
			End Do
          End Do			
		End If
!Lattice Parameter		
	    If(fileLineBuffer(1:32).eq."     lattice parameter (alat)  =")Then
		  tempBuffer = fileLineBuffer(33:100)
		  read(tempBuffer,*) bufferA, bufferB
		  read(bufferA,*) aLat
		  aLat = UnitConvert(aLat,"Bohr","A")
		End If
!Crystal Axes/Unit Vector - scf
        If(readType.eq.0.and.fileLineBuffer(1:18).eq."     crystal axes:")Then
		  Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  tempBuffer = fileLineBuffer(24:100)
		  read(tempBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) crystalAxes(1,1)
		  read(bufferB,*) crystalAxes(1,2)
		  read(bufferC,*) crystalAxes(1,3)
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  tempBuffer = fileLineBuffer(24:100)
		  read(tempBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) crystalAxes(2,1)
		  read(bufferB,*) crystalAxes(2,2)
		  read(bufferC,*) crystalAxes(2,3)
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  tempBuffer = fileLineBuffer(24:100)
		  read(tempBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) crystalAxes(3,1)
		  read(bufferB,*) crystalAxes(3,2)
		  read(bufferC,*) crystalAxes(3,3)
		End If
!Crystal Axes/Unit Vector - vc-relax
        If(readType.eq.1.and.fileLineBuffer(1:15).eq."CELL_PARAMETERS")Then
		  Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  tempBuffer = fileLineBuffer(1:100)
		  read(tempBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) crystalAxes(1,1)
		  read(bufferB,*) crystalAxes(1,2)
		  read(bufferC,*) crystalAxes(1,3)
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  tempBuffer = fileLineBuffer(1:100)
		  read(tempBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) crystalAxes(2,1)
		  read(bufferB,*) crystalAxes(2,2)
		  read(bufferC,*) crystalAxes(2,3)
	      Read(1,"(A255)",IOSTAT=ios) fileLineBuffer
		  tempBuffer = fileLineBuffer(1:100)
		  read(tempBuffer,*) bufferA, bufferB, bufferC
		  read(bufferA,*) crystalAxes(3,1)
		  read(bufferB,*) crystalAxes(3,2)
		  read(bufferC,*) crystalAxes(3,3)
		End If		
!Symbol/Coords - scf
        If(readType.eq.0.and.fileLineBuffer(1:22).eq."     site n.     atom ")Then
		  Do j=1,numberOfAtoms
		    Read(1,"(A255)",IOSTAT=ios) fileLineBuffer	
			atomType(j) = TrimSpaces(fileLineBuffer(11:23))
			Do k=1,size(dftInReplaceSym,1)
			  If(Adjustl(StrToUpper(dftInReplaceSym(k,1))).eq.&
			    Adjustl(StrToUpper(atomType(j))))Then
				atomType(j) = dftInReplaceSym(k,2)
				Exit
			  End If
		    End Do
			tempBuffer = fileLineBuffer(39:75)
		    read(tempBuffer,*) bufferA, bufferB, bufferC
		    read(bufferA,*) atomCoords(j,1)
		    read(bufferB,*) atomCoords(j,2)
		    read(bufferC,*) atomCoords(j,3)
		  End Do
		End If
!Symbol/Coords - vcrelax
		If(fileLineBuffer(1:16).eq."ATOMIC_POSITIONS")Then
		  Do j=1,numberOfAtoms
		    Read(1,"(A255)",IOSTAT=ios) fileLineBuffer	
			atomType(j) = TrimSpaces(fileLineBuffer(1:7))
			Do k=1,size(dftInReplaceSym,1)
			  If(Adjustl(StrToUpper(dftInReplaceSym(k,1))).eq.&
			    Adjustl(StrToUpper(atomType(j))))Then
				atomType(j) = dftInReplaceSym(k,2)
				Exit
			  End If
		    End Do
			tempBuffer = fileLineBuffer(8:48)
		    read(tempBuffer,*) bufferA, bufferB, bufferC
		    read(bufferA,*) atomCoords(j,1)
		    read(bufferB,*) atomCoords(j,2)
		    read(bufferC,*) atomCoords(j,3)
		  End Do
		End If		
!atom forces
        If(fileLineBuffer(1:27).eq."     Forces acting on atoms")Then
		  Read(1,"(A255)",IOSTAT=ios) fileLineBuffer !read blank line
		  Do j=1,numberOfAtoms
		    Read(1,"(A255)",IOSTAT=ios) fileLineBuffer	
			tempBuffer = fileLineBuffer(33:76)
		    read(tempBuffer,*) bufferA, bufferB, bufferC
		    read(bufferA,*) atomForcess(j,1)
		    read(bufferB,*) atomForcess(j,2)
		    read(bufferC,*) atomForcess(j,3)
		  End Do
		  Do j=1,numberOfAtoms
		    Do k=1,3
			  atomForcess(j,k) = 1.0E0*UnitConvert(1.0D0*atomForcess(j,k),"RYBOHR","EVANG")
			End Do  
		  End Do
		End If
	  End If
!break out if end of file
	  If (ios /= 0) then
	    EXIT 
	  End If
	End Do
!close file	
	CLOSE(1) 
!Write to file	
    If(mpiProcessID.eq.0)Then	
!open temp config file	
	  open(unit=801,file=trim(configurationsFilePathTemp),&
	  status="old",position="append",action="write")	
!write to output file
	  write(801,"(A4)") "#new"	
	  write(801,"(A21,I8,A9,E20.10)") "#added config, atoms ",numberOfAtoms,&
	    ", energy ",configTotalEnergy 	
	  write(801,"(A4,E20.10)") "#LP ",aLat
	  write(801,"(A3,F12.7,F12.7,F12.7)") "#X ",crystalAxes(1,1),crystalAxes(1,2),crystalAxes(1,3)
	  write(801,"(A3,F12.7,F12.7,F12.7)") "#Y ",crystalAxes(2,1),crystalAxes(2,2),crystalAxes(2,3)
	  write(801,"(A3,F12.7,F12.7,F12.7)") "#Z ",crystalAxes(3,1),crystalAxes(3,2),crystalAxes(3,3)
	  write(801,"(A4,F12.7,F12.7,F12.7)") "#SX ",stress(1,1),stress(1,2),stress(1,3)
	  write(801,"(A4,F12.7,F12.7,F12.7)") "#SY ",stress(2,1),stress(2,2),stress(2,3)
	  write(801,"(A4,F12.7,F12.7,F12.7)") "#SZ ",stress(3,1),stress(3,2),stress(3,3)
	  write(801,"(A9)") "#CC 1 1 1"
	  write(801,"(A4,F10.7)") "#RC ",dftInRadiusCutoff
	  write(801,"(A5,F10.7,A3)") "#EPA ",configEnergyPerAtom," eV"
	  write(801,"(A4)") "#F Y"
	  Do i=1,size(atomType,1)
	    write(801,"(A3,A2,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10,A1,F16.10)") &
		atomType(i),"  ",atomCoords(i,1)," ",atomCoords(i,2)," ",atomCoords(i,3),&
		" ",atomForcess(i,1)," ",atomForcess(i,2)," ",atomForcess(i,3)
	  End Do
	  write(801,"(A4)") "#end"	
	  Close(801)
    End If
!-----------------------------
! Deallocate arrays
!-----------------------------
	If(Allocated(atomCoords))Then
	  Deallocate(atomCoords)
	End If	
	If(Allocated(atomForcess))Then
	  Deallocate(atomForcess)
	End If			
  End Subroutine readPWSCFFile
!================================================================================================================================================  
!================================================================================================================================================
!Read in configuration file
  Subroutine readConfiguration()	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, fileRows
	Integer(kind=StandardInteger) :: configurationCount, atomCount, ecCount
	Integer(kind=StandardInteger) :: startConfig, configLength, readForces
	Integer(kind=StandardInteger) :: absolutePosition
	Real(kind=DoubleReal) :: bufferDa
	Character(len=128) :: buffera, bufferb, bufferc, bufferd, buffere
	Character(len=128) :: bufferf, bufferg, bufferh, bufferi, bufferj
	Character(len=255) :: bufferLongA
	Character(len=255) :: fileRowBuffer    
	Character(len=255) :: configurationsFilePathTemp
!Use the temp file as the configuration input file
    configurationsFilePath = trim(configurationsFilePath)//".temp"	
!open output file	
    If(mpiProcessID.eq.0)Then
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
!write to output file
	  write(999,"(F8.4,A2,A18,A60)") ProgramTime(),"  ",&
	  "Reading from file ",configurationsFilePath
	End If
!Count file rows
	fileRows = 0
  	Open(UNIT=1,FILE=configurationsFilePath) 
    do i=1,maxFileRows 
!count file rows
	  fileRows = fileRows + 1
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
	enddo
!close file	
	CLOSE(1) 	
!count data amounts + get array sizes
	configurationCount = 0
	atomCount = 0
  	Open(UNIT=1,FILE=configurationsFilePath) 
    Do i=1,fileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera	  
!skip blank or commented lines
	  If(buffera(1:1).eq."!".or.buffera(1:10).eq."          ")Then
	    !skip
	  Else
!Check for the header of a configuration
        If(StrToUpper(buffera(1:4)).eq."#NEW")then
!count new configuration
		  configurationCount = configurationCount + 1		
	    End If
        If(buffera(1:1).ne."#".and.CheckIfElement(buffera(1:2)).eq.1)then
	      atomCount = atomCount + 1
	    End If
	  End If
	End Do
!close file	
	CLOSE(1) 
!allocate arrays
    configCount = configurationCount
	Allocate(configHeaderI(1:configurationCount,1:headerWidth))
	Allocate(configHeaderR(1:configurationCount,1:50)) !36 used
	Allocate(configCoordsI(1:atomCount))
	Allocate(configCoordsR(1:atomCount,1:3))
	Allocate(configForcesR(1:atomCount,1:4))
!Initialise arrays
	Do i=1,atomCount
	  configCoordsI(i) = 0
	  configCoordsR(i,1) = 0.0D0
	  configCoordsR(i,2) = 0.0D0
	  configCoordsR(i,3) = 0.0D0
	  configForcesR(i,1) = 0.0D0
	  configForcesR(i,2) = 0.0D0
	  configForcesR(i,3) = 0.0D0
	  configForcesR(i,4) = -1.0D0
	End Do	
!Defaults	
	readForces = 0
	absolutePosition = 0
!Default
    Do i=1,configurationCount
!References
	  configHeaderR(i,12) = -2.1E20		!Set no ref energy
	  configHeaderR(i,13) = -2.1E20		!Set no ref forces
	  configHeaderR(i,14) = -2.1E20		!
	  configHeaderR(i,15) = -2.1E20		!Set no bulk modulus
!Weighting
      configHeaderR(i,16) = 1.0E0		!configuration energy weighting
      configHeaderR(i,17) = 1.0E0		!configuration stresses weighting
      configHeaderR(i,18) = 1.0E0		!configuration forces weighting
      configHeaderR(i,19) = 1.0E0		!configuration stresses weighting
!Reference elastic constants
	  configHeaderR(i,20) = -2.1E20		!No. of elastic constants
	  configHeaderR(i,21) = -2.1E20		!C11   C11   C11
	  configHeaderR(i,22) = -2.1E20		!C12   C12   C12
	  configHeaderR(i,23) = -2.1E20		!C44   C13   C13
	  configHeaderR(i,24) = -2.1E20		!      C33   C33
	  configHeaderR(i,25) = -2.1E20		!      C55   C44
	  configHeaderR(i,26) = -2.1E20		!            C66
!References - Equilibrium volume	  
	  configHeaderR(i,27) = -2.1E20
!References - Equilibrium volume	  
	  configHeaderR(i,28) = -2.1D20		!SXX
	  configHeaderR(i,29) = -2.1D20		!SXY
	  configHeaderR(i,30) = -2.1D20		!SXZ
	  configHeaderR(i,31) = -2.1D20		!SYZ
	  configHeaderR(i,32) = -2.1D20		!SYY
	  configHeaderR(i,33) = -2.1D20		!SYZ
	  configHeaderR(i,34) = -2.1D20		!SZX
	  configHeaderR(i,35) = -2.1D20		!SZY
	  configHeaderR(i,36) = -2.1D20		!SZZ
	End Do
!Load Data
    configurationCount = 0
	atomCount = 0
	startConfig = 1
	configLength = 0
  	Open(UNIT=1,FILE=configurationsFilePath) 
    Do i=1,fileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!skip blank or commented lines
	  If(buffera(1:1).eq."!".or.buffera(1:10).eq."          ")Then
	    !skip
	  Else
!Check for the header of a configuration
        if(StrToUpper(buffera(1:4)).eq."#NEW")then
!Set defaults for new configuration
	      readForces = 0
!count new configuration, adjust startConfig, reset configLength
		  configurationCount = configurationCount + 1	
		  startConfig = startConfig + configLength
		  configLength = 0	
	    End If
	    if(StrToUpper(buffera(1:4)).eq."#END")then
!store header information start/length
          configHeaderI(configurationCount,headerWidth-1) = startConfig
          configHeaderI(configurationCount,headerWidth) = configLength	
	    End If
	    if(StrToUpper(buffera(1:3)).eq."#LP")then     !Lattice parameter
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb
		  Read(bufferb,*) configHeaderR(configurationCount,1)
	    End If
	    if(StrToUpper(buffera(1:4)).eq."#EPA")then     !Energy per atom
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		  Read(bufferb,*) bufferDa 
		  configHeaderR(configurationCount,12) = &
		    UnitConvert(bufferDa,TrimSpaces(bufferc),"EV")
	    End If
	    if(StrToUpper(buffera(1:2)).eq."#F")then     !Read forces
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb
		  bufferb = StrToUpper(bufferb)
		  If(bufferb(1:1).eq."Y")Then
		    readForces = 1
			configHeaderR(configurationCount,13) = 1.0E0
		  End If
	    End If
		if(StrToUpper(buffera(1:3)).eq."#AP")then     !Absolute positions not crystal 
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb
		  bufferb = StrToUpper(bufferb)
		  If(bufferb(1:1).eq."Y")Then
		    absolutePosition = 1
		  End If
	    End If
!Unit vectors
	    If(StrToUpper(buffera(1:2)).eq."#X")then     !x unit vector
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		  Read(bufferb,*) configHeaderR(configurationCount,3)
		  Read(bufferc,*) configHeaderR(configurationCount,4)
		  Read(bufferd,*) configHeaderR(configurationCount,5)
	    End If
	    if(StrToUpper(buffera(1:2)).eq."#Y")then     !Y unit vector
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		  Read(bufferb,*) configHeaderR(configurationCount,6)
		  Read(bufferc,*) configHeaderR(configurationCount,7)
		  Read(bufferd,*) configHeaderR(configurationCount,8)
	    End If
	    if(StrToUpper(buffera(1:2)).eq."#Z")then     !Z unit vector
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		  Read(bufferb,*) configHeaderR(configurationCount,9)
		  Read(bufferc,*) configHeaderR(configurationCount,10)
		  Read(bufferd,*) configHeaderR(configurationCount,11)
	    End If
	    if(StrToUpper(buffera(1:3)).eq."#CC")then     !Cell Copies
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		  Read(bufferb,*) configHeaderI(configurationCount,10)
		  Read(bufferc,*) configHeaderI(configurationCount,11)
		  Read(bufferd,*) configHeaderI(configurationCount,12)
	    End If
!Read Stresses 
		If(StrToUpper(buffera(1:3)).eq."#SX")Then 
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd		  
		  Read(bufferb,*) configHeaderR(configurationCount,28)
		  Read(bufferc,*) configHeaderR(configurationCount,29)
		  Read(bufferd,*) configHeaderR(configurationCount,30)
		End If		
		If(StrToUpper(buffera(1:3)).eq."#SY")Then 
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd		  
		  Read(bufferb,*) configHeaderR(configurationCount,31)
		  Read(bufferc,*) configHeaderR(configurationCount,32)
		  Read(bufferd,*) configHeaderR(configurationCount,33)
		End If	
		If(StrToUpper(buffera(1:3)).eq."#SZ")Then 
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd		  
		  Read(bufferb,*) configHeaderR(configurationCount,34)
		  Read(bufferc,*) configHeaderR(configurationCount,35)
		  Read(bufferd,*) configHeaderR(configurationCount,36)
		End If
		
!Read equilibrium volume
	    If(StrToUpper(buffera(1:3)).eq."#EV")Then     !equilibrium volume		
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		  Read(bufferb,*) configHeaderR(configurationCount,27)
          configHeaderR(configurationCount,27) = &
		  UnitConvert(1.0D0*configHeaderR(configurationCount,27),bufferc,"ANG3")		  
		End If		
		
		
!Read bulk modulus
	    If(StrToUpper(buffera(1:3)).eq."#BM")Then     !bulk modulus		
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc
		  Read(bufferb,*) configHeaderR(configurationCount,15)
          configHeaderR(configurationCount,15) = &
		  UnitConvert(1.0D0*configHeaderR(configurationCount,15),bufferc,"GPA")		  
		End If
		
!Read elastic constants
	    If(StrToUpper(buffera(1:3)).eq."#EC")Then     !bulk modulus		
		  ecCount = 3  !Default
		  Do j=21,26
		    configHeaderR(configurationCount,j) = -2.1D20
		  End Do
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb
		  Read(bufferb,*) ecCount !count of elastic constants	
		  configHeaderR(configurationCount,20) = 1.0D0 * ecCount
!Cubic         C11 C12 C14
          If(ecCount.eq.3)Then
!re-read row
		    Backspace(1)		
	        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere, bufferf		  
		    Read(bufferc,*) configHeaderR(configurationCount,21)   !C11
            configHeaderR(configurationCount,21) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,21),bufferf,"GPA")			  
		    Read(bufferd,*) configHeaderR(configurationCount,22)   !C12
            configHeaderR(configurationCount,22) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,22),bufferf,"GPA")				  
		    Read(buffere,*) configHeaderR(configurationCount,23)   !C44
            configHeaderR(configurationCount,23) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,23),bufferf,"GPA")		
		  End If  
!Hexagonal     C11 C12 C13 C33 C55
          If(ecCount.eq.5)Then
!re-read row
		    Backspace(1)		
	        Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere, &
			bufferf, bufferg, bufferh		  
		    Read(bufferc,*) configHeaderR(configurationCount,21)   !C11
            configHeaderR(configurationCount,21) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,21),bufferh,"GPA")			  
		    Read(bufferd,*) configHeaderR(configurationCount,22)   !C12
            configHeaderR(configurationCount,22) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,22),bufferh,"GPA")				  
		    Read(buffere,*) configHeaderR(configurationCount,23)   !C13
            configHeaderR(configurationCount,23) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,23),bufferh,"GPA")					  
		    Read(bufferf,*) configHeaderR(configurationCount,24)   !C33
            configHeaderR(configurationCount,24) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,24),bufferh,"GPA")					  
		    Read(bufferg,*) configHeaderR(configurationCount,25)   !C55
            configHeaderR(configurationCount,25) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,25),bufferh,"GPA")	
		  End If  
!Tetragonal     C11 C12 C13 C33 C44 C66
          If(ecCount.eq.6)Then
!re-read row
		    Backspace(1)		
			Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere, &
			bufferf, bufferg, bufferh, bufferi		  		  
		    Read(bufferc,*) configHeaderR(configurationCount,21)   !C11
            configHeaderR(configurationCount,21) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,21),bufferi,"GPA")			  
		    Read(bufferd,*) configHeaderR(configurationCount,22)   !C12
            configHeaderR(configurationCount,22) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,22),bufferi,"GPA")				  
		    Read(buffere,*) configHeaderR(configurationCount,23)   !C13
            configHeaderR(configurationCount,23) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,23),bufferi,"GPA")					  
		    Read(bufferf,*) configHeaderR(configurationCount,24)   !C33
            configHeaderR(configurationCount,24) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,24),bufferi,"GPA")					  
		    Read(bufferg,*) configHeaderR(configurationCount,25)   !C44
            configHeaderR(configurationCount,25) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,25),bufferi,"GPA")					  
		    Read(bufferh,*) configHeaderR(configurationCount,26)   !C66
            configHeaderR(configurationCount,26) = &
		    UnitConvert(1.0D0*configHeaderR(configurationCount,26),bufferi,"GPA")	
		  End If  
		End If
		
		
if(StrToUpper(buffera(1:3)).eq."#CW")then     !Configuration weighting energy stress force bproperties
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd, buffere	
		  Read(bufferb,*) configHeaderR(configurationCount,16)
		  Read(bufferc,*) configHeaderR(configurationCount,17)
		  Read(bufferd,*) configHeaderR(configurationCount,18)
		  Read(buffere,*) configHeaderR(configurationCount,19)
	    End If
		
		
!Read potential radius cutoff
	    if(StrToUpper(buffera(1:3)).eq."#RC")Then     !radius cutoff
!re-read row
		  Backspace(1)		
	      Read(1,*,IOSTAT=ios) buffera, bufferb
		  Read(bufferb,*) configHeaderR(configurationCount,2)
	    End If
        if(buffera(1:1).ne."#".and.CheckIfElement(buffera(1:2)).eq.1)then
	      atomCount = atomCount + 1
		  configLength = configLength + 1
!re-read row
		  Backspace(1)
		  If(readForces.eq.0)Then	
		    Read(1,"(A255)",IOSTAT=ios) fileRowBuffer
		    Read(fileRowBuffer,*) buffera, bufferb, bufferc, bufferd
		    configCoordsI(atomCount) = QueryUniqueElement(buffera)
		    Read(bufferb,*) configCoordsR(atomCount,1)
		    Read(bufferc,*) configCoordsR(atomCount,2)
		    Read(bufferd,*) configCoordsR(atomCount,3)
		  End If
		  If(readForces.eq.1)Then	
		    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd,&
		    buffere, bufferf, bufferg
		    configCoordsI(atomCount) = QueryUniqueElement(buffera)
		    Read(bufferb,*) configCoordsR(atomCount,1)
		    Read(bufferc,*) configCoordsR(atomCount,2)
		    Read(bufferd,*) configCoordsR(atomCount,3)
		    Read(buffere,*) configForcesR(atomCount,1)
	        Read(bufferf,*) configForcesR(atomCount,2)
	        Read(bufferg,*) configForcesR(atomCount,3)
	        configForcesR(atomCount,4) = 1.0E0
		  End If
		  If(absolutePosition.eq.1)Then
!Adjust if absolute positions are given (divide by alat)
		    configCoordsR(atomCount,1) = configCoordsR(atomCount,1)/&
			  configHeaderR(configurationCount,1)
		    configCoordsR(atomCount,2) = configCoordsR(atomCount,2)/&
			  configHeaderR(configurationCount,1)
		    configCoordsR(atomCount,3) = configCoordsR(atomCount,3)/&
			  configHeaderR(configurationCount,1)
		  End If
		End If  
	  End If
	End Do
!close file	
	CLOSE(1) 	
	If(mpiProcessID.eq.0)Then
!write to output file
	  write(999,"(A6,A21,I8)") "      ","Configurations read: ",configurationCount
	  write(999,"(A6,A21,I8)") "      ","Atoms read:          ",atomCount
!close output file
      close(999)	
	End If  	
  End Subroutine readConfiguration
!================================================================================================================================================  
!================================================================================================================================================
!---------------------------------------------------------------------------------------------------
!
!    Other Module Subroutines
!
!--------------------------------------------------------------------------------------------------- 
  Subroutine synchMpiProcesses()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i,send,receive,tag 
    Integer(kind=StandardInteger) :: processID,processCount,error,status 
!call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!Send out data to all processes
    send = 0
    If(processID.eq.0) Then
      Do i=1,(processCount-1)
        tag = 1000 + i	
        Call MPI_send(send,1,MPI_integer,i,tag,&
		MPI_comm_world,error)
      End Do
    Else
      tag = 1000 + processID
	  Call MPI_recv(receive,1,MPI_integer,0,tag,&
      MPI_comm_world,status,error)  
    End If  
  End Subroutine synchMpiProcesses 
	
	
	
  
  
  
!================================================================================================================================================  
!================================================================================================================================================
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!    
  Function AddUniqueElement (element) RESULT (output)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: output 
	Integer(kind=StandardInteger) :: i, j, k, found
!convert to uppercase
    element = StrToUpper(element)  
!loop through elements array	
    k = 0    
	found = 0
	do i=1,size(elements)
	  k = k + 1
	  if(elements(i).eq."ZZ")then
	    exit
	  endif
	  if(elements(i).eq.element)then
	    found = 1
		exit
	  endif
	enddo
!save element if not found
    if(found.eq.0)then
      elements(k) = element
	  output = k
    endif
  End function AddUniqueElement    
!================================================================================================================================================  
!================================================================================================================================================  
  Function QueryUniqueElement (element) RESULT (output)
    Character(len=2) :: element
    Integer(kind=StandardInteger) :: output 
	Integer(kind=StandardInteger) :: i, j, k, found
!convert to uppercase
    element = StrToUpper(element)  
!loop through elements array	
    k = 0    
	found = 0
	do i=1,size(elements)
	  k = k + 1
	  if(elements(i).eq.element)then
	    found = 1
		exit
	  endif
	enddo
!save element if not found
    if(found.eq.1)then
	  output = k
	else
	  output = 0
    endif
  End function QueryUniqueElement  
!================================================================================================================================================  
!================================================================================================================================================  
  Function CheckIfElement (element) RESULT (output)
    Character(len=2) :: element
	Integer(kind=StandardInteger) :: output 
	Integer(kind=StandardInteger) :: i, j, k
	element = StrToUpper(element)
	output = 0
	if(iachar(element(1:1)).ge.65.and.iachar(element(1:1)).le.90&
	.and.iachar(element(2:2)).ge.65.and.iachar(element(2:2)).le.90)then
	  output = 1
	endif	
  End Function CheckIfElement
  

End Module input