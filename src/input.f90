Module input

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths
  Use initialise


!force declaration of all variables
  Implicit None
  
!declare global variables  
  Character(len=255)  :: inputFileName
  Character(len=255)  :: potentialFilePath
  Character(len=255)  :: configurationsFilePath
  Character(len=2), Dimension( : ), Allocatable :: elements
  Integer(kind=StandardInteger) :: elementCount
  Integer(kind=StandardInteger) :: eamType
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKey
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamData
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: unitVector
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: configHeaderI
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: configHeaderR
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configCoordsI
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: configCoordsR
  Integer(kind=StandardInteger), Parameter :: headerWidth = 14 
  Integer(kind=StandardInteger) :: numberPotentials
  Integer(kind=StandardInteger) :: pairCount, densCount, dendCount, embeCount, embdCount
  Character(len=1) :: saveFileCoords
  Character(len=1) :: saveFileNeighbourList
  Character(len=1) :: saveFilePot
  Integer(kind=StandardInteger) :: calcRunType
  Integer(kind=StandardInteger) :: configCount
  Integer(kind=StandardInteger) :: eamInterpType
  Integer(kind=StandardInteger) :: calcForcesOnOff
  
  
!Privacy of functions/subroutines/variables
  Private
  Public :: runInput				    !Subroutine  
  Public :: inputFileName			    !Variable
  Public :: potentialFilePath			!Variable
  Public :: configurationsFilePath		!Variable
  Public :: elements             		!Variable
  Public :: unitVector           		!Variable
  Public :: configHeaderI          		!Variable
  Public :: configHeaderR             	!Variable
  Public :: configCoordsI          		!Variable
  Public :: configCoordsR             	!Variable
  Public :: headerWidth             	!Variable
  Public :: numberPotentials           	!Variable
  Public :: eamType           	        !Variable
  Public :: pairCount           	    !Variable
  Public :: densCount           	    !Variable
  Public :: dendCount           	    !Variable
  Public :: embeCount           	    !Variable
  Public :: embdCount           	    !Variable
  Public :: elementCount           	    !Variable
  Public :: eamKey               	    !Variable
  Public :: eamData               	    !Variable
  Public :: saveFileCoords              !Variable
  Public :: saveFileNeighbourList       !Variable
  Public :: saveFilePot                 !Variable
  Public :: configCount                 !Variable
  Public :: calcRunType                 !Variable
  Public :: eamInterpType               !Variable
  Public :: calcForcesOnOff             !Variable

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
	Call readUserInput()
	Call readInputFile()
	Call readEamPot()
	Call readConfiguration()
		

  End Subroutine runInput
  
  
!read in user input data, input from the command line
  Subroutine readUserInput()
  !force declaration of all variables
	Implicit None
  
!Read in command line arguments
    call get_command_argument(1,inputFileName)
	
  
  End Subroutine readUserInput
  
  
  
  
!read in input file
  Subroutine readInputFile()      
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, elementCounter, headerRow
	Character(len=32) :: buffera, bufferb, bufferc, bufferd
	Character(len=255) :: bufferLongA
	
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
!write to output file
    If(mpiProcessID.eq.0)Then
	  write(999,"(F8.4,A2,A24,A60)") ProgramTime(),"  ",&
	  "Reading user input file ",inputFileName
	End If
	
!set defaults, overridden by user input
	Allocate(unitVector(1:3,1:3))
	unitVector(1,1) = 1			!x1
	unitVector(2,1) = 0			!x2
	unitVector(3,1) = 0			!x3
	unitVector(1,2) = 0			!y1
	unitVector(2,2) = 1			!y2
	unitVector(3,2) = 0			!y3
	unitVector(1,3) = 0			!z1
	unitVector(2,3) = 0			!z2
	unitVector(3,3) = 1			!z3
	eamInterpType = 1
	calcRunType = 1
	calcForcesOnOff = 1
	
!open & read in file	
  	Open(UNIT=1,FILE=inputFileName) 
    do i=1,maxFileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!break out if end of file
	  if (ios /= 0) then
	    EXIT 
	  end if
	  if(buffera(1:10).eq."#potential")then
!read next line
	    Read(1,*,IOSTAT=ios) potentialFilePath
	  endif
	  if(buffera(1:15).eq."#configurations")then
!read next line
	    Read(1,*,IOSTAT=ios) configurationsFilePath
	  endif
	  if(buffera(1:11).eq."#savecoords")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFileCoords = StrToUpper(buffera(1:1))
	  endif
	  if(buffera(1:18).eq."#saveneighbourlist")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFileNeighbourList = StrToUpper(buffera(1:1))
	  endif
	  if(buffera(1:8).eq."#savepot")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		saveFilePot = StrToUpper(buffera(1:1))
	  endif
	  
!Run types
	  if(buffera(1:8).eq."#runtype")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)
		if(buffera(1:3).eq."ENE")then		!ENERGY
		  calcRunType = 1
		endif
		if(buffera(1:6).eq."FOR")then		!FORCES
		  calcRunType = 2
		endif
		if(buffera(1:3).eq."BUL")then		!BULKMODULUS
		  calcRunType = 3
		endif
		if(buffera(1:3).eq."ELA")then		!ELASTICCONSTANTS
		  calcRunType = 4
		endif
		if(buffera(1:3).eq."OPT")then		!OPTIMISE
		  calcRunType = 5
		endif
	  endif
	  
!Unit vector components
	  if(buffera(1:11).eq."#unitvector")then
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
	  if(buffera(1:10).eq."#eaminterp")then
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
	  
	  
!Calc Forces
	  if(buffera(1:11).eq."#calcforces")then
!read next line
	    Read(1,*,IOSTAT=ios) buffera
		buffera = StrToUpper(buffera)	  
		If(buffera(1:1).eq."Y")Then
	      calcForcesOnOff = 1
		Else
		  calcForcesOnOff = 0
		End If
	  endif
	  
	  
    enddo
!close file	
	CLOSE(1) 
	  
!close output file
    close(999)	  
  
  End Subroutine readInputFile
  
  
  
  
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
	  "Reading EAM potential from ",potentialFilePath
	endif
	
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
  	Open(UNIT=1,FILE=potentialFilePath) 
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
  	Open(UNIT=1,FILE=potentialFilePath) 
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
Open(UNIT=1,FILE=potentialFilePath) 
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
	Open(UNIT=1,FILE=potentialFilePath) 
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
	
!output/store potentials	
    If(mpiProcessID.eq.0)Then
      write(999,"(A6,A32)") "      ","Summary of potential functions "
      do i=1,size(eamKey)/5
	    write(999,"(A6,I8,I8,I8,I8,I8,I8)") "      ",i,eamKey(i,1),eamKey(i,2),&
	    eamKey(i,3),eamKey(i,4),eamKey(i,5)
      enddo
	End If
	

	
!close output file
    close(999)	
	
  End Subroutine readEamPot
	
	
	
	
  
!Read in configuration file
  Subroutine readConfiguration()
	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, j, k, fileRows
	Integer(kind=StandardInteger) :: configurationCount, atomCount
	Integer(kind=StandardInteger) :: startConfig, configLength
	Character(len=32) :: buffera, bufferb, bufferc, bufferd
	Character(len=255) :: bufferLongA
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
    do i=1,fileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!Check for the header of a configuration
      if(StrToUpper(buffera(1:4)).eq."#NEW")then
!count new configuration
		configurationCount = configurationCount + 1		
	  endif
      if(buffera(1:1).ne."#".and.CheckIfElement(buffera(1:2)).eq.1)then
	    atomCount = atomCount + 1
	  endif
	enddo
!close file	
	CLOSE(1) 

!allocate arrays
    configCount = configurationCount
	Allocate(configHeaderI(1:configurationCount,1:headerWidth))
	Allocate(configHeaderR(1:configurationCount,1:2))
	Allocate(configCoordsI(1:atomCount))
	Allocate(configCoordsR(1:atomCount,1:3))
	
!Load Data
    configurationCount = 0
	atomCount = 0
	startConfig = 1
	configLength = 0
  	Open(UNIT=1,FILE=configurationsFilePath) 
    do i=1,fileRows 
!Read in line
	  Read(1,*,IOSTAT=ios) buffera
!Check for the header of a configuration
      if(StrToUpper(buffera(1:4)).eq."#NEW")then
!count new configuration, adjust startConfig, reset configLength
		configurationCount = configurationCount + 1	
		startConfig = startConfig + configLength
		configLength = 0	
	  endif
	  if(StrToUpper(buffera(1:4)).eq."#END")then
!store header information start/length
        configHeaderI(configurationCount,headerWidth-1) = startConfig
        configHeaderI(configurationCount,headerWidth) = configLength	
	  endif
	  if(StrToUpper(buffera(1:3)).eq."#LP")then     !Lattice parameter
!re-read row
		Backspace(1)		
	    Read(1,*,IOSTAT=ios) buffera, bufferb
		Read(bufferb,*) configHeaderR(configurationCount,1)
	  endif
	  if(StrToUpper(buffera(1:2)).eq."#X")then     !x unit vector
!re-read row
		Backspace(1)		
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		Read(bufferb,*) configHeaderI(configurationCount,1)
		Read(bufferc,*) configHeaderI(configurationCount,2)
		Read(bufferd,*) configHeaderI(configurationCount,3)
	  endif
	  if(StrToUpper(buffera(1:2)).eq."#Y")then     !Y unit vector
!re-read row
		Backspace(1)		
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		Read(bufferb,*) configHeaderI(configurationCount,4)
		Read(bufferc,*) configHeaderI(configurationCount,5)
		Read(bufferd,*) configHeaderI(configurationCount,6)
	  endif
	  if(StrToUpper(buffera(1:2)).eq."#Z")then     !Z unit vector
!re-read row
		Backspace(1)		
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		Read(bufferb,*) configHeaderI(configurationCount,7)
		Read(bufferc,*) configHeaderI(configurationCount,8)
		Read(bufferd,*) configHeaderI(configurationCount,9)
	  endif
	  if(StrToUpper(buffera(1:3)).eq."#CC")then     !Cell Copies
!re-read row
		Backspace(1)		
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd	
		Read(bufferb,*) configHeaderI(configurationCount,10)
		Read(bufferc,*) configHeaderI(configurationCount,11)
		Read(bufferd,*) configHeaderI(configurationCount,12)
	  endif
	  if(StrToUpper(buffera(1:3)).eq."#RC")then     !radius cutoff
!re-read row
		Backspace(1)		
	    Read(1,*,IOSTAT=ios) buffera, bufferb
		Read(bufferb,*) configHeaderR(configurationCount,2)
	  endif
      if(buffera(1:1).ne."#".and.CheckIfElement(buffera(1:2)).eq.1)then
	    atomCount = atomCount + 1
		configLength = configLength + 1
!re-read row
		Backspace(1)
	    Read(1,*,IOSTAT=ios) buffera, bufferb, bufferc, bufferd		
		configCoordsI(atomCount) = QueryUniqueElement(buffera)
		Read(bufferb,*) configCoordsR(atomCount,1)
		Read(bufferc,*) configCoordsR(atomCount,2)
		Read(bufferd,*) configCoordsR(atomCount,3)
	  endif
	enddo
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