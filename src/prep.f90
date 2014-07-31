Module prep

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use general
  Use strings		!string functions
  Use maths
  Use initialise
  Use input

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!declare global variables  
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: atomTypeKey                 !atom types by atom id
  Integer(kind=StandardInteger), Dimension(1:1024) :: configAtoms					!count of total atoms for each config
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtomsUnitCell			!
  Integer(kind=StandardInteger), Dimension(1:1024) :: configAtomsStart	
  Integer(kind=StandardInteger) :: configAtomsTotal											!Total atoms in all configs
  Integer(kind=StandardInteger) :: configAtomsMax											!Max atoms in 1 config
  Integer(kind=StandardInteger) :: configAtomsAdd											!Max rounded up to nearest mult of 10
  !Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: configAtomsMap			!Start-End points for atoms by config 
  Integer(kind=StandardInteger), Dimension(1:1024,1:50) :: configAtomsMap			!Start-End points for atoms by config  
  Integer(kind=StandardInteger), Dimension(1:10) :: globalCounter				!Array of global counts	
  Real(kind=SingleReal), Dimension(1:1024,1:3) :: configLatticeParameters
  Real(kind=DoubleReal), Dimension(1:50) :: globalTimer
  
!Neighbour list arrays
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKey		!key for atom neighbour list
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListI          !neighbour list integer data
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: neighbourListCoords
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListR
  
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKeyStore
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListIStore
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: neighbourListCoordsStore
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListRStore
  
!calculated values
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationVolume  
  !Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationEnergy
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationEnergy
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationBM
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configurationEC
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationEquVolume
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationEquLat
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptVolume
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptEnergy
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationForceX
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationForceY
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationForceZ
  Real(kind=DoubleReal), Dimension(1:9216) :: configurationStress
  Real(kind=DoubleReal), Dimension(1:9216) :: configurationUnitVector
  Real(kind=DoubleReal), Dimension(1:1000000,1:3) :: configurationForce
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configurationStresses
  
  
!Reference energies, forces, stresses etc  
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationRefEnergy
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationRefBM
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configurationRefEC
  Real(kind=DoubleReal), Dimension(1:1024) :: configurationRefEquVolume
  Real(kind=DoubleReal), Dimension(1:1000000,1:3) :: configurationRefForce
  Real(kind=DoubleReal), Dimension(1:1024,1:9) :: configurationRefStresses
  Real(kind=DoubleReal), Dimension(1:9216) :: configurationRefStress
  
!RSS configuration array    
  Real(kind=DoubleReal), Dimension(1:1024,1:10) :: configurationRSS
    
!RSS configuration opt weights    
  Real(kind=DoubleReal), Dimension(1:1024,1:10) :: configurationOptWeights
  
!Cutoff
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationRadiusCutoff  
  
!Reduced set of data points for optimization
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyReduced
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataReduced
  
!Temporary arrays for trial eam data points for optimization
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyTrial
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataTrial
    
!Optimum eam reduced data points so far
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOpt
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOpt
  
!Wobbliness
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: eamFunctionWobbliness
  
!Choice of potential data set
  Integer(kind=StandardInteger) :: eamDataSet
  
!Misc data variables
  Real(kind=DoubleReal) :: trialResidualSquareSum  
  Real(kind=DoubleReal) :: rssEnergyDifference, rssForceDifference
  Real(kind=DoubleReal) :: rssStressDifference, rssBPDifference
  
  
!MPI key for bringing forces back together  
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: configurationForceMPIKey
  

!Privacy of functions/subroutines/variables
  Private
  
!Variables - Define majority of global variables here
  Public :: runPrep	
  Public :: neighbourListKey,neighbourListI
  Public :: neighbourListR,neighbourListCoords	
  Public :: neighbourListKeyStore,neighbourListIStore
  Public :: neighbourListRStore,neighbourListCoordsStore		
  Public :: configAtoms	
  Public :: configAtomsUnitCell  
  Public :: configAtomsStart
  Public :: configAtomsMap	
  Public :: configLatticeParameters	
  Public :: configurationOptVolume	
  Public :: configurationOptEnergy	
  Public :: configAtomsTotal
  Public :: configAtomsMax
  Public :: configAtomsAdd
  Public :: configurationForceMPIKey
  Public :: configurationForce
  Public :: configurationStress   !remove in the future
  Public :: configurationStresses
!Configuration calculated properties 
  Public :: configurationVolume
  Public :: configurationEnergy	 
  Public :: configurationForceX
  Public :: configurationForceY
  Public :: configurationForceZ
  Public :: configurationBM	
  Public :: configurationEC	
  Public :: configurationEquVolume
  Public :: configurationEquLat
!Configuration reference properties
  Public :: configurationRefEnergy
  !Public :: configurationRefForceX
  !Public :: configurationRefForceY
  !Public :: configurationRefForceZ
  Public :: configurationRefBM
  Public :: configurationRefEC
  Public :: configurationRefEquVolume
  Public :: configurationRefForce
  Public :: configurationRefStress
  Public :: configurationRefStresses
!potential arrays
  Public :: atomTypeKey
  Public :: eamKeyReduced, eamDataReduced
  Public :: eamKeyTrial, eamDataTrial
  Public :: eamKeyOpt, eamDataOpt
  Public :: eamDataSet
  Public :: trialResidualSquareSum
  Public :: rssEnergyDifference, rssForceDifference
  Public :: rssStressDifference, rssBPDifference
  Public :: globalCounter, globalTimer
  Public :: configurationUnitVector
  Public :: configurationRadiusCutoff
  Public :: configurationRSS
  Public :: eamFunctionWobbliness
  Public :: configurationOptWeights  
!Subroutines  
  Public :: applyUnitVector
  Public :: applyUnitVectorAction
  Public :: makeTrialEAMSet
  Public :: makeReducedEAMSet
  Public :: makeReducedEAMSetTrial
  Public :: makeExpandedEAMSet
  Public :: storeReducedToOpt
  Public :: loadOptToReduced
  Public :: storeReducedToFile
  Public :: setPotentialDerivatives
  Public :: storeEAMToFile
  Public :: storeEAMToFileMaster
  Public :: printEAM
  Public :: storeToReduced
  Public :: synchMpiProcesses
  Public :: storeNeighboutList
  Public :: loadNeighboutList
  Public :: applyDistortionVector
  Public :: calcConfigurationVolume
  Public :: clearEAM,clearOptEAM,clearReducedEAM,clearTrialEAM
  Public :: splinePotential
  Public :: eamForceZBLCore
  Public :: storeReducedEAM
  Public :: loadReducedEAM
  
!Functions  
  Public :: eamWobbliness
  Public :: eamTotalWobbliness
  Public :: eamCurvature
  Public :: eamCurveLength
  

!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runPrep()
!force declaration of all variables
	Implicit None	
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
	Real(kind=DoubleReal) :: prepTime, prepTimeReal
!Start timer    
	Call MPI_timer(prepTime,prepTimeReal,0)
!--------------------------
! Prep data
!--------------------------
!Prep data 
    Call prepDataArrays()
	Call prepDataVariables()
!--------------------------
! Neighbour list
!--------------------------
!Make neighbour list
	Call makeNeighbourList()
!Calculate configuration volumes	
	Call synchMpiProcesses()
!--------------------------
! Potential
!--------------------------
!Sort out potential
	Call orderEamPotentials()
	If(mpiProcessID.eq.0)Then
	  Call storeEAM()
    End If
!--------------------------
! Synch and output
!--------------------------
	Call synchMpiProcesses()
	If(mpiProcessID.eq.0)Then
	  Call summariseConfigurations()
    End If	
	Call synchMpiProcesses()
    Call MPI_timer(prepTime,prepTimeReal,1)
	globalTimer(1) = prepTime
	globalTimer(2) = prepTimeReal
  End Subroutine runPrep  
!------------------------------------------------------------------------!
!                                                                        
! Prepare Data                                             
!                                                                        
!------------------------------------------------------------------------!     
  Subroutine prepDataArrays() 
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger) :: inputConfigAtomStart
	Integer(kind=StandardInteger) :: energyCounter, bulkModulusCounter, ecCounter
	Integer(kind=StandardInteger) :: evCounter
!-
!This sets up important data arrays needed here and elsewhere in the program
!-  
!Allocate arrays	  
	Allocate(configurationRadiusCutoff(1:configCount))		
! ConfigAtomsMap(i,1)		Start atom
! ConfigAtomsMap(i,2)		End atom
! ConfigAtomsMap(i,3)		Assigned MPI process
! ConfigAtomsMap(i,4)		Reference forces -1 no 1 yes
! ConfigAtomsMap(i,5)		Reference energy -1 no 1 yes
! ConfigAtomsMap(i,6)		Reference stress -1 no 1 yes
! ConfigAtomsMap(i,11)		Atom start in configCoordsR array
! ConfigAtomsMap(i,12)		Length/count in configCoordsR array
! ConfigAtomsMap(i,21)		Assigned MPI process
		
!zero out/fill arrays
    configAtomsTotal = 0		!total atoms in all configs
	configAtomsMax = 0			!max atoms in one config
	configAtomsStart(1) = 0		!start atom number for each config
	configurationRefEC = -2.1D20
	configurationEnergy = 0.0D0
	configurationRefEnergy = 0.0D0
	configurationRefForce = -2.1D20
	configurationRefStresses = -2.1D20
	configurationRSS = 0.0D0  !configurationRSS(i,j)
	configurationOptWeights = 0.0D0
!counters
    energyCounter = 0
	bulkModulusCounter = 0
	ecCounter = 0
	evCounter = 0
!Loop through configurations
    Do i=1,configCount
	  configurationVolume(i) = 0.0D0
      configurationEnergy(i) = 0.0D0	
	  configurationRefEnergy(i) = configHeaderR(i,12)
	  configurationBM(i) = -2.1D20
	  configurationRefBM(i) = configHeaderR(i,15)		!No reference data if less than -2.0D20
	  Do j=1,6
	    configurationRefEC(i,j) = configHeaderR(i,20+j)		!No reference data if less than -2.0D20
	  End Do	
!Store stresses
	  Do j=1,9
	    configurationRefStresses(i,j) = configHeaderR(i,27+j)
	  End Do
	  
	  
	  configurationRefEquVolume(i) = configHeaderR(i,27)
	  configurationEquVolume(i) = -2.1D20
	  configurationEquLat(i) = -2.1D20
	  Do j=1,6
	    configurationEC(i,j) = -2.1D20 
	  End Do
	  configAtoms(i) = configHeaderI(i,10) * configHeaderI(i,11) *&
	  configHeaderI(i,12) * configHeaderI(i,headerWidth)
	  configAtomsTotal = configAtomsTotal + configAtoms(i)
	  If(configAtoms(i).gt.configAtomsMax)Then
	    configAtomsMax = configAtoms(i)
	  End If
!config atoms start	
      If(i.lt.configCount)Then
	    configAtomsStart(i+1) = configAtomsStart(i) + configAtoms(i)
	  End If
!config atoms start-end map
	  If(i.eq.1)Then
        configAtomsMap(i,1) = 1					!start atom
        configAtomsMap(i,2) = configAtoms(i)    !end atom
	  Else
	    configAtomsMap(i,1) = configAtomsMap(i-1,2)+1
		configAtomsMap(i,2) = configAtomsMap(i,1)+configAtoms(i)-1		
	  End If
	  configAtomsMap(i,3) = -1    			!assigned process
	  configAtomsMap(i,11) =  configHeaderI(i,13)	!start atom in configCoordsR
	  configAtomsMap(i,12) =  configHeaderI(i,14)	!length/number of atoms in configCoordsR
!Set ref forces/energy/stress exist
	  If(configForcesR(configAtomsMap(i,11),4).lt.0)Then
        configAtomsMap(i,4) = -1    			!config has forces -1 no 1 yes
	  Else
	    configAtomsMap(i,4) = 1    			    !config has forces -1 no 1 yes
	  End If	  
	  If(configHeaderR(i,14).gt.-2.0D20)Then !If ref stress set
	    configAtomsMap(i,6) = 1 
	  Else
	    configAtomsMap(i,6) = -1 
	  End If  
	  If(configHeaderR(i,12).gt.-2.0D20)Then !If ref energy set
	    configAtomsMap(i,5) = 1 
	  Else
	    configAtomsMap(i,5) = -1 
	  End If	  
!MPI details 21-30
!Increment counters for MPI
	  energyCounter = energyCounter + 1
	  configAtomsMap(i,21) = mod(energyCounter-1,mpiProcessCount)		     !Assigned MPI process
	  If(configurationRefBM(i).gt.-2.0D20.or.calcBMCalcOnOff.eq.1)Then       !BM calc MPI process
	    bulkModulusCounter = bulkModulusCounter + 1
		configAtomsMap(i,22) = mod(bulkModulusCounter-1,mpiProcessCount)
	  Else
	    configAtomsMap(i,22) = -1
	  End If
	  If(configurationRefEC(i,1).gt.-2.0D20.or.calcECCalcOnOff.eq.1)Then       !EC calc MPI process
	    ecCounter = ecCounter + 1
		configAtomsMap(i,23) = mod(ecCounter-1,mpiProcessCount)
	  Else
	    configAtomsMap(i,23) = -1
	  End If
	  If(configurationRefEquVolume(i).gt.-2.0D20.or.calcEVCalcOnOff.eq.1)Then       !EV calc MPI process
	    evCounter = evCounter + 1
		configAtomsMap(i,24) = mod(evCounter-1,mpiProcessCount)
	  Else
	    configAtomsMap(i,24) = -1
	  End If	  
!Lattice parameters
	  configLatticeParameters(i,1) = configHeaderR(i,1) * configHeaderI(i,10)
	  configLatticeParameters(i,2) = configHeaderR(i,1) * configHeaderI(i,11)
	  configLatticeParameters(i,3) = configHeaderR(i,1) * configHeaderI(i,12)	
!Potential radius cutoff
      configurationRadiusCutoff(i) = configHeaderR(i,2)	  
!Configuration rates	  
      configurationOptWeights(i,1) = 1.0D0*configHeaderR(i,16)   !Configuration Energy Weight
      configurationOptWeights(i,2) = 1.0D0*configHeaderR(i,17)   !Configuration Force Weight
      configurationOptWeights(i,3) = 1.0D0*configHeaderR(i,18)   !Configuration Stress Weight  
      configurationOptWeights(i,4) = 1.0D0*configHeaderR(i,19)   !Configuration Bulk Properties Weight
!End loop through configurations
	End Do	
!If stress is to be calculated
	If(calcStressOnOff.eq.1)Then
	  Do i=1,(9*configCount)
	    configurationStress(i) = 0.0D0
		
	  End Do
	End If	
!Atom type key array
	Allocate(atomTypeKey(1:configAtomsTotal))
!Global counter - counting energy calculations, optimisation steps etc
!1 energy calc counter, 2 , 3 rss difference calc
	globalCounter = 0
!Global timer 	
	!1 prep time, 2 prep time real
	!3 energy calc time, 4 energy calc time real
	!5 bm calc, 6 bm calc real
	globalTimer = 0.0D0
!-----------------------------------------------
! open output file	
!-----------------------------------------------  
	open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
	status="old",position="append",action="write")
!write to output file
    If(mpiProcessID.eq.0)Then
	  write(999,"(F8.4,A2,A11)") ProgramTime(),"  ","Process Map"
	  write(999,"(A8,A5,A7,A8,A8,A8)") "        ","Cfg  ","E_Proc ","BM_Proc ","EC_Proc ",&
	  "EV_Proc "
	  Do i=1,configCount
	    write(999,"(A8,I4,A1,I4,A3,I4,A4,I4,A4,I4,A4)")&
		"        ",i," ",configAtomsMap(i,21),"   ",&
		configAtomsMap(i,22),"    ",configAtomsMap(i,23),"    ",configAtomsMap(i,24),"    "
	  End Do
	End If
!Close file
    close(999)  
  End Subroutine prepDataArrays
  
  
  
  
  Subroutine prepDataVariables() 
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k 
    eamDataSet = 1 !Default value, sets full EAM potential
	trialResidualSquareSum = -2.1D20	!Last calculation error measurement
  End Subroutine prepDataVariables

  
!------------------------------------------------------------------------!
!                                                                        
! Neighbour List Subroutines                                              
!                                                                        
!------------------------------------------------------------------------!  
!-------------------------- 
! Make Neighbour List
!--------------------------  
  Subroutine makeNeighbourList()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k, l, m, mm, n, x, y, z
	Integer(kind=StandardInteger) :: atoms
	Integer(kind=StandardInteger) :: atomA, atomB, nlKey
	Integer(kind=StandardInteger) :: neighbourListLength, tempNeighbourListLength
	Integer(kind=StandardInteger) :: neighbourListCount, configStart, configLength
	Integer(kind=StandardInteger) :: coordsStart, coordsLength
	Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
	Integer(kind=StandardInteger) :: atomCounter
	Real(kind=DoubleReal) :: rCutoff, rCutoffNeg, rCutoffSq
	Real(kind=DoubleReal) :: xcoord, ycoord, zcoord
	Real(kind=DoubleReal) :: xshift, yshift, zshift
	Real(kind=DoubleReal) :: xlat, ylat, zlat, alat
	Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB
	Real(kind=DoubleReal) :: rdSq, xdSq, ydSq, zdSq
	Integer(kind=StandardInteger), Dimension(1:100000) :: nlUniqueKeys
	Integer(kind=StandardInteger), Dimension(1:1024,1:2) :: neighbourListKeyTemp
	Integer(kind=StandardInteger), Dimension(1:1000000,1:6) :: neighbourListITemp
	Real(kind=DoubleReal), Dimension(1:1000000) :: neighbourListRTemp
    Real(kind=DoubleReal), Dimension(1:1000000,1:6) :: neighbourListCoordsTemp
	Integer(kind=StandardInteger), Dimension(1:100000) :: coordsITemp
	Real(kind=DoubleReal), Dimension(1:100000,1:3) :: coordsRTemp	
!Unit vectors
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: configUnitVector
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: workingUnitVector
!open files
	open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
	status="old",position="append",action="write")
	If(mpiProcessID.eq.0.and.saveFileCoords.eq."Y")Then
	  outputFile = trim(currentWorkingDirectory)//"/output/"//"coords.dat"
	  open(unit=10,file=trim(outputFile))	
	End If
	If(mpiProcessID.eq.0.and.saveFileNeighbourList.eq."Y")Then
	  outputFile = trim(currentWorkingDirectory)//"/output/"//"nl.dat"
	  open(unit=11,file=trim(outputFile))	
	End If
!write to output file
    If(mpiProcessID.eq.0)Then
      write(999,"(F8.4,A2,A37)") ProgramTime(),"  ",&
	  "Start building config neighbour lists"
	End If
!configHeaderR: 1 LP, 2 RC, configHeaderI: 1 X1, 2 X2, 3 X3, 4 Y1, 5 Y2, 6 Y3, 7 Z1, 8 Z2, 9 Z3, 10 CC1, 11 CC2, 12 CC3, 13 Start, 14 Length    
!Allocate temp neighbour list array - assume 100 neighbours per atom
	tempNeighbourListLength = 0
	Do i=1,configCount
      tempNeighbourListLength = tempNeighbourListLength + &
	  configHeaderI(i,10) * configHeaderI(i,11)* configHeaderI(i,12) * &
	  configHeaderI(i,headerWidth) * 100
	Enddo
	Allocate(configAtomsUnitCell(1:configCount)) 	
!loop through configurations - estimate size of
	configStart = 1
	neighbourListCount = 0
	atomCounter = 0
	neighbourListLength = 0
	Do i=1,configCount
!set variables      
      alat = configHeaderR(i,1)
	  xCopy = configHeaderI(i,10)
	  yCopy = configHeaderI(i,11)
	  zCopy = configHeaderI(i,12)	
!Set config unit vector
      configUnitVector(1,1) = configHeaderR(i,3)
      configUnitVector(1,2) = configHeaderR(i,4)
      configUnitVector(1,3) = configHeaderR(i,5)
      configUnitVector(2,1) = configHeaderR(i,6)
      configUnitVector(2,2) = configHeaderR(i,7)
      configUnitVector(2,3) = configHeaderR(i,8)
      configUnitVector(3,1) = configHeaderR(i,9)
      configUnitVector(3,2) = configHeaderR(i,10)
      configUnitVector(3,3) = configHeaderR(i,11)
!Make working vector
	  workingUnitVector = MatMul(globalUnitVector,configUnitVector)
!Store this as the configuration vector  
      n = 0
      Do j=1,3
	    Do k=1,3
		  n = n+1
	      configurationUnitVector((i-1)*9+n) = workingUnitVector(j,k)
		End Do
	  End Do
	  configurationVolume(i) = calcConfigurationVolume(i)
!expand unit cell	    
	  atoms = configAtoms(i)
	  configAtomsUnitCell(i) = configHeaderI(i,headerWidth)
!Blank temp arrays
      coordsITemp = 0
	  coordsRTemp = 0.0D0
	  m = 0
!loop over copies to make full crystal	  
	  coordsStart = configHeaderI(i,headerWidth-1)
	  coordsLength = configHeaderI(i,headerWidth)	 
	  If(mpiProcessID.eq.0.and.saveFileCoords.eq."Y")Then
	    write(10,"(A40)") "Configuration "//trim(intToString(i))
	  End If
	  do x=1,xCopy
	    do y=1,yCopy
	      do z=1,zCopy
		    do n=coordsStart,(coordsStart+coordsLength-1)  
			  m = m + 1
			  atomCounter = atomCounter + 1
!get co-ords with lattice parameter and copy applied
			  coordsITemp(m) = configCoordsI(n)
			  coordsRTemp(m,1) = alat*(x + configCoordsR(n,1) - 1)
			  coordsRTemp(m,2) = alat*(y + configCoordsR(n,2) - 1)
			  coordsRTemp(m,3) = alat*(z + configCoordsR(n,3) - 1)
!Apply configuration unit vector to these co-ordinates
!x-coord
              coordsRTemp(m,1) = & 		
			  coordsRTemp(m,1) * workingUnitVector(1,1) + &
              coordsRTemp(m,2) * workingUnitVector(1,2) + &
              coordsRTemp(m,3) * workingUnitVector(1,3)
!y-coord
              coordsRTemp(m,2) = & 		
			  coordsRTemp(m,1) * workingUnitVector(2,1) + &
              coordsRTemp(m,2) * workingUnitVector(2,2) + &
              coordsRTemp(m,3) * workingUnitVector(2,3)
!z-coord
              coordsRTemp(m,3) = & 		
			  coordsRTemp(m,1) * workingUnitVector(3,1) + &
              coordsRTemp(m,2) * workingUnitVector(3,2) + &
              coordsRTemp(m,3) * workingUnitVector(3,3)
!round coordinates
			  coordsRTemp(m,1) = rounding(coordsRTemp(m,1),6)
			  coordsRTemp(m,2) = rounding(coordsRTemp(m,2),6)
			  coordsRTemp(m,3) = rounding(coordsRTemp(m,3),6)
!Write to file
	          If(mpiProcessID.eq.0.and.saveFileCoords.eq."Y")Then
	            write(10,"(I8,I8,I8,F16.8,F16.8,F16.8)") i,m,coordsITemp(m),&
			    coordsRTemp(m,1),coordsRTemp(m,2),coordsRTemp(m,3)
	          End If
			  atomTypeKey(atomCounter) = configCoordsI(n)
!store atom-force data
              If(configForcesR(n,4).gt.0)Then
                configurationRefForce(atomCounter,1) = configForcesR(n,1)
                configurationRefForce(atomCounter,2) = configForcesR(n,2)
                configurationRefForce(atomCounter,3) = configForcesR(n,3)
			  End If
			enddo
		  enddo
		enddo
      enddo
!build neighbour list
      rCutoff = configHeaderR(i,2)
      rCutoffNeg = -1 * configHeaderR(i,2)
      rCutoffSq = configHeaderR(i,2)**2
	  xlat = alat * xCopy
	  ylat = alat * yCopy
	  zlat = alat * zCopy
!loop through atom pairs    Atom A in inner cell, Atom B in 3x3x3 supercell      
	  configLength = 0
!loop through Atom B 3x3x3
	  do l=-1,1
	    do m=-1,1
	      do n=-1,1	              			
!Set co-ordinate shift
			xshift = xlat * l
			yshift = ylat * m
			zshift = zlat * n
!Reset unique key list
			nlUniqueKeys = 0
            Do atomA=1,atoms
	          Do atomB=1,atoms 			  
!don't self count atom
			    if(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)then
			    else
				  If(atomA.lt.atomB)Then
				    nlKey = (atomB-1)*(atomB-2)/2.0D0+atomA
				  Else
				    nlKey = (atomA-1)*(atomA-2)/2.0D0+atomB
				  End If
				  If(nlUniqueKeys(nlKey).eq.0)Then
				    nlUniqueKeys(nlKey) = 1
!check range one co-ord at a time then distance squared
                    xA = 1.0D0*coordsRTemp(atomA,1)
                    xB = 1.0D0*(xshift + coordsRTemp(atomB,1))
                    yA = 1.0D0*coordsRTemp(atomA,2)
                    yB = 1.0D0*(yshift + coordsRTemp(atomB,2))
                    zA = 1.0D0*coordsRTemp(atomA,3)
                    zB = 1.0D0*(zshift + coordsRTemp(atomB,3))
                    xdSq = (xA-xB)**2
				    If(xdSq.le.rCutoffSq)then
				      ydSq = (yA-yB)**2
				      If(ydSq.le.rCutoffSq)then
				        zdSq = (zA-zB)**2
				        If(zdSq.le.rCutoffSq)then
					      rdSq = xdSq + ydSq + zdSq						
					      If(rdSq.le.rCutoffSq)then
					        neighbourListCount = neighbourListCount + 1
						    configLength = configLength + 1
!atom type data
						    neighbourListITemp(neighbourListCount,1) = coordsITemp(atomA) !Atom A Type
						    neighbourListITemp(neighbourListCount,2) = coordsITemp(atomB)	!Atom B Type
						    neighbourListITemp(neighbourListCount,5) = nlKey	!Pair key
						    If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
						      neighbourListITemp(neighbourListCount,6) = 1
						    Else
						      neighbourListITemp(neighbourListCount,6) = 0
						    End If
!atom id
						    neighbourListITemp(neighbourListCount,3) = atomA 	!Atom A ID
						    neighbourListITemp(neighbourListCount,4) = atomB 	!Atom B ID
!atom A and B co-ords (to recalculate rd if strain matrix applied)
						    neighbourListCoordsTemp(neighbourListCount,1) = xA	
						    neighbourListCoordsTemp(neighbourListCount,2) = yA
						    neighbourListCoordsTemp(neighbourListCount,3) = zA
						    neighbourListCoordsTemp(neighbourListCount,4) = xB
						    neighbourListCoordsTemp(neighbourListCount,5) = yB
						    neighbourListCoordsTemp(neighbourListCount,6) = zB
						    neighbourListRTemp(neighbourListCount) = rdSq**0.5						  
!Save neighbour list file
					      End If
					    End If
					  End If
			        End If
			      End If
				endif
			  enddo
			enddo
          enddo			
!end loop through Atom B 3x3x3	  
		enddo
	  enddo
!store configuration key information and update start point
      neighbourListKeyTemp(i,1) = configStart	  
      neighbourListKeyTemp(i,2) = configLength	  
	  configStart = configStart + configLength
	enddo
!end loop through configurations   
!Move from temp arrays - Allocate arrays
	Allocate(neighbourListKey(1:size(configHeaderI)/headerWidth,1:3))
	Allocate(neighbourListI(1:neighbourListCount,1:6))
	Allocate(neighbourListR(1:neighbourListCount)) 
	Allocate(neighbourListCoords(1:neighbourListCount,1:12)) !Ax Ay Az, Bx By Bz, R_ABi R_ABj R_ABk 
!Transfer key data 
    !do i=1,size(configHeaderI)/headerWidth
    do i=1,configCount
	  neighbourListKey(i,1) = neighbourListKeyTemp(i,1)	!Start
	  neighbourListKey(i,2) = neighbourListKeyTemp(i,2) !Length
	enddo
!Transfer neighbour list data	
    Do i=1,neighbourListCount
	  neighbourListI(i,1) = neighbourListITemp(i,1)  !Atom A Type
	  neighbourListI(i,2) = neighbourListITemp(i,2)  !Atom B Type
	  neighbourListI(i,3) = neighbourListITemp(i,3)  !Atom A ID
	  neighbourListI(i,4) = neighbourListITemp(i,4)  !Atom B ID
	  neighbourListI(i,5) = neighbourListITemp(i,5)  !nl key
	  neighbourListI(i,6) = neighbourListITemp(i,6)  !atom b in inner vol
	  neighbourListR(i) = neighbourListRTemp(i) !Distance/separation
	  neighbourListCoords(i,1) = neighbourListCoordsTemp(i,1) !xA
	  neighbourListCoords(i,2) = neighbourListCoordsTemp(i,2) !yA
	  neighbourListCoords(i,3) = neighbourListCoordsTemp(i,3) !zA
	  neighbourListCoords(i,4) = neighbourListCoordsTemp(i,4) !xB
	  neighbourListCoords(i,5) = neighbourListCoordsTemp(i,5) !yB
	  neighbourListCoords(i,6) = neighbourListCoordsTemp(i,6) !zB 
	  neighbourListCoords(i,7) = neighbourListCoordsTemp(i,4)-neighbourListCoordsTemp(i,1) !xB-xA
	  neighbourListCoords(i,8) = neighbourListCoordsTemp(i,5)-neighbourListCoordsTemp(i,2) !yA-yB
	  neighbourListCoords(i,9) = neighbourListCoordsTemp(i,6)-neighbourListCoordsTemp(i,3) !zA-zB
      neighbourListCoords(i,10) = neighbourListCoords(i,7)/neighbourListR(i)
      neighbourListCoords(i,11) = neighbourListCoords(i,8)/neighbourListR(i)
      neighbourListCoords(i,12) = neighbourListCoords(i,9)/neighbourListR(i)
!Round	  
	  neighbourListR(i) = Rounding(neighbourListR(i),8)
	  Do j=1,12
	    neighbourListCoords(i,j) = Rounding(neighbourListCoords(i,j),8)
	  End Do
	End Do		
!write neighbour list
    If(mpiProcessID.eq.0.and.saveFileNeighbourList.eq."Y")Then
      Do i=1,configCount
	    Do j=neighbourListKey(i,1),neighbourListKey(i,1)+neighbourListKey(i,2)-1
		  write(11,"(I8,I8,I8,I8,I8,I8,I8,I8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)") i,j,&
		  neighbourListI(j,1),neighbourListI(j,2),neighbourListI(j,3),neighbourListI(j,4),&
		  neighbourListI(j,5),neighbourListI(j,6),&
		  neighbourListCoords(j,1),neighbourListCoords(j,2),neighbourListCoords(j,3),&
		  neighbourListCoords(j,4),neighbourListCoords(j,5),neighbourListCoords(j,6),&
		  neighbourListR(j)
		End Do
	  End Do
	End If
!output to file	
    If(mpiProcessID.eq.0)Then
      write(999,"(A6,A31,I8)") "      ","Total configurations prepared: ",configCount
	  write(999,"(F8.4,A2,A53,I8)") ProgramTime(),"  ",&
	  "Finish building config neighbour lists, total pairs: ",neighbourListCount
	End If
!Close files
    close(999)
	If(mpiProcessID.eq.0.and.saveFileCoords.eq."Y")Then
	  close(10)		
	End If
	If(mpiProcessID.eq.0.and.saveFileNeighbourList.eq."Y")Then
	  close(11)	
	End If		
  End Subroutine makeNeighbourList    
  
!-------------------------- 
! Store neighbour list for recall
!--------------------------
  Subroutine storeNeighboutList()
!force declaration of all variables
	Implicit None		
!Deallocate arrays
	If(Allocated(neighbourListKeyStore))Then
	  Deallocate(neighbourListKeyStore)
	End If
	If(Allocated(neighbourListIStore))Then
	  Deallocate(neighbourListIStore)
	End If
	If(Allocated(neighbourListRStore))Then
	  Deallocate(neighbourListRStore)
	End If
	If(Allocated(neighbourListCoordsStore))Then
	  Deallocate(neighbourListCoordsStore)
	End If
!Store neighbour list arrays	
	neighbourListKeyStore = neighbourListKey
	neighbourListIStore = neighbourListI
	neighbourListRStore = neighbourListR
	neighbourListCoordsStore = neighbourListCoords
  End Subroutine storeNeighboutList  
  Subroutine loadNeighboutList()
!force declaration of all variables
	Implicit None		
!Deallocate arrays
	If(Allocated(neighbourListKey))Then
	  Deallocate(neighbourListKey)
	End If
	If(Allocated(neighbourListI))Then
	  Deallocate(neighbourListI)
	End If
	If(Allocated(neighbourListR))Then
	  Deallocate(neighbourListR)
	End If
	If(Allocated(neighbourListCoords))Then
	  Deallocate(neighbourListCoords)
	End If
!Store neighbour list arrays	
	neighbourListKey = neighbourListKeyStore
	neighbourListI = neighbourListIStore
	neighbourListR = neighbourListRStore
	neighbourListCoords = neighbourListCoordsStore
  End Subroutine loadNeighboutList 
!------------------------------------------------------------------------!
! ***calcConfigurationVolume*** 
!------------------------------------------------------------------------!    
  Subroutine applyDistortionVector(configurationID,distortionUnitVector)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: distortionUnitVector
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger) :: configStart, configLength
	Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, rD		
	configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)	
!Loop through neighbour list for the configuration
	Do i=configStart,configStart+configLength-1 
      xA = neighbourListCoords(i,1) * distortionUnitVector(1,1) + &
	  neighbourListCoords(i,2) * distortionUnitVector(2,1) + &
	  neighbourListCoords(i,3) * distortionUnitVector(3,1)
	  yA = neighbourListCoords(i,1) * distortionUnitVector(1,2) + &
	  neighbourListCoords(i,2) * distortionUnitVector(2,2) + &
	  neighbourListCoords(i,3) * distortionUnitVector(3,2)
	  zA = neighbourListCoords(i,1) * distortionUnitVector(1,3) + &
	  neighbourListCoords(i,2) * distortionUnitVector(2,3) + &
	  neighbourListCoords(i,3) * distortionUnitVector(3,3)
	  xB = neighbourListCoords(i,4) * distortionUnitVector(1,1) + &
	  neighbourListCoords(i,5) * distortionUnitVector(2,1) + &
	  neighbourListCoords(i,6) * distortionUnitVector(3,1)
	  yB = neighbourListCoords(i,4) * distortionUnitVector(1,2) + &
	  neighbourListCoords(i,5) * distortionUnitVector(2,2) + &
	  neighbourListCoords(i,6) * distortionUnitVector(3,2)
	  zB = neighbourListCoords(i,4) * distortionUnitVector(1,3) + &
	  neighbourListCoords(i,5) * distortionUnitVector(2,3) + &
	  neighbourListCoords(i,6) * distortionUnitVector(3,3)	  
!Update co-ordinate values
	  neighbourListCoords(i,1) = xA !xA
	  neighbourListCoords(i,2) = yA !yA
	  neighbourListCoords(i,3) = zA !zA
	  neighbourListCoords(i,4) = xB !xB
	  neighbourListCoords(i,5) = yB !yB
	  neighbourListCoords(i,6) = zB !zB 
	  neighbourListCoords(i,7) = xB-xA !xB-xA
	  neighbourListCoords(i,8) = yB-yA !yA-yB
	  neighbourListCoords(i,9) = zB-zA !zA-zB
	  rD = (neighbourListCoords(i,7)**2+&
	  neighbourListCoords(i,8)**2+neighbourListCoords(i,9)**2)**0.5	
	  !print *,rD,neighbourListR(i)
	  neighbourListR(i) = rD
      neighbourListCoords(i,10) = neighbourListCoords(i,7)/neighbourListR(i)
      neighbourListCoords(i,11) = neighbourListCoords(i,8)/neighbourListR(i)
      neighbourListCoords(i,12) = neighbourListCoords(i,9)/neighbourListR(i)	
    End Do	
  End Subroutine applyDistortionVector
!------------------------------------------------------------------------!
! applyUnitVector 
!------------------------------------------------------------------------!
  Subroutine applyUnitVector(configurationID)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k				
!loop through configurations
    if(configurationID.eq.0)then
	  do i=1,configCount
!write to output file
        Call applyUnitVectorAction(i)
	  enddo
	else
	  Call applyUnitVectorAction(configurationID)
	endif	
  End Subroutine applyUnitVector
!------------------------------------------------------------------------!
! applyUnitVectorAction 
!------------------------------------------------------------------------!
  Subroutine applyUnitVectorAction(configurationID)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k
	Integer(kind=StandardInteger) :: configStart, configLength
	Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, rD		
	configStart = neighbourListKey(configurationID,1)
    configLength = neighbourListKey(configurationID,2)	
!Loop through neighbour list for the configuration
	do i=configStart,configStart+configLength-1 
!calculate atom positions with unitVector applied
	  xA = neighbourListCoords(i,1) * unitVector(1,1) + &
	  neighbourListCoords(i,2) * unitVector(2,1) + &
	  neighbourListCoords(i,3) * unitVector(3,1)
	  yA = neighbourListCoords(i,1) * unitVector(1,2) + &
	  neighbourListCoords(i,2) * unitVector(2,2) + &
	  neighbourListCoords(i,3) * unitVector(3,2)
	  zA = neighbourListCoords(i,1) * unitVector(1,3) + &
	  neighbourListCoords(i,2) * unitVector(2,3) + &
	  neighbourListCoords(i,3) * unitVector(3,3)
	  xB = neighbourListCoords(i,4) * unitVector(1,1) + &
	  neighbourListCoords(i,5) * unitVector(2,1) + &
	  neighbourListCoords(i,6) * unitVector(3,1)
	  yB = neighbourListCoords(i,4) * unitVector(1,2) + &
	  neighbourListCoords(i,5) * unitVector(2,2) + &
	  neighbourListCoords(i,6) * unitVector(3,2)
	  zB = neighbourListCoords(i,4) * unitVector(1,3) + &
	  neighbourListCoords(i,5) * unitVector(2,3) + &
	  neighbourListCoords(i,6) * unitVector(3,3)	  	  
!store new displacement vector between atom A and B
	  neighbourListCoords(i,7) = xB-xA
	  neighbourListCoords(i,8) = yB-yA
	  neighbourListCoords(i,9) = zB-zA
!store new distance between atoms	
	  neighbourListR(i) = ((xA-xB)**2+(yA-yB)**2+(zA-zB)**2)**0.5	  
	enddo	
!recalculate configuration volume
    !Call calcConfigurationVolume(configurationID)	
  End Subroutine applyUnitVectorAction
!------------------------------------------------------------------------!
! ***calcConfigurationVolumes*** 
!------------------------------------------------------------------------!
  Subroutine calcConfigurationVolumes()
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i
!Calculate and store volumes 
    Do i=1,configCount
	  configurationVolume(i) = 1.0D0* calcConfigurationVolume(i)
	End Do  
  End Subroutine calcConfigurationVolumes
!------------------------------------------------------------------------!
! ***calcConfigurationVolume*** 
!------------------------------------------------------------------------!
  Function calcConfigurationVolume(configurationID, distortionArrayIn) result (volume)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k, n
    Real(kind=SingleReal) :: xA, yA, zA
    Real(kind=SingleReal) :: xB, yB, zB
    Real(kind=SingleReal) :: xC, yC, zC
	Real(kind=SingleReal) :: crossProductI,crossProductJ,crossProductK
	Real(kind=DoubleReal) :: volume
	Real(kind=DoubleReal), Optional, Dimension( : , : ), Allocatable :: distortionArrayIn
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: multiplierArray, distortionArray
!Set distortion array
    Allocate(multiplierArray(1:3,1:3))
!Set the configuration unit vector
    i = configurationID
	multiplierArray(1,1) = 1.0D0*configurationUnitVector((i-1)*9+1)
	multiplierArray(1,2) = 1.0D0*configurationUnitVector((i-1)*9+2)
	multiplierArray(1,3) = 1.0D0*configurationUnitVector((i-1)*9+3)
	multiplierArray(2,1) = 1.0D0*configurationUnitVector((i-1)*9+4)
	multiplierArray(2,2) = 1.0D0*configurationUnitVector((i-1)*9+5)
	multiplierArray(2,3) = 1.0D0*configurationUnitVector((i-1)*9+6)
	multiplierArray(3,1) = 1.0D0*configurationUnitVector((i-1)*9+7)
	multiplierArray(3,2) = 1.0D0*configurationUnitVector((i-1)*9+8)
	multiplierArray(3,3) = 1.0D0*configurationUnitVector((i-1)*9+9)
    If(present(distortionArrayIn))Then
	  distortionArray = distortionArrayIn
	  multiplierArray = matmul(multiplierArray,distortionArray)
	End If
!loop through configurations if 0 and calc volume for all configurations
    i = configurationID
	xA = configLatticeParameters(i,1) * multiplierArray(1,1)
    xB = configLatticeParameters(i,1) * multiplierArray(1,2)
    xC = configLatticeParameters(i,1) * multiplierArray(1,3)
    yA = configLatticeParameters(i,2) * multiplierArray(2,1)
    yB = configLatticeParameters(i,2) * multiplierArray(2,2)
    yC = configLatticeParameters(i,2) * multiplierArray(2,3)
    zA = configLatticeParameters(i,3) * multiplierArray(3,1)
    zB = configLatticeParameters(i,3) * multiplierArray(3,2)
    zC = configLatticeParameters(i,3) * multiplierArray(3,3)
	crossProductI = (xB*yC-xC*yB)
	crossProductJ = (xC*yA-xA*yC)
	crossProductK = (xA*yB-xB*yA)
	volume = 1.0D0*zA*crossProductI+zB*crossProductJ+zC*crossProductK	
  End Function calcConfigurationVolume
  Subroutine summariseConfigurations()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i
	Character(2) :: refForce
!If master process
	If(mpiProcessID.eq.0)Then
!open output file	
	  open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
	  status="old",position="append",action="write")	
!save summary of configurations to output file
      write(999,"(A6,A40)") "      ","----------------------------------------"
      write(999,"(A6,A25)") "      ","Summary of Configurations"
      write(999,"(A6,A40)") "      ","----------------------------------------"
      Do i=1,configCount
	    If(configAtomsMap(i,4).eq.1)Then
          refForce = "RF"
        Else
          refForce = "  "
        End If		
	  
        write(999,"(A6,A15,I8,A12,F14.6,A11,I8,A15,I8,A3,I8,I8,A2,A2,A2,I8)") &
		"      ",&
	    "Configuration: ",i,"    Volume: ",configurationVolume(i),&
	    "    Atoms: ",configAtoms(i),"    NL Length: ",neighbourListKey(i,2),"   ",&
		configAtomsMap(i,1),configAtomsMap(i,2),"  ",refForce,"  ",configAtomsMap(i,21)
	  End Do	
!close the output file
      close(999) 	
	End If
  End Subroutine summariseConfigurations    
!------------------------------------------------------------------------!
!                                                                        
! Potential Preparation Subroutines                                              
!                                                                        
!------------------------------------------------------------------------!  
!Order eam potentials
  Subroutine orderEamPotentials()  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, sortLoop, sorted, reordered
    Integer(kind=StandardInteger) :: potStart, potLength
	Real(kind=SingleReal) :: tempXA, tempYA, tempXB, tempYB	
!write if master process
    If(mpiProcessID.eq.0)Then	
!open output file	
	  open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
	  status="old",position="append",action="write")	
!write to output file
	  write(999,"(F8.4,A2,A20)") ProgramTime(),"  ","Order EAM Potentials"	
    End If
	reordered = 0
	do i=1,size(eamKey,1)
	  potStart = eamKey(i,4)
	  potLength = eamKey(i,5)
	  If(mpiProcessID.eq.0)Then
	    write(999,"(I8,I8)") potStart, potLength	
	  End If
	  sorted = 1
	  do while(sorted.eq.1)
	    sorted = 0
	    do j=potStart,potStart+potLength-2	    
	      if(eamData(j,1).gt.eamData(j+1,1))then
	        sorted = 1
			reordered = reordered + 1
!temporarily store variables
			tempXA = eamData(j,1)
			tempYA = eamData(j,2)
			tempXB = eamData(j+1,1)
			tempYB = eamData(j+1,2)
!move variables
			eamData(j,1) = tempXB
			eamData(j,2) = tempYB
			eamData(j+1,1) = tempXA
			eamData(j+1,2) = tempYA
	      endif
	    enddo
	  enddo		
	enddo	
	If(mpiProcessID.eq.0)Then
!close output file
      close(999)	
    End If	  
  End Subroutine orderEamPotentials

  
  
  
  
  Subroutine storeEAM()
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k  
!store eam data and close file
    if(saveFilePot.eq."Y".and.mpiProcessID.eq.0)then
	  do i=1,size(eamKey,1)
	    if(eamKey(i,3).eq.1)then
	      write(21,"(A5,A2,A1,A2,A4,I4,I4)") "PAIR ",elements(eamKey(i,1)),&
		  " ",elements(eamKey(i,2)),"    ",eamKey(i,1),eamKey(i,2)
		endif  
		if(eamKey(i,3).eq.2)then
		  if(eamType.eq.1)then
	        write(21,"(A5,A2,A4,I4)") "DENS ",elements(eamKey(i,1)),"    ",&
			eamKey(i,1)
		  endif
		  if(eamType.eq.2)then
	        write(21,"(A5,A2,A1,A2,A4,I4,I4)") "SDEN ",elements(eamKey(i,1))," ",&
			elements(eamKey(i,2)),"    ",eamKey(i,1),eamKey(i,2)
		  endif
		endif 
		If(eamKey(i,3).eq.3)then
		  If(eamType.eq.1)then
	        write(21,"(A5,A2,A4,I4)") "EMBE ",elements(eamKey(i,1)),"    ",&
		    eamKey(i,1)
		  End If
		  If(eamType.eq.2)then
		    write(21,"(A5,A2,A4,I4)") "SEMB ",elements(eamKey(i,1)),"    ",&
		    eamKey(i,1)
		  End If
		endif 
		if(eamKey(i,4).eq.4)then
	      write(21,"(A5,A2,A4,I4)") "DEND ",elements(eamKey(i,1)),"    ",&
		  eamKey(i,1)
		endif 
		if(eamKey(i,4).eq.5)then
	      write(21,"(A5,A2,A4,I4)") "DEMB ",elements(eamKey(i,1)),"    ",&
		  eamKey(i,1)
		endif 		
		k = 1
		do j=eamKey(i,4),(eamKey(i,4)+eamKey(i,5)-1)
		  write(21,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") eamData(j,1),"  ",&
		  eamData(j,2),"  ",eamData(j,3),"    ",i,k,j
		  k = k + 1
		enddo
      enddo
      close(21)
    endif	
  End Subroutine storeEAM
  
  
!Reduced set of points from input eam
  Subroutine makeReducedEAMSet(eamKeyArrayIn,eamDataArrayIn,eamKeyArrayOut,eamDataArrayOut) 	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArrayIn 
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArrayOut
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArrayIn
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArrayOut
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalReducedDataPoints = 0
    Do potKey=1,size(eamKeyArrayIn,1)
!Get number of points to reduce to
	  If(eamKeyArrayIn(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKeyArrayIn(potKey,3).eq.2.or.eamKeyArrayIn(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKeyArrayIn(potKey,3).eq.3.or.eamKeyArrayIn(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If	 
      totalReducedDataPoints = totalReducedDataPoints + numberOfPoints
    End Do	
!Allocate arrays
    If(Allocated(eamKeyArrayOut))Then
	  Deallocate(eamKeyArrayOut)
	End If
    If(Allocated(eamDataArrayOut))Then
	  Deallocate(eamDataArrayOut)
	End If
    Allocate(eamKeyArrayOut(1:size(eamKeyArrayIn,1),1:5))
    Allocate(eamDataArrayOut(1:totalReducedDataPoints,1:3))
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyArrayIn,1)
!Get number of points to reduce to
	  If(eamKeyArrayIn(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKeyArrayIn(potKey,3).eq.2.or.eamKeyArrayIn(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKeyArrayIn(potKey,3).eq.3.or.eamKeyArrayIn(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If
!make reduces set of data points
	  potStart = eamKeyArrayIn(potKey,4)
      potLength = eamKeyArrayIn(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyArrayOut(potKey,1) = eamKeyArrayIn(potKey,1) !Atom i
      eamKeyArrayOut(potKey,2) = eamKeyArrayIn(potKey,2) !Atom j
	  eamKeyArrayOut(potKey,3) = eamKeyArrayIn(potKey,3) !Pot type
	  eamKeyArrayOut(potKey,4) = dataPointCounter
	  eamKeyArrayOut(potKey,5) = numberOfPoints
!loop through reduced data points
	  Do i=1,numberOfPoints
	    x = eamDataArrayIn(potStart,1)+&
		    1.0D0*(i-1)*((eamDataArrayIn(potEnd,1)-eamDataArrayIn(potStart,1))/(numberOfPoints-1))
		!yArray = PointInterpolationArr(eamData,x,5,potStart,potLength,"N")
		yArray = PointInterpolationArr(eamDataArrayIn,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
!store data
        eamDataArrayOut(dataPointCounter,1) = x
        eamDataArrayOut(dataPointCounter,2) = y
        eamDataArrayOut(dataPointCounter,3) = dy
!increment counter
		dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
  End Subroutine makeReducedEAMSet
  
!Make expanded set of points (expanded from the reduced set, using 2nd order spline)
  Subroutine makeExpandedEAMSet(eamKeyArrayIn,eamDataArrayIn,eamKeyArrayOut,eamDataArrayOut)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalTrialDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splineXY
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArrayIn 
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArrayOut
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArrayIn
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArrayOut
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalTrialDataPoints = 0
    Do potKey=1,size(eamKeyArrayIn,1)
!Get number of points to reduce to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF
	  End If	 
      totalTrialDataPoints = totalTrialDataPoints + numberOfPoints
    End Do
!Allocate arrays
    If(Allocated(eamKeyArrayOut))Then
	  Deallocate(eamKeyArrayOut)
	End If
    If(Allocated(eamDataArrayOut))Then
	  Deallocate(eamDataArrayOut)
	End If
    Allocate(eamKeyArrayOut(1:size(eamKeyArrayIn,1),1:5))
    Allocate(eamDataArrayOut(1:totalTrialDataPoints,1:3))
!
! Spline nodes with matching zero, first and second order derivatives (6 points, 5th order poly)
!
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyArrayIn,1)
!Get number of points to expand to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF
	  End If
!make expanded set of data points
	  potStart = eamKeyArrayIn(potKey,4)
      potLength = eamKeyArrayIn(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyArrayOut(potKey,1) = eamKeyArrayIn(potKey,1) !Atom i
      eamKeyArrayOut(potKey,2) = eamKeyArrayIn(potKey,2) !Atom j
	  eamKeyArrayOut(potKey,3) = eamKeyArrayIn(potKey,3) !Pot type		
	  eamKeyArrayOut(potKey,4) = dataPointCounter
	  eamKeyArrayOut(potKey,5) = numberOfPoints
!loop through reduced data points
      splineXY = Spline(eamDataArrayIn,numberOfPoints,potStart,potLength)
	  Do i=1,numberOfPoints
!store data
        eamDataArrayOut(dataPointCounter,1) = splineXY(i,1)
        eamDataArrayOut(dataPointCounter,2) = splineXY(i,2)
        eamDataArrayOut(dataPointCounter,3) = splineXY(i,3)
!increment counter
	    dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
  End Subroutine makeExpandedEAMSet
  
  
  
  
  
  Subroutine storeReducedToOpt()  	
!force declaration of all variables
	Implicit None
!declare private variables  
    If(Allocated(eamKeyOpt))Then
	  Deallocate(eamKeyOpt)
	End If
    If(Allocated(eamDataOpt))Then
	  Deallocate(eamDataOpt)
	End If
    eamKeyOpt = eamKeyReduced  
    eamDataOpt = eamDataReduced
  End Subroutine storeReducedToOpt
  
  
  
  Subroutine loadOptToReduced()  	
!force declaration of all variables
	Implicit None
!declare private variables  
    If(Allocated(eamKeyReduced))Then
	  Deallocate(eamKeyReduced)
	End If
    If(Allocated(eamDataReduced))Then
	  Deallocate(eamDataReduced)
	End If
    eamKeyReduced = eamKeyOpt
    eamDataReduced = eamDataOpt
  End Subroutine loadOptToReduced
    
  Subroutine storeToReduced(eamKeyInput,eamDataInput)  	
!force declaration of all variables
	Implicit None	
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyInput
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataInput
!declare private variables  
    If(Allocated(eamKeyReduced))Then
	  Deallocate(eamKeyReduced)
	End If
    If(Allocated(eamDataReduced))Then
	  Deallocate(eamDataReduced)
	End If
    eamKeyReduced = eamKeyInput
    eamDataReduced = eamDataInput
  End Subroutine storeToReduced
  
  Subroutine clearReducedEAM()  	
!force declaration of all variables
	Implicit None	
!declare private variables  
    If(Allocated(eamKeyReduced))Then
	  Deallocate(eamKeyReduced)
	End If
    If(Allocated(eamDataReduced))Then
	  Deallocate(eamDataReduced)
	End If
  End Subroutine clearReducedEAM
  
  Subroutine clearOptEAM()  	
!force declaration of all variables
	Implicit None	
!declare private variables  
    If(Allocated(eamKeyOpt))Then
	  Deallocate(eamKeyOpt)
	End If
    If(Allocated(eamDataOpt))Then
	  Deallocate(eamDataOpt)
	End If
  End Subroutine clearOptEAM
  
  Subroutine clearTrialEAM()  	
!force declaration of all variables
	Implicit None	
!declare private variables  
    If(Allocated(eamKeyTrial))Then
	  Deallocate(eamKeyTrial)
	End If
    If(Allocated(eamDataTrial))Then
	  Deallocate(eamDataTrial)
	End If
  End Subroutine clearTrialEAM
  
  Subroutine clearEAM()  	
!force declaration of all variables
	Implicit None	
!declare private variables  
    If(Allocated(eamKey))Then
	  Deallocate(eamKey)
	End If
    If(Allocated(eamData))Then
	  Deallocate(eamData)
	End If
  End Subroutine clearEAM
  
  
  
  Subroutine setPotentialDerivatives(eamKeyArray,eamDataArray)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
!Loop through potential functions
	Do potKey=1,size(eamKeyArray,1)
!make reduces set of data points
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
!loop over data points
	  Do i=potStart,potEnd
	    x = eamDataArray(i,1)
		yArray = PointInterpolationArr(eamDataArray,x,4,potStart,potLength) 
		eamDataArray(i,3) = yArray(2)		  
	  End Do	
    End Do  
  End Subroutine setPotentialDerivatives
  

  Subroutine storeReducedToFile()  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
!Allocate arrays	
!Store reduced potential to file
	if(mpiProcessID.eq.0)then
	  outputFile = trim(currentWorkingDirectory)//"/"//"outputreducedtemp.pot"
	  open(unit=24,file=trim(outputFile))
	  Do potKey=1,size(eamKeyReduced,1)
!Get number of points to reduce to
	    If(eamKeyReduced(potKey,3).eq.1)Then     !Pair
	      write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyReduced(potKey,1))," ",&
		  elements(eamKeyReduced(potKey,2))," "
		Else  
		  If(eamKeyReduced(potKey,3).eq.2.and.eamType.eq.1)Then
		    eamTypeText = "DENS "
		  End If
		  If(eamKeyReduced(potKey,3).eq.2.and.eamType.eq.2)Then
		    eamTypeText = "SDEN "
		  End If
		  If(eamKeyReduced(potKey,3).eq.3.and.eamType.eq.1)Then
		    eamTypeText = "EMBE "
		  End If
		  If(eamKeyReduced(potKey,3).eq.3.and.eamType.eq.2)Then
		    eamTypeText = "SEMB "
		  End If
		  If(eamKeyReduced(potKey,3).eq.4)Then
		    eamTypeText = "DDEN "
		  End If
		  If(eamKeyReduced(potKey,3).eq.5)Then
		    eamTypeText = "DEMB "
		  End If
		  write(24,"(A5,A2)") eamTypeText,elements(eamKeyReduced(potKey,1))
	    End If	
!make reduces set of data points
	    potStart = eamKeyReduced(potKey,4)
        potLength = eamKeyReduced(potKey,5)
	    potEnd = potStart + potLength - 1
!loop over data points
        k = 0
	    Do i=potStart,potEnd
		  k = k + 1
	      x = eamDataReduced(i,1)
		  y = eamDataReduced(i,2)
		  dy = eamDataReduced(i,3)
		  dy = 0.0D0
!write to file
          write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		  y,"  ",dy,"    ",potKey,k,i
	    End Do	
      End Do
!Close file
	  close(24)
	End If
  End Subroutine storeReducedToFile  
  
  Subroutine storeEAMToFile(eamKeyArray, eamDataArray, fileName, numberOfPointsIn, processIn)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints, processFlag
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	Character(*) :: fileName
!optional variables	
	Integer(kind=StandardInteger), optional :: numberOfPointsIn
	Integer(kind=StandardInteger), optional :: processIn
	  If(Present(numberOfPointsIn))Then
	    numberOfPoints = numberOfPointsIn
	  Else
        numberOfPoints = 0	  
	  End If
	  If(Present(processIn))Then
	    processFlag = processIn
	  Else
	    processFlag = -1
	  End If	
!Store reduced potential to file
	If(mpiProcessID.eq.0)Then
	  open(unit=24,file=trim(fileName))
	  If(numberOfPoints.lt.10)Then
	    Do potKey=1,size(eamKeyArray,1)
!Get number of points to reduce to
	      If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	        write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyArray(potKey,1))," ",&
		    elements(eamKeyArray(potKey,2))," "
		  Else  
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		      eamTypeText = "DENS "
		    End If
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		      eamTypeText = "SDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		      eamTypeText = "EMBE "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		      eamTypeText = "SEMB "
		    End If
		    If(eamKeyArray(potKey,3).eq.4)Then
		      eamTypeText = "DDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.5)Then
		      eamTypeText = "DEMB "
		    End If
		    write(24,"(A5,A2)") eamTypeText,elements(eamKeyArray(potKey,1))
	      End If	
!pot positions
	      potStart = eamKeyArray(potKey,4)
          potLength = eamKeyArray(potKey,5)
	      potEnd = potStart + potLength - 1
!loop over data points
          k = 0
	      Do i=potStart,potEnd
		    k = k + 1
	        x = eamDataArray(i,1)
		    y = eamDataArray(i,2)
		    dy = eamDataArray(i,3)
!write to file
            write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		    y,"  ",dy,"    ",potKey,k,i
	      End Do	
        End Do
	  ElseIf(numberOfPoints.ge.10)Then
!Interpolate between points
!Loop through potential functions
        dataPointCounter = 1
	    Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	      If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	        write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyArray(potKey,1))," ",&
		    elements(eamKeyArray(potKey,2))," "
		  Else  
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		      eamTypeText = "DENS "
		    End If
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		      eamTypeText = "SDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		      eamTypeText = "EMBE "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		      eamTypeText = "SEMB "
		    End If
		    If(eamKeyArray(potKey,3).eq.4)Then
		      eamTypeText = "DDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.5)Then
		      eamTypeText = "DEMB "
		    End If
		    write(24,"(A5,A2)") eamTypeText,elements(eamKeyArray(potKey,1))
	      End If	
	      potStart = eamKeyArray(potKey,4)
          potLength = eamKeyArray(potKey,5)
	      potEnd = potStart + potLength - 1
!loop through reduced data points
!Interpolate 4th order polynomial (5 data points)
	      Do i=1,numberOfPoints
	        x = eamDataArray(potStart,1)+&
		    1.0D0*(i-1)*((eamDataArray(potEnd,1)-&
			eamDataArray(potStart,1))/(numberOfPoints-1))
		    !yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength,"N")
		    yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength)
		    y = yArray(1)
		    dy = yArray(2)
!write to file
            write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		    y,"  ",dy,"    ",potKey,i,dataPointCounter	
!increment counter
		    dataPointCounter = dataPointCounter + 1		
	      End Do	
        End Do
	  End If
!Close file
	  close(24)
	End If
!Wait for all processes to catch up	 
    If(processFlag.lt.0)Then
	  Call synchMpiProcesses()
	End If  
  End Subroutine storeEAMToFile  
  Subroutine storeEAMToFileMaster(eamKeyArray, eamDataArray, fileName, numberOfPointsIn)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints, processFlag
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	Character(*) :: fileName
!optional variables	
	Integer(kind=StandardInteger), optional :: numberOfPointsIn
	  If(Present(numberOfPointsIn))Then
	    numberOfPoints = numberOfPointsIn
	  Else
        numberOfPoints = 0	  
	  End If
!Open output file	  
	  outputFile = trim(currentWorkingDirectory)//"/"//trim(fileName)
	  open(unit=24,file=trim(outputFile))
	  If(numberOfPoints.lt.10)Then
	    Do potKey=1,size(eamKeyArray,1)
!Get number of points to reduce to
	      If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	        write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyArray(potKey,1))," ",&
		    elements(eamKeyArray(potKey,2))," "
		  Else  
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		      eamTypeText = "DENS "
		    End If
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		      eamTypeText = "SDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		      eamTypeText = "EMBE "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		      eamTypeText = "SEMB "
		    End If
		    If(eamKeyArray(potKey,3).eq.4)Then
		      eamTypeText = "DDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.5)Then
		      eamTypeText = "DEMB "
		    End If
		    write(24,"(A5,A2)") eamTypeText,elements(eamKeyArray(potKey,1))
	      End If	
!pot positions
	      potStart = eamKeyArray(potKey,4)
          potLength = eamKeyArray(potKey,5)
	      potEnd = potStart + potLength - 1
!loop over data points
          k = 0
	      Do i=potStart,potEnd
		    k = k + 1
	        x = eamDataArray(i,1)
		    y = eamDataArray(i,2)
		    dy = eamDataArray(i,3)
!write to file
            write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		    y,"  ",dy,"    ",potKey,k,i
	      End Do	
        End Do
	  ElseIf(numberOfPoints.ge.10)Then
!Interpolate between points
!Loop through potential functions
        dataPointCounter = 1
	    Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	      If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	        write(24,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyArray(potKey,1))," ",&
		    elements(eamKeyArray(potKey,2))," "
		  Else  
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		      eamTypeText = "DENS "
		    End If
		    If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		      eamTypeText = "SDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		      eamTypeText = "EMBE "
		    End If
		    If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		      eamTypeText = "SEMB "
		    End If
		    If(eamKeyArray(potKey,3).eq.4)Then
		      eamTypeText = "DDEN "
		    End If
		    If(eamKeyArray(potKey,3).eq.5)Then
		      eamTypeText = "DEMB "
		    End If
		    write(24,"(A5,A2)") eamTypeText,elements(eamKeyArray(potKey,1))
	      End If	
	      potStart = eamKeyArray(potKey,4)
          potLength = eamKeyArray(potKey,5)
	      potEnd = potStart + potLength - 1
!loop through reduced data points
!Interpolate 4th order polynomial (5 data points)
	      Do i=1,numberOfPoints
	        x = eamDataArray(potStart,1)+&
		    1.0D0*(i-1)*((eamDataArray(potEnd,1)-&
			eamDataArray(potStart,1))/(numberOfPoints-1))
		    !yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength,"N")
		    yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength)
		    y = yArray(1)
		    dy = yArray(2)
!write to file
            write(24,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		    y,"  ",dy,"    ",potKey,i,dataPointCounter	
!increment counter
		    dataPointCounter = dataPointCounter + 1		
	      End Do	
        End Do
	  End If
!Close file
	  close(24) 
  End Subroutine storeEAMToFileMaster  
  
  Subroutine printEAM(eamKeyArray, eamDataArray)
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
!Print
	If(mpiProcessID.eq.0)Then
	  Do potKey=1,size(eamKeyArray,1)
!Get number of points to reduce to
	    If(eamKeyArray(potKey,3).eq.1)Then     !Pair
	      print "(A5,A2,A1,A2,A1)","PAIR ",elements(eamKeyArray(potKey,1))," ",&
		  elements(eamKeyArray(potKey,2))," "
		Else  
		  If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.1)Then
		    eamTypeText = "DENS "
		  End If
		  If(eamKeyArray(potKey,3).eq.2.and.eamType.eq.2)Then
		    eamTypeText = "SDEN "
		  End If
		  If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.1)Then
		    eamTypeText = "EMBE "
		  End If
		  If(eamKeyArray(potKey,3).eq.3.and.eamType.eq.2)Then
		    eamTypeText = "SEMB "
		  End If
		  If(eamKeyArray(potKey,3).eq.4)Then
		    eamTypeText = "DDEN "
		  End If
		  If(eamKeyArray(potKey,3).eq.5)Then
		    eamTypeText = "DEMB "
		  End If
		  print "(A5,A2)",eamTypeText,elements(eamKeyArray(potKey,1))
	    End If	
!pot positions
	    potStart = eamKeyArray(potKey,4)
        potLength = eamKeyArray(potKey,5)
	    potEnd = potStart + potLength - 1
!loop over data points
        k = 0
	    Do i=potStart,potEnd
		  k = k + 1
	      x = eamDataArray(i,1)
		  y = eamDataArray(i,2)
		  dy = eamDataArray(i,3)
!write to file
          print "(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)",x,"  ",&
		  y,"  ",dy,"    ",potKey,k,i
	    End Do	
      End Do
	End If
  End Subroutine printEAM
  
!Convert input eam potential into a splined version
  Subroutine splinePotential(eamKeyWorking,eamDataWorking,pointToVaryIn,varyAmountIn) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints, totalSplineDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splineXY
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyWorking
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyReducedArray
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOut 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataWorking
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataReducedArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOut
!optional variables
	Integer(kind=StandardInteger), Optional :: pointToVaryIn
	Integer(kind=StandardInteger) :: pointToVary
	Real(kind=DoubleReal), Optional :: varyAmountIn
	Real(kind=DoubleReal) :: varyAmount	
!Set optional variables
    pointToVary = 0
	varyAmount = 0.0D0
    If(Present(pointToVaryIn))Then
	  pointToVary = pointToVaryIn
	End If		
    If(Present(varyAmountIn))Then
	  varyAmount = varyAmountIn
	End If		
!
! Reduced points of input potential
!		
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalReducedDataPoints = 0
    Do potKey=1,size(eamKeyWorking,1)
!Get number of points to reduce to
	  If(eamKeyWorking(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV	!Default 50
	  End If	
	  If(eamKeyWorking(potKey,3).eq.2.or.eamKeyWorking(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP	!Default 50
	  End If		
	  If(eamKeyWorking(potKey,3).eq.3.or.eamKeyWorking(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF	!Default 25
	  End If	 
      totalReducedDataPoints = totalReducedDataPoints + numberOfPoints
    End Do	
!Allocate arrays
    If(Allocated(eamKeyReducedArray))Then
	  Deallocate(eamKeyReducedArray)
	End If
    If(Allocated(eamDataReducedArray))Then
	  Deallocate(eamDataReducedArray)
	End If
    Allocate(eamKeyReducedArray(1:size(eamKeyWorking,1),1:5))
    Allocate(eamDataReducedArray(1:totalReducedDataPoints,1:3))
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyWorking,1)
!Get number of points to reduce to
	  If(eamKeyWorking(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKeyWorking(potKey,3).eq.2.or.eamKeyWorking(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKeyWorking(potKey,3).eq.3.or.eamKeyWorking(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If
!make reduces set of data points
	  potStart = eamKeyWorking(potKey,4)
      potLength = eamKeyWorking(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyReducedArray(potKey,1) = eamKeyWorking(potKey,1) !Atom i
      eamKeyReducedArray(potKey,2) = eamKeyWorking(potKey,2) !Atom j
	  eamKeyReducedArray(potKey,3) = eamKeyWorking(potKey,3) !Pot type
	  eamKeyReducedArray(potKey,4) = dataPointCounter
	  eamKeyReducedArray(potKey,5) = numberOfPoints
!loop through reduced data points
	  Do i=1,numberOfPoints
	    x = eamDataWorking(potStart,1)+&
		    1.0D0*(i-1)*((eamDataWorking(potEnd,1)-eamDataWorking(potStart,1))/(numberOfPoints-1))
		yArray = PointInterpolationArr(eamDataWorking,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
!store data
        eamDataReducedArray(dataPointCounter,1) = x
        eamDataReducedArray(dataPointCounter,2) = y
        eamDataReducedArray(dataPointCounter,3) = dy
!increment counter
		dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
!
! Vary point if required
!	
	If(pointToVary.gt.0.and.varyAmount.ne.0.0D0)Then
      pointToVary = mod(pointToVary,size(eamDataReducedArray,1))+1
	  eamDataReducedArray = VaryPoints(eamDataReducedArray,varyAmount,pointToVary)
    End If	
!
! Spline reduced points
!	
!Sum total spline points
    totalSplineDataPoints = 0
    Do potKey=1,size(eamKeyReducedArray,1)
!Get number of points to reduce to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV	!default 1001
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP	!default 1001
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF	!default 1001
	  End If	 
      totalSplineDataPoints = totalSplineDataPoints + numberOfPoints
    End Do	
!Allocate arrays
    If(Allocated(eamKeyWorking))Then
	  Deallocate(eamKeyWorking)
	End If
    If(Allocated(eamDataWorking))Then
	  Deallocate(eamDataWorking)
	End If
    Allocate(eamKeyWorking(1:size(eamKeyReducedArray,1),1:5))
    Allocate(eamDataWorking(1:totalSplineDataPoints,1:3))
!
! Spline nodes with matching zero, first and second order derivatives (6 points, 5th order poly)
!
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyReducedArray,1)
!Get number of points to expand to
	  If(eamKeyReducedArray(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV	!default 1001
	  End If	
	  If(eamKeyReducedArray(potKey,3).eq.2.or.eamKeyReducedArray(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP	!default 1001
	  End If		
	  If(eamKeyReducedArray(potKey,3).eq.3.or.eamKeyReducedArray(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF	!default 1001
	  End If
!make expanded set of data points
	  potStart = eamKeyReducedArray(potKey,4)
      potLength = eamKeyReducedArray(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyWorking(potKey,1) = eamKeyReducedArray(potKey,1) !Atom i
      eamKeyWorking(potKey,2) = eamKeyReducedArray(potKey,2) !Atom j
	  eamKeyWorking(potKey,3) = eamKeyReducedArray(potKey,3) !Pot type
	  eamKeyWorking(potKey,4) = dataPointCounter
	  eamKeyWorking(potKey,5) = numberOfPoints	  
!loop through reduced data points
      splineXY = Spline(eamDataReducedArray,numberOfPoints,potStart,potLength)
	  Do i=1,numberOfPoints
!store data	  
		eamDataWorking(dataPointCounter,1) = splineXY(i,1)
        eamDataWorking(dataPointCounter,2) = splineXY(i,2)
        eamDataWorking(dataPointCounter,3) = splineXY(i,3)	
!increment counter
	    dataPointCounter = dataPointCounter + 1
	  End Do
    End Do	
  End Subroutine splinePotential
  

  
  Subroutine eamForceZBLCore(eamKeyArray, eamDataArray)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey, potType
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints, processFlag
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd
	Integer(kind=StandardInteger) :: atomA, atomB
	Integer(kind=StandardInteger) :: zA, zB
	Real(kind=DoubleReal) :: x,y,dy
	Real(kind=DoubleReal) :: xLower, xUpper
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	Real(kind=DoubleReal), Dimension(1:3) :: yArray, yArrayA, yArrayB
	Real(kind=DoubleReal), Dimension(1:8) :: splinePoints
	Real(kind=DoubleReal), Dimension(0:5) :: pairSplineCoefficients, densitySplineCoefficients
	Real(kind=DoubleReal), Dimension(0:5) :: embeddingSplineCoefficients
	Do potKey=1,size(eamKeyArray,1)
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
	  !pot type
	  potType = eamKeyArray(potKey,3) !1 PAIR, 2 DENS/SDEN, 3 EMBE/SEMB, 4 DDEN, 5 DEMB
!atom details
	  atomA = eamKeyArray(potKey,1)
	  atomB = eamKeyArray(potKey,2)
	  zA = elementsCharge(atomA)
	  zB = elementsCharge(atomB)
!pot start and end
	  xLower = eamDataArray(potStart,1)
	  xUpper = eamDataArray(potEnd,1)
!----------------------------------------------------
!Pair Potential Spline start-end
!----------------------------------------------------
	  If(potType.eq.1.and.&
	  eamZBLPairLower.gt.-2.0D20.and.&
	  eamZBLPairUpper.gt.-2.0D20.and.&
	  eamZBLPairUpper.gt.eamZBLPairLower)Then
!Start point
        yArrayA = ZblFull(eamZBLPairLower,zA,zB)
!End point  
        yArrayB = PointInterpolationArr(eamDataArray,eamZBLPairUpper,5,potStart,potLength,2)
!Spline coefficients
        splinePoints(1) = 1.0D0*eamZBLPairLower 
        splinePoints(2) = 1.0D0*yArrayA(1) 
        splinePoints(3) = 1.0D0*yArrayA(2) 
        splinePoints(4) = 1.0D0*yArrayA(3)
        splinePoints(5) = 1.0D0*eamZBLPairUpper 
        splinePoints(6) = 1.0D0*yArrayB(1) 
        splinePoints(7) = 1.0D0*yArrayB(2) 
        splinePoints(8) = 1.0D0*yArrayB(3)  
!splinePoints
		pairSplineCoefficients = SplineAB(splinePoints) 
!loop over data points
	    Do i=potStart,potEnd
	      x = eamDataArray(i,1)
          If(x.le.eamZBLPairLower.and.potType.eq.1)Then
		    eamDataArray(i,2) = Zbl(x,zA,zB)
		  End If
          If(x.gt.eamZBLPairLower.and.x.lt.eamZBLPairUpper.and.potType.eq.1)Then
		    eamDataArray(i,2) = &
			computePolynomialA(pairSplineCoefficients,size(pairSplineCoefficients),x)
		  End If	
	    End Do	  
	  End If
!----------------------------------------------------
!Density Spline start-end
!----------------------------------------------------
	  If(potType.eq.2)Then
!End point  
        yArrayB = PointInterpolationArr(eamDataArray,eamZBLDensCutoff,5,potStart,potLength,2)
!Spline coefficients
        splinePoints(1) = 0.0D0
        splinePoints(2) = 1.0D0*eamZBLDensZero
        splinePoints(3) = 0.0D0 
        splinePoints(4) = 0.0D0
        splinePoints(5) = 1.0D0*eamZBLDensCutoff 
        splinePoints(6) = 1.0D0*yArrayB(1) 
        splinePoints(7) = 1.0D0*yArrayB(2) 
        splinePoints(8) = 1.0D0*yArrayB(3)  		  
!splinePoints
		densitySplineCoefficients = SplineAB(splinePoints) 
!loop over data points
	    Do i=potStart,potEnd
	      x = eamDataArray(i,1)
          If(x.lt.eamZBLDensCutoff.and.(potType.eq.2.or.potType.eq.4))Then
		    eamDataArray(i,2) = &
			computePolynomialA(densitySplineCoefficients,size(densitySplineCoefficients),x)
		  End If
	    End Do	  
      End If
!----------------------------------------------------
!Embedding Spline start-end
!----------------------------------------------------
	  If(potType.eq.3)Then
!End point  
        yArrayB = PointInterpolationArr(eamDataArray,eamZBLEmbeCutoff,5,potStart,potLength,2)
!Spline coefficients
        splinePoints(1) = 0.0D0
        splinePoints(2) = 1.0D0*eamZBLEmbeZero
        splinePoints(3) = 0.0D0 
        splinePoints(4) = 0.0D0
        splinePoints(5) = 1.0D0*eamZBLEmbeCutoff 
        splinePoints(6) = 1.0D0*yArrayB(1) 
        splinePoints(7) = 1.0D0*yArrayB(2) 
        splinePoints(8) = 1.0D0*yArrayB(3)  		  
!splinePoints
		embeddingSplineCoefficients = SplineAB(splinePoints) 
		Do i=potStart,potEnd
	      x = eamDataArray(i,1)
          If(x.lt.eamZBLEmbeCutoff.and.(potType.eq.3.or.potType.eq.5))Then
		    eamDataArray(i,2) = &
			computePolynomialA(embeddingSplineCoefficients,size(embeddingSplineCoefficients),x)
		  End If
	    End Do
      End If   
	End Do
  End Subroutine eamForceZBLCore  
  
  
  Subroutine storeReducedEAM(eamKeyArrayR,eamDataArrayR,fileName,headerIn)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArrayR 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArrayR
	Character(*) :: fileName
    Integer(kind=StandardInteger) :: send,receive,tag 
    Integer(kind=StandardInteger) :: processID,processCount,error,status 
!optional
    Character(64), optional :: headerIn
!call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)	
!Save to file   
    If(processID.eq.0)Then
	  outputFile = trim(currentWorkingDirectory)//"/"//trim(fileName)
	  open(unit=32,file=trim(outputFile))
!Write header
      If(Present(headerIn))Then
	    write(32,"(A64)") headerIn
	  End If
!Write key
      write(32,"(A9)") "<EAM Key>"
	  Do i=1,size(eamKeyArrayR,1)
	    write(32,"(I8,A2,I8,A2,I8,A2,I8,A2,I8)")&
		eamKeyArrayR(i,1),"  ",eamKeyArrayR(i,2),"  ",eamKeyArrayR(i,3),&
		"  ",eamKeyArrayR(i,4),"  ",eamKeyArrayR(i,5)
      End Do
      write(32,"(A10)") "</EAM Key>"
	  write(32,"(A10)") "<EAM Data>"
	  Do i=1,size(eamDataArrayR,1)
	    write(32,"(E24.16E3,A2,E24.16E3)") eamDataArrayR(i,1),"  ",eamDataArrayR(i,2)
	  End Do	
	  write(32,"(A11)") "</EAM Data>"
	  close(32)
	End If
!Synch Processes
	Call synchMpiProcesses()
  End Subroutine storeReducedEAM
   
  
  Subroutine loadReducedEAM(eamKeyArrayR,eamDataArrayR,fileName)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
    Integer(kind=StandardInteger) :: ios, i, j
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArrayR 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArrayR
	Character(*) :: fileName
	Character(len=255) :: fileRow
    Integer(kind=StandardInteger) :: send,receive,tag 
    Integer(kind=StandardInteger) :: processID,processCount,error,status 
!call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
!Read in file
	Open(unit=34,file=trim(trim(currentWorkingDirectory)//"/"//trim(fileName)))
    Do i=1,maxFileRows 
!Read in line
	  Read(34,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
	  If (ios /= 0) Then
	    EXIT 
	  End If 
!Read Key
      If(fileRow(1:9).eq."<EAM Key>")Then
	    Do j=1,maxFileRows
	      Read(34,"(A255)",IOSTAT=ios) fileRow
		  If(fileRow(1:2).eq."</".or.j.gt.size(eamKeyArrayR,1))Then
		    EXIT
		  End If
!read row
          Read(fileRow,*) eamKeyArrayR(j,1), eamKeyArrayR(j,2),&
		  eamKeyArrayR(j,3),eamKeyArrayR(j,4),eamKeyArrayR(j,5)
        End Do		
	  End If
!Read Data
      If(fileRow(1:10).eq."<EAM Data>")Then
	    Do j=1,maxFileRows
	      Read(34,"(A255)",IOSTAT=ios) fileRow
		  If(fileRow(1:2).eq."</".or.j.gt.size(eamDataArrayR,1))Then
		    EXIT
		  End If
!read row
          Read(fileRow,*) eamDataArrayR(j,1), eamDataArrayR(j,2)
	      eamDataArrayR(j,1) = 1.0D0 * eamDataArrayR(j,1)
	      eamDataArrayR(j,2) = 1.0D0 * eamDataArrayR(j,2)
        End Do		
	  End If
	  
    End Do
	CLOSE(34) !close file
!Synch Processes
	Call synchMpiProcesses()
  End Subroutine loadReducedEAM
  
  
  !Public :: storeReducedEAM
  !Public :: loadReducedEAM
  
  
  
  
!Misc Functions/Subroutines
  
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !------------------------------------------------------------------------------------------
  !
  !  Old functions
  !
  !------------------------------------------------------------------------------------------  
  
  
  
  
!Convert input eam potential into a splined version
  Subroutine splinePotentialOld(eamKeyIn,eamDataIn,eamKeyOut,eamDataOut,pointToVaryIn,varyAmountIn) 
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints, totalSplineDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splineXY
	Character(len=5)  :: eamTypeText
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyIn 
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyReducedArray
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyOut 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataIn
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataReducedArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataOut
!optional variables
	Integer(kind=StandardInteger), Optional :: pointToVaryIn
	Integer(kind=StandardInteger) :: pointToVary
	Real(kind=DoubleReal), Optional :: varyAmountIn
	Real(kind=DoubleReal) :: varyAmount	
!Set optional variables
    pointToVary = 0
	varyAmount = 0.0D0
    If(Present(pointToVaryIn))Then
	  pointToVary = pointToVaryIn
	End If		
    If(Present(varyAmountIn))Then
	  varyAmount = varyAmountIn
	End If		
!
! Reduced points of input potential
!		
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalReducedDataPoints = 0
    Do potKey=1,size(eamKeyIn,1)
!Get number of points to reduce to
	  If(eamKeyIn(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV	!Default 50
	  End If	
	  If(eamKeyIn(potKey,3).eq.2.or.eamKeyIn(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP	!Default 50
	  End If		
	  If(eamKeyIn(potKey,3).eq.3.or.eamKeyIn(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF	!Default 25
	  End If	 
      totalReducedDataPoints = totalReducedDataPoints + numberOfPoints
    End Do	
!Allocate arrays
    If(Allocated(eamKeyReducedArray))Then
	  Deallocate(eamKeyReducedArray)
	End If
    If(Allocated(eamDataReducedArray))Then
	  Deallocate(eamDataReducedArray)
	End If
    Allocate(eamKeyReducedArray(1:size(eamKeyIn,1),1:5))
    Allocate(eamDataReducedArray(1:totalReducedDataPoints,1:3))
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyIn,1)
!Get number of points to reduce to
	  If(eamKeyIn(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKeyIn(potKey,3).eq.2.or.eamKeyIn(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKeyIn(potKey,3).eq.3.or.eamKeyIn(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If
!make reduces set of data points
	  potStart = eamKeyIn(potKey,4)
      potLength = eamKeyIn(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyReducedArray(potKey,1) = eamKeyIn(potKey,1) !Atom i
      eamKeyReducedArray(potKey,2) = eamKeyIn(potKey,2) !Atom j
	  eamKeyReducedArray(potKey,3) = eamKeyIn(potKey,3) !Pot type
	  eamKeyReducedArray(potKey,4) = dataPointCounter
	  eamKeyReducedArray(potKey,5) = numberOfPoints
!loop through reduced data points
	  Do i=1,numberOfPoints
	    x = eamDataIn(potStart,1)+&
		    1.0D0*(i-1)*((eamDataIn(potEnd,1)-eamDataIn(potStart,1))/(numberOfPoints-1))
		yArray = PointInterpolationArr(eamDataIn,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
!store data
        eamDataReducedArray(dataPointCounter,1) = x
        eamDataReducedArray(dataPointCounter,2) = y
        eamDataReducedArray(dataPointCounter,3) = dy
!increment counter
		dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
!
! Vary point if required
!	
	If(pointToVary.gt.0.and.varyAmount.ne.0.0D0)Then
      pointToVary = mod(pointToVary,size(eamDataReducedArray,1))+1
	  If(mpiProcessID.eq.0)Then
	    print *,varyAmount,pointToVary
	    print *,eamDataReducedArray(pointToVary,1),eamDataReducedArray(pointToVary,2)
	  End If
	  eamDataReducedArray = VaryPoints(eamDataReducedArray,varyAmount,pointToVary)
	  If(mpiProcessID.eq.0)Then
	    print *,eamDataReducedArray(pointToVary,1),eamDataReducedArray(pointToVary,2)
	  End If
    End If
!
! Spline reduced points
!	
!Sum total spline points
    totalSplineDataPoints = 0
    Do potKey=1,size(eamKeyReducedArray,1)
!Get number of points to reduce to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV	!default 1001
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP	!default 1001
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF	!default 1001
	  End If	 
      totalSplineDataPoints = totalSplineDataPoints + numberOfPoints
    End Do
!Allocate arrays
    If(Allocated(eamKeyOut))Then
	  Deallocate(eamKeyOut)
	End If
    If(Allocated(eamDataOut))Then
	  Deallocate(eamDataOut)
	End If
    Allocate(eamKeyOut(1:size(eamKeyReducedArray,1),1:5))
    Allocate(eamDataOut(1:totalSplineDataPoints,1:3))
!
! Spline nodes with matching zero, first and second order derivatives (6 points, 5th order poly)
!
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyReducedArray,1)
!Get number of points to expand to
	  If(eamKeyReducedArray(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV	!default 1001
	  End If	
	  If(eamKeyReducedArray(potKey,3).eq.2.or.eamKeyReducedArray(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP	!default 1001
	  End If		
	  If(eamKeyReducedArray(potKey,3).eq.3.or.eamKeyReducedArray(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF	!default 1001
	  End If
!make expanded set of data points
	  potStart = eamKeyReducedArray(potKey,4)
      potLength = eamKeyReducedArray(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyOut(potKey,1) = eamKeyReducedArray(potKey,1) !Atom i
      eamKeyOut(potKey,2) = eamKeyReducedArray(potKey,2) !Atom j
	  eamKeyOut(potKey,3) = eamKeyReducedArray(potKey,3) !Pot type
	  eamKeyOut(potKey,4) = dataPointCounter
	  eamKeyOut(potKey,5) = numberOfPoints
!loop through reduced data points
      splineXY = Spline(eamDataReducedArray,numberOfPoints,potStart,potLength)
	  Do i=1,numberOfPoints
!store data
        eamDataOut(dataPointCounter,1) = splineXY(i,1)
        eamDataOut(dataPointCounter,2) = splineXY(i,2)
        eamDataOut(dataPointCounter,3) = splineXY(i,3)
!increment counter
	    dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
  End Subroutine splinePotentialOld
!-----------------------------
!Make a trial set of points (expanded from the reduced set, using interp or 0,1st,2nd order spline)
  Function eamWobbliness(eamKeyArray, eamDataArray, numberOfPoints) result (wobbliness)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k
	Integer(kind=StandardInteger) :: potKey, numberOfPoints, potStart, potLength, potEnd
	Real(kind=DoubleReal) :: x,y,dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: wobbliness
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
!Allocate arrays	
	Allocate(yArray(1:2)) 
	Allocate(wobbliness(1:size(eamKeyArray,1)))
!Loop through potential functions
	Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
!loop through reduced data points
	  wobbliness(potKey) = 0.0D0
	  Do i=1,numberOfPoints
	    x = eamDataArray(potStart,1)+&
		    1.0D0*(i-1)*((eamDataArray(potEnd,1)-&
			eamDataArray(potStart,1))/(numberOfPoints-1))
		!yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength,"N")
		yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
		wobbliness(potKey) = 1.0D0*wobbliness(potKey)+abs(1.0D0*dy)
	  End Do	
    End Do    
  End Function eamWobbliness
!-----------------------------
  Function eamTotalWobbliness(eamKeyArray, eamDataArray, numberOfPoints) result (wobbliness)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, numberOfPoints
	Real(kind=DoubleReal) :: wobbliness
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: wobblinessArray
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	wobblinessArray = eamWobbliness(eamKeyArray, eamDataArray, numberOfPoints)
 !Loop
    wobbliness = 0.0D0
    Do i=1,size(wobblinessArray,1)
	  wobbliness = wobbliness + wobblinessArray(i)
	End Do    
  End Function eamTotalWobbliness
!-----------------------------
!Make a trial set of points (expanded from the reduced set, using interp or 0,1st,2nd order spline)
  Function eamCurvature(eamKeyArray, eamDataArray, numberOfPointsIn) result (curvature)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, posOffset,totalDataPoints
	Integer(kind=StandardInteger) :: potKey, numberOfPoints, potStart, potLength, potEnd
	Real(kind=DoubleReal) :: x,y,xA,xB,h,yA,yB,dy,dyA,ddy,totalCurvature
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: curvature
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPoints
	Integer(kind=StandardInteger), optional :: numberOfPointsIn
!optional variables	
	If(Present(numberOfPointsIn))Then
	  numberOfPoints = numberOfPointsIn
	Else
	  numberOfPoints = 1000
	End If
	totalDataPoints = numberOfPoints * size(eamKeyArray,1)
!Allocate arrays	
	Allocate(yArray(1:2)) 
	Allocate(dataPoints(1:totalDataPoints,1:2))
	Allocate(curvature(1:totalDataPoints))
!Loop through potential functions
    k = 0
	Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
!loop through reduced data points
      h = ((eamDataArray(potEnd,1)-eamDataArray(potStart,1))/(numberOfPoints-1))
	  Do i=1,numberOfPoints
	    k = k + 1
	    x = eamDataArray(potStart,1)+1.0D0*(i-1)*h
		!yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength,"N")
		yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
		xA = x+0.2D0*h
		!yArray = PointInterpolationArr(eamDataArray,xA,5,potStart,potLength,"N")
		yArray = PointInterpolationArr(eamDataArray,xA,5,potStart,potLength)
		yA = yArray(1)
		dyA = yArray(2)
		ddy = (yA-y)/(0.2D0*h)		
		If(mpiProcessID.eq.0)Then
		  print *,x,y,dy,ddy
		End If
		dataPoints(k,1) = x
		dataPoints(k,2) = y
	  End Do	
    End Do   
!Loop through potential functions
    k = 0
	Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
!loop through reduced data points
      h = ((eamDataArray(potEnd,1)-eamDataArray(potStart,1))/(numberOfPoints-1))
	  Do i=1,numberOfPoints
		k = k + 1
		posOffset = 0
	    If(k.eq.potStart)Then
		  posOffset = 1
		End If		
	    If(k.eq.potEnd)Then
		  posOffset = -1
		End If
		x = dataPoints(k,1)
		yB = dataPoints(k-1+posOffset,2)
		y = dataPoints(k+posOffset,2)
		yA = dataPoints(k+1+posOffset,2)
		curvature(k) = (yA-2*y+yB)/(h**2)
	  End Do	
    End Do 	
  End Function eamCurvature
!-----------------------------
  Function eamCurveLength(eamKeyArray, eamDataArray, numberOfPointsIn) result (curveLength) 
 !force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, posOffset,totalDataPoints
	Integer(kind=StandardInteger) :: potKey, numberOfPoints, potStart, potLength, potEnd
	Real(kind=DoubleReal) :: x,y,xA,xB,h,yA,yB,dy,dyA,ddy,totalCurvature
	Real(kind=DoubleReal) :: xLast, yLast
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: curveLength
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray 
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: dataPoints
	Integer(kind=StandardInteger), optional :: numberOfPointsIn
!optional variables	
	If(Present(numberOfPointsIn))Then
	  numberOfPoints = numberOfPointsIn
	Else
	  numberOfPoints = 1000
	End If
	totalDataPoints = numberOfPoints * size(eamKeyArray,1)
!Allocate arrays	
	Allocate(yArray(1:2)) 
	Allocate(curveLength(1:(size(eamKeyArray,1)+1)))
!Loop through potential functions
    k = 0
	curveLength(size(eamKeyArray,1)+1) = 0.0D0
	Do potKey=1,size(eamKeyArray,1)
!make expanded set of data points
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
!loop through reduced data points
      h = ((eamDataArray(potEnd,1)-eamDataArray(potStart,1))/(numberOfPoints-1))
	  xLast = eamDataArray(potStart,1)
	  yLast = eamDataArray(potStart,2)
	  curveLength(potKey) = 0.0D0
	  Do i=2,numberOfPoints
	    x = eamDataArray(potStart,1)+1.0D0*(i-1)*h
		!yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength,"N")
		yArray = PointInterpolationArr(eamDataArray,x,5,potStart,potLength)
		y = yArray(1)
		curveLength(potKey)=curveLength(potKey)+((xLast-x)**2+(yLast-y)**2)**0.5
		xLast = x
		yLast = y
	  End Do	
	  curveLength(size(eamKeyArray,1)+1)=curveLength(size(eamKeyArray,1)+1)+curveLength(potKey) 
    End Do 	
  End Function eamCurveLength 
  
  
   
  
  
  
!Reduced set of points from input eam
  Subroutine makeReducedEAMSetOld()  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalReducedDataPoints = 0
    Do potKey=1,size(eamKey,1)
!Get number of points to reduce to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If	 
      totalReducedDataPoints = totalReducedDataPoints + numberOfPoints
    End Do	
!Allocate arrays
    If(Allocated(eamKeyReduced))Then
	  Deallocate(eamKeyReduced)
	End If
    If(Allocated(eamDataReduced))Then
	  Deallocate(eamDataReduced)
	End If
    Allocate(eamKeyReduced(1:size(eamKey,1),1:5))
    Allocate(eamDataReduced(1:totalReducedDataPoints,1:3))
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKey,1)
!Get number of points to reduce to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If
!make reduces set of data points
      !eamKeyReduced, eamDataReduced
	  potStart = eamKey(potKey,4)
      potLength = eamKey(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyReduced(potKey,1) = eamKey(potKey,1) !Atom i
      eamKeyReduced(potKey,2) = eamKey(potKey,2) !Atom j
	  eamKeyReduced(potKey,3) = eamKey(potKey,3) !Pot type
	  eamKeyReduced(potKey,4) = dataPointCounter
	  eamKeyReduced(potKey,5) = numberOfPoints
!loop through reduced data points
	  Do i=1,numberOfPoints
	    x = eamData(potStart,1)+&
		    1.0D0*(i-1)*((eamData(potEnd,1)-eamData(potStart,1))/(numberOfPoints-1))
		!yArray = PointInterpolationArr(eamData,x,5,potStart,potLength,"N")
		yArray = PointInterpolationArr(eamData,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
!store data
        eamDataReduced(dataPointCounter,1) = x
        eamDataReduced(dataPointCounter,2) = y
        eamDataReduced(dataPointCounter,3) = dy
!increment counter
		dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
  End Subroutine makeReducedEAMSetOld
  
  
  
!Reduced set of points from trial eam
  Subroutine makeReducedEAMSetTrial()  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalReducedDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Character(len=5)  :: eamTypeText
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalReducedDataPoints = 0
    Do potKey=1,size(eamKeyTrial,1)
!Get number of points to reduce to
	  If(eamKeyTrial(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKeyTrial(potKey,3).eq.2.or.eamKeyTrial(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKeyTrial(potKey,3).eq.3.or.eamKeyTrial(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If	 
      totalReducedDataPoints = totalReducedDataPoints + numberOfPoints
    End Do	
!Allocate arrays
    If(Allocated(eamKeyReduced))Then
	  Deallocate(eamKeyReduced)
	End If
    If(Allocated(eamDataReduced))Then
	  Deallocate(eamDataReduced)
	End If
    Allocate(eamKeyReduced(1:size(eamKeyTrial,1),1:5))
    Allocate(eamDataReduced(1:totalReducedDataPoints,1:3))
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyTrial,1)
!Get number of points to reduce to
	  If(eamKeyTrial(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamReducedPointsV
	  End If	
	  If(eamKeyTrial(potKey,3).eq.2.or.eamKeyTrial(potKey,3).eq.4)Then
	    numberOfPoints = eamReducedPointsP
	  End If		
	  If(eamKeyTrial(potKey,3).eq.3.or.eamKeyTrial(potKey,3).eq.5)Then
	    numberOfPoints = eamReducedPointsF
	  End If
!make reduces set of data points
	  potStart = eamKeyTrial(potKey,4)
      potLength = eamKeyTrial(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyReduced(potKey,1) = eamKeyTrial(potKey,1) !Atom i
      eamKeyReduced(potKey,2) = eamKeyTrial(potKey,2) !Atom j
	  eamKeyReduced(potKey,3) = eamKeyTrial(potKey,3) !Pot type
	  eamKeyReduced(potKey,4) = dataPointCounter
	  eamKeyReduced(potKey,5) = numberOfPoints
!loop through reduced data points
	  Do i=1,numberOfPoints
	    x = eamDataTrial(potStart,1)+&
		    1.0D0*(i-1)*((eamDataTrial(potEnd,1)-eamDataTrial(potStart,1))/(numberOfPoints-1))
		yArray = PointInterpolationArr(eamDataTrial,x,5,potStart,potLength)
		y = yArray(1)
		dy = yArray(2)
!store data
        eamDataReduced(dataPointCounter,1) = x
        eamDataReduced(dataPointCounter,2) = y
        eamDataReduced(dataPointCounter,3) = dy
!increment counter
		dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
!Store reduced potential to file
	if(mpiProcessID.eq.0)then
	  outputFile = trim(currentWorkingDirectory)//"/"//"outputreduced.pot"
	  open(unit=22,file=trim(outputFile))
	  Do potKey=1,size(eamKeyReduced,1)
!Get number of points to reduce to
	    If(eamKeyReduced(potKey,3).eq.1)Then     !Pair
	      write(22,"(A5,A2,A1,A2,A1)") "PAIR ",elements(eamKeyReduced(potKey,1))," ",&
		  elements(eamKeyReduced(potKey,2))," "
		Else  
		  If(eamKeyReduced(potKey,3).eq.2.and.eamType.eq.1)Then
		    eamTypeText = "DENS "
		  End If
		  If(eamKeyReduced(potKey,3).eq.2.and.eamType.eq.2)Then
		    eamTypeText = "SDEN "
		  End If
		  If(eamKeyReduced(potKey,3).eq.3.and.eamType.eq.1)Then
		    eamTypeText = "EMBE "
		  End If
		  If(eamKeyReduced(potKey,3).eq.3.and.eamType.eq.2)Then
		    eamTypeText = "SEMB "
		  End If
		  If(eamKeyReduced(potKey,3).eq.4)Then
		    eamTypeText = "DDEN "
		  End If
		  If(eamKeyReduced(potKey,3).eq.5)Then
		    eamTypeText = "DEMB "
		  End If
		  write(22,"(A5,A2)") eamTypeText,elements(eamKeyReduced(potKey,1))
	    End If	
!make reduces set of data points
	    potStart = eamKeyReduced(potKey,4)
        potLength = eamKeyReduced(potKey,5)
	    potEnd = potStart + potLength - 1
!loop over data points
        k = 0
	    Do i=potStart,potEnd
		  k = k + 1
	      x = eamDataReduced(i,1)
		  y = eamDataReduced(i,2)
		  dy = eamDataReduced(i,3)
!write to file
          write(22,"(E24.16E3,A2,E24.16E3,A2,E24.16E3,A4,I8,I8,I8)") x,"  ",&
		  y,"  ",dy,"    ",potKey,k,i
	    End Do	
      End Do
!Close file
	  close(22)
	endif
  End Subroutine makeReducedEAMSetTrial
  
  
  
!Make a trial set of points (expanded from the reduced set, using interp or 0,1st,2nd order spline)
  Subroutine makeTrialEAMSet(fitTypeIn,eamKeyArrayIn,eamDataArrayIn)  	
!force declaration of all variables
	Implicit None
!declare private variables
	Integer(kind=StandardInteger) :: ios, i, j, k, potKey
	Integer(kind=StandardInteger) :: numberOfPoints, totalTrialDataPoints
	Integer(kind=StandardInteger) :: potStart, potLength, potEnd, dataPointCounter
	Real(kind=DoubleReal) :: x, y, dy
	Real(kind=DoubleReal), Dimension( : ), Allocatable :: yArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: splineXY
	Character(len=5)  :: eamTypeText
!Declare optional variables
    Integer(kind=StandardInteger), Optional :: fitTypeIn
    Integer(kind=StandardInteger) :: fitType
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable, Optional :: eamKeyArrayIn 
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: eamKeyArray
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable, Optional :: eamDataArrayIn
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: eamDataArray
!Set optional variables
    If(Present(fitTypeIn))Then
	  fitType = fitTypeIn
	Else  
	  fitType = 1
	End If
    If(Present(eamKeyArrayIn))Then
	  eamKeyArray = eamKeyArrayIn
	Else  
	  eamKeyArray = eamKeyReduced
	End If
    If(Present(eamDataArrayIn))Then
	  eamDataArray = eamDataArrayIn
	Else  
	  eamDataArray = eamDataReduced
	End If
!Allocate arrays	
	Allocate(yArray(1:2))
!Sum number of reduced data points
    totalTrialDataPoints = 0
    Do potKey=1,size(eamKeyArray,1)
!Get number of points to reduce to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF
	  End If	 
      totalTrialDataPoints = totalTrialDataPoints + numberOfPoints
    End Do
!Allocate arrays
    If(Allocated(eamKeyTrial))Then
	  Deallocate(eamKeyTrial)
	End If
    If(Allocated(eamDataTrial))Then
	  Deallocate(eamDataTrial)
	End If
    Allocate(eamKeyTrial(1:size(eamKeyArray,1),1:5))
    Allocate(eamDataTrial(1:totalTrialDataPoints,1:3))
!
! Spline nodes with matching zero, first and second order derivatives (6 points, 5th order poly)
!
!Loop through potential functions
    dataPointCounter = 1
	Do potKey=1,size(eamKeyTrial,1)
!Get number of points to expand to
	  If(eamKey(potKey,3).eq.1)Then     !Pair
	    numberOfPoints = eamTrialPointsV
	  End If	
	  If(eamKey(potKey,3).eq.2.or.eamKey(potKey,3).eq.4)Then
	    numberOfPoints = eamTrialPointsP
	  End If		
	  If(eamKey(potKey,3).eq.3.or.eamKey(potKey,3).eq.5)Then
	    numberOfPoints = eamTrialPointsF
	  End If
!make expanded set of data points
	  potStart = eamKeyArray(potKey,4)
      potLength = eamKeyArray(potKey,5)
	  potEnd = potStart + potLength - 1
!store key details
      eamKeyTrial(potKey,1) = eamKeyArray(potKey,1) !Atom i
      eamKeyTrial(potKey,2) = eamKeyArray(potKey,2) !Atom j
	  eamKeyTrial(potKey,3) = eamKeyArray(potKey,3) !Pot type
	  eamKeyTrial(potKey,4) = dataPointCounter
	  eamKeyTrial(potKey,5) = numberOfPoints
!loop through reduced data points
      splineXY = Spline(eamDataArray,numberOfPoints,potStart,potLength)
	  Do i=1,numberOfPoints
!store data
        eamDataTrial(dataPointCounter,1) = splineXY(i,1)
        eamDataTrial(dataPointCounter,2) = splineXY(i,2)
        eamDataTrial(dataPointCounter,3) = splineXY(i,3)
!increment counter
	    dataPointCounter = dataPointCounter + 1
	  End Do	
    End Do
  End Subroutine makeTrialEAMSet

  
  

End Module prep