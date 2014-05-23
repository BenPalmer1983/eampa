Module prep

! Setup Modules
  Use kinds
  Use constants
  Use mpif
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
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtoms					!count of total atoms for each config
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtomsUnitCell			!
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: configAtomsStart	
  Integer(kind=StandardInteger) :: configAtomsTotal											!Total atoms in all configs
  Integer(kind=StandardInteger) :: configAtomsMax											!Max atoms in 1 config
  Integer(kind=StandardInteger) :: configAtomsAdd											!Max rounded up to nearest mult of 10
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: configAtomsMap			!Start-End points for atoms by config  
  Integer(kind=StandardInteger), Dimension( : ), Allocatable :: globalCounter				!Array of global counts	
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: configLatticeParameters
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: globalTimer
  
!Neighbour list arrays
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKey		!key for atom neighbour list
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListI          !neighbour list integer data
  Real(kind=SingleReal), Dimension( : , : ), Allocatable :: neighbourListCoords
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListR
  
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKeyStore
  Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListIStore
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: neighbourListCoordsStore
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListRStore
  
!calculated values
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationVolume  
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationEnergy
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationBM
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptVolume
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationOptEnergy
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationForceX
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationForceY
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationForceZ
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationStress
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationUnitVector
  
!Reference energies, forces, stresses etc  
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationRefEnergy
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationRefForceX
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationRefForceY
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationRefForceZ
  Real(kind=DoubleReal), Dimension( : ), Allocatable :: configurationRefBM
  
!RSS configuration array    
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configurationRSS
    
!RSS configuration opt weights    
  Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configurationOptWeights
  
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
  Public :: configurationVolume
  Public :: configurationEnergy	
  Public :: configurationBM	
  Public :: configurationOptVolume	
  Public :: configurationOptEnergy	
  Public :: configAtomsTotal
  Public :: configAtomsMax
  Public :: configAtomsAdd
  Public :: configurationForceX
  Public :: configurationForceY
  Public :: configurationForceZ
  Public :: configurationForceMPIKey
  Public :: configurationRefEnergy
  Public :: configurationRefForceX
  Public :: configurationRefForceY
  Public :: configurationRefForceZ
  Public :: configurationRefBM
  Public :: configurationStress
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
  Public :: storeEAMToFile
  Public :: printEAM
  Public :: storeToReduced
  Public :: synchMpiProcesses
  Public :: storeNeighboutList
  Public :: loadNeighboutList
  Public :: applyDistortionVector
  Public :: calcConfigurationVolume
  Public :: clearEAM,clearOptEAM,clearReducedEAM,clearTrialEAM
  Public :: splinePotential
  
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
!Prep data 
    Call prepDataArrays()
	Call prepDataVariables()
!Make neighbour list
	Call makeNeighbourList()
	!Call applyUnitVector(0)
	Call orderEamPotentials()
	If(mpiProcessID.eq.0)Then
	  Call storeEAM()
    End If
	!Call makeReducedEAMSet()
	!Call makeTrialEAMSet()
	Call synchMpiProcesses()
	Call prepDataStore()
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
	Integer(kind=StandardInteger) :: energyCounter, bulkModulusCounter
!-
!This sets up important data arrays needed here and elsewhere in the program
!-  
!Allocate arrays	  
    Allocate(configurationRefEnergy(1:configCount))
    Allocate(configurationEnergy(1:configCount))		!Store calculated configuration energies	  
	Allocate(configAtoms(1:configCount)) 
    Allocate(configAtomsStart(1:configCount))	  
    Allocate(configAtomsMap(1:configCount,1:50))		!information stored about each config
	Allocate(configurationStress(1:(9*configCount)))
	Allocate(configurationVolume(1:configCount))
	Allocate(configurationBM(1:configCount))
	Allocate(configurationRefBM(1:configCount))
	Allocate(configurationUnitVector(1:(9*configCount)))
	Allocate(configLatticeParameters(1:configCount,1:3))
	Allocate(configurationRadiusCutoff(1:configCount))
	Allocate(configurationRSS(1:configCount,1:10))
	Allocate(configurationOptWeights(1:configCount,1:4))
	
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
!counters
    energyCounter = 0
	bulkModulusCounter = 0
!Loop through configurations
    Do i=1,configCount
	  configurationVolume(i) = 0.0D0
      configurationEnergy(i) = 0.0D0	
	  configurationRefEnergy(i) = configHeaderR(i,12)
	  configurationBM(i) = -2.1D20
	  configurationRefBM(i) = configHeaderR(i,15)		!No reference data if less than -2.0D20
	  !print *,i,configurationRefBM(i)
	  configAtoms(i) = configHeaderI(i,10) * configHeaderI(i,11) *&
	  configHeaderI(i,12) * configHeaderI(i,headerWidth)
	  configAtomsTotal = configAtomsTotal + configAtoms(i)
	  Do j=1,10
	    configurationRSS(i,j) = 0.0D0
	  End Do
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
	  If(configHeaderR(i,14).gt.-2.0E20)Then !If ref stress set
	    configAtomsMap(i,6) = 1 
	  Else
	    configAtomsMap(i,6) = -1 
	  End If  
	  If(configHeaderR(i,12).gt.-2.0E20)Then !If ref energy set
	    configAtomsMap(i,5) = 1 
	  Else
	    configAtomsMap(i,5) = -1 
	  End If	  
!MPI details 21-30
!Increment counters for MPI
	  energyCounter = energyCounter + 1
	  configAtomsMap(i,21) = mod(energyCounter-1,mpiProcessCount)		!Assigned MPI process
	  If(configurationRefBM(i).gt.-2.0D20)Then
	    bulkModulusCounter = bulkModulusCounter + 1
		configAtomsMap(i,22) = mod(bulkModulusCounter-1,mpiProcessCount)
	  Else
	    configAtomsMap(i,22) = -1
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
!	  
!End loop through configurations
	End Do	
!If forces are to be calculated
	If(calcForcesOnOff.eq.1)Then	  
	  !Keeps track of which MPIprocess calculates which set of forces
      Allocate(configurationForceMPIKey(1:configCount,1:3))	 
      Do i=1,configCount
        configurationForceMPIKey(i,1) = -1  !Process
        configurationForceMPIKey(i,2) = -1  !Start
        configurationForceMPIKey(i,3) = -1  !End
	  End Do
      Allocate(configurationForceX(1:configAtomsTotal))
      Allocate(configurationForceY(1:configAtomsTotal))
      Allocate(configurationForceZ(1:configAtomsTotal))
	  !zero out forces
      Do i=1,configAtomsTotal
        configurationForceX(i) = 0.0D0
        configurationForceY(i) = 0.0D0
        configurationForceZ(i) = 0.0D0	
	  End Do
    End If		
!If stress is to be calculated
	If(calcStressOnOff.eq.1)Then
	  Do i=1,(9*configCount)
	    configurationStress(i) = 0.0D0
	  End Do
	End If	
!Atom type key array
	Allocate(atomTypeKey(1:configAtomsTotal))
!Global counter - counting energy calculations, optimisation steps etc
    Allocate(globalCounter(1:10)) 	!1 energy calc counter, 2
	Do i=1,10
	  globalCounter(i) = 0
	End Do
!Global timer 
    Allocate(globalTimer(1:22)) 	
	!1 prep time, 2 prep time real
	!3 energy calc time, 4 energy calc time real
	!5 bm calc, 6 bm calc real
	!
    Do i=1,20
	  globalTimer(i) = 0.0D0
	End Do
  
  
  
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
	Integer(kind=StandardInteger) :: atomA, atomB
	Integer(kind=StandardInteger) :: neighbourListLength, tempNeighbourListLength
	Integer(kind=StandardInteger) :: neighbourListCount, configStart, configLength
	Integer(kind=StandardInteger) :: coordsStart, coordsLength
	Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
	Integer(kind=StandardInteger) :: atomCounter
	Real(kind=SingleReal) :: rCutoff, rCutoffNeg, rCutoffSq
	Real(kind=SingleReal) :: xcoord, ycoord, zcoord
	Real(kind=SingleReal) :: xshift, yshift, zshift
	Real(kind=SingleReal) :: xlat, ylat, zlat, alat
	Real(kind=SingleReal) :: xA, xB, yA, yB, zA, zB
	Real(kind=SingleReal) :: rdSq, xdSq, ydSq, zdSq
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListKeyTemp
	Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: neighbourListITemp
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: neighbourListRTemp
    Real(kind=SingleReal), Dimension( : , : ), Allocatable :: neighbourListCoordsTemp
	Integer(kind=StandardInteger), Dimension( : ), Allocatable :: coordsITemp
	Real(kind=SingleReal), Dimension( : , : ), Allocatable :: coordsRTemp	
!Unit vectors
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: configUnitVector
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: workingUnitVector
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
!open coords file	
    if(saveFileCoords.eq."Y")then
	  outputFile = trim(currentWorkingDirectory)//"/"//"coords.dat"
	  open(unit=10,file=trim(outputFile))	
	endif		
!write to output file
    If(mpiProcessID.eq.0)Then
      write(999,"(F8.4,A2,A37)") ProgramTime(),"  ",&
	  "Start building config neighbour lists"
	End If
	!configHeaderR
	!1 LP, 2 RC
	!configHeaderI
	!1 X1, 2 X2, 3 X3, 4 Y1, 5 Y2, 6 Y3, 7 Z1, 8 Z2, 9 Z3, 10 CC1, 11 CC2, 12 CC3, 13 Start, 14 Length    
	!configCoordsI
	!configCoordsR
	
!Unit vector used to transpose coords    
    Allocate(configUnitVector(1:3,1:3))	

!Allocate temp neighbour list array - assume 100 neighbours per atom
	tempNeighbourListLength = 0
	Do i=1,configCount
      tempNeighbourListLength = tempNeighbourListLength + &
	  configHeaderI(i,10) * configHeaderI(i,11)* configHeaderI(i,12) * &
	  configHeaderI(i,headerWidth) * 100
	Enddo
	Allocate(neighbourListKeyTemp(1:size(configHeaderI)/headerWidth,1:2)) 
	Allocate(neighbourListITemp(1:tempNeighbourListLength,1:4))
	Allocate(neighbourListRTemp(1:tempNeighbourListLength)) 
	Allocate(neighbourListCoordsTemp(1:tempNeighbourListLength,1:6)) 
	Allocate(configAtomsUnitCell(1:configCount)) 
	
!Allocate config reference data arrays
    Allocate(configurationRefForceX(1:configAtomsTotal))
    Allocate(configurationRefForceY(1:configAtomsTotal))
    Allocate(configurationRefForceZ(1:configAtomsTotal))
	
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
	      !print *,((i-1)*9+n),configurationUnitVector((i-1)*9+n)
		End Do
	  End Do
	  configurationVolume(i) = calcConfigurationVolume(i)
	  !print *,"Vol",i,configurationVolume(i) 
!expand unit cell	    
	  !atoms = xCopy * yCopy * zCopy * configHeaderI(i,headerWidth)
	  !configAtoms(i) = atoms
	  atoms = configAtoms(i)
	  configAtomsUnitCell(i) = configHeaderI(i,headerWidth)
!Allocate temp array
	  if(Allocated(coordsITemp))then
	    Deallocate(coordsITemp)
	  endif
	  if(Allocated(coordsRTemp))then
	    Deallocate(coordsRTemp)
	  endif
	  Allocate(coordsITemp(1:atoms))
	  Allocate(coordsRTemp(1:atoms,1:3))
	  m = 0
!loop over copies to make full crystal	  
	  coordsStart = configHeaderI(i,headerWidth-1)
	  coordsLength = configHeaderI(i,headerWidth)	  
	  do x=1,xCopy
	    do y=1,yCopy
	      do z=1,zCopy
		    do n=coordsStart,(coordsStart+coordsLength-1)  
			  m = m + 1
			  atomCounter = atomCounter + 1
!get co-ords with lattice parameter and copy applied
			  coordsITemp(m) = configCoordsI(n)
			  coordsRTemp(m,1) = alat * (x + configCoordsR(n,1) - 1)
			  coordsRTemp(m,2) = alat * (y + configCoordsR(n,2) - 1)
			  coordsRTemp(m,3) = alat * (z + configCoordsR(n,3) - 1)
!Apply configuration unit vector to these co-ordinates
!x-coord
              coordsRTemp(m,1) = & 		
			  coordsRTemp(m,1) * configHeaderR(i,3) + &
              coordsRTemp(m,2) * configHeaderR(i,4) + &
              coordsRTemp(m,3) * configHeaderR(i,5)
!y-coord
              coordsRTemp(m,2) = & 		
			  coordsRTemp(m,1) * configHeaderR(i,6) + &
              coordsRTemp(m,2) * configHeaderR(i,7) + &
              coordsRTemp(m,3) * configHeaderR(i,8)
!z-coord
              coordsRTemp(m,3) = & 		
			  coordsRTemp(m,1) * configHeaderR(i,9) + &
              coordsRTemp(m,2) * configHeaderR(i,10) + &
              coordsRTemp(m,3) * configHeaderR(i,11)
!save coords to file if required
              if(saveFileCoords.eq."Y")then
	            write(10,"(I8,I8,I8,F16.8,F16.8,F16.8)") i,m,coordsITemp(m),&
			    coordsRTemp(m,1),coordsRTemp(m,2),coordsRTemp(m,3)
			  endif
!store atom-force data
              configurationRefForceX(atomCounter) = configForcesR(n,1)
              configurationRefForceY(atomCounter) = configForcesR(n,2)
              configurationRefForceZ(atomCounter) = configForcesR(n,3)
			  atomTypeKey(atomCounter) = configCoordsI(n)
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
	  !configLatticeParameters(i,1) = xlat
	  !configLatticeParameters(i,2) = ylat
	  !configLatticeParameters(i,3) = zlat
!loop through atom pairs    Atom A in inner cell, Atom B in 3x3x3 supercell      
	  configLength = 0
      do atomA=1,atoms
!loop through Atom B 3x3x3
	    do l=-1,1
	      do m=-1,1
	        do n=-1,1	              			
			  xshift = xlat * l
			  yshift = ylat * m
			  zshift = zlat * n
	          do atomB=1,atoms 			  
!don't self count atom
			    if(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)then
			    else
!check range one co-ord at a time then distance squared
                  xA = 1.0D0*coordsRTemp(atomA,1)
                  xB = 1.0D0*(xshift + coordsRTemp(atomB,1))
                  yA = 1.0D0*coordsRTemp(atomA,2)
                  yB = 1.0D0*(yshift + coordsRTemp(atomB,2))
                  zA = 1.0D0*coordsRTemp(atomA,3)
                  zB = 1.0D0*(zshift + coordsRTemp(atomB,3))
                  xdSq = (xA-xB)**2
				  if(xdSq.le.rCutoffSq)then
				    ydSq = (yA-yB)**2
				    if(ydSq.le.rCutoffSq)then
				      zdSq = (zA-zB)**2
				      if(zdSq.le.rCutoffSq)then
					    rdSq = xdSq + ydSq + zdSq						
					    if(rdSq.le.rCutoffSq)then
					      neighbourListCount = neighbourListCount + 1
						  configLength = configLength + 1
!atom type data
						  neighbourListITemp(neighbourListCount,1) = coordsITemp(atomA) !Atom A Type
						  neighbourListITemp(neighbourListCount,2) = coordsITemp(atomB)	!Atom B Type
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
					    endif
					  endif
					endif
			      endif
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
	Allocate(neighbourListKey(1:size(configHeaderI)/headerWidth,1:2))
	Allocate(neighbourListI(1:neighbourListCount,1:4))
	Allocate(neighbourListR(1:neighbourListCount)) 
	Allocate(neighbourListCoords(1:neighbourListCount,1:12)) !Ax Ay Az, Bx By Bz, R_ABi R_ABj R_ABk 
!Transfer key data 
    !do i=1,size(configHeaderI)/headerWidth
    do i=1,configCount
	  neighbourListKey(i,1) = neighbourListKeyTemp(i,1)
	  neighbourListKey(i,2) = neighbourListKeyTemp(i,2)
	enddo
!Transfer neighbour list data	
    Do i=1,neighbourListCount
	  neighbourListI(i,1) = neighbourListITemp(i,1)  !Atom A Type
	  neighbourListI(i,2) = neighbourListITemp(i,2)  !Atom B Type
	  neighbourListI(i,3) = neighbourListITemp(i,3)  !Atom A ID
	  neighbourListI(i,4) = neighbourListITemp(i,4)  !Atom B ID
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
	End Do
!Deallocate arrays
    Deallocate(neighbourListITemp)
    Deallocate(neighbourListRTemp)
    Deallocate(neighbourListCoordsTemp)

!close coords file
    if(saveFileCoords.eq."Y")then
      close(10)	 	 
    endif	
!output to file	
    If(mpiProcessID.eq.0)Then
      write(999,"(A6,A31,I8)") "      ","Total configurations prepared: ",configCount
	  write(999,"(F8.4,A2,A53,I8)") ProgramTime(),"  ",&
	  "Finish building config neighbour lists, total pairs: ",neighbourListCount
	End If
!close output file
    close(999)	
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
! ***calcConfigurationVolume*** 
!------------------------------------------------------------------------!

  Function calcConfigurationVolume(configurationID, distortionArrayIn) result (volume)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k
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
	multiplierArray(1,1) = configurationUnitVector((i-1)*9+1)
	multiplierArray(1,2) = configurationUnitVector((i-1)*9+2)
	multiplierArray(1,3) = configurationUnitVector((i-1)*9+3)
	multiplierArray(2,1) = configurationUnitVector((i-1)*9+4)
	multiplierArray(2,2) = configurationUnitVector((i-1)*9+5)
	multiplierArray(2,3) = configurationUnitVector((i-1)*9+6)
	multiplierArray(3,1) = configurationUnitVector((i-1)*9+7)
	multiplierArray(3,2) = configurationUnitVector((i-1)*9+8)
	multiplierArray(3,3) = configurationUnitVector((i-1)*9+9)
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
   
  Function calcConfigurationVolumeB(configurationID, distortionArrayIn) result (volume)
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: configurationID
	Integer(kind=StandardInteger) :: i, j, k
    Real(kind=SingleReal) :: xA, yA, zA
    Real(kind=SingleReal) :: xB, yB, zB
    Real(kind=SingleReal) :: xC, yC, zC
	Real(kind=SingleReal) :: crossProductI,crossProductJ,crossProductK
	Real(kind=DoubleReal) :: volume
	Real(kind=DoubleReal), Optional, Dimension( : , : ), Allocatable :: distortionArrayIn
	Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: distortionArray, multiplierArray
!Set distortion array
    Allocate(multiplierArray(1:3,1:3))
	multiplierArray(1,1) = configurationUnitVector((i-1)*9+1)
	multiplierArray(1,2) = configurationUnitVector((i-1)*9+2)
	multiplierArray(1,3) = configurationUnitVector((i-1)*9+3)
	multiplierArray(2,1) = configurationUnitVector((i-1)*9+4)
	multiplierArray(2,2) = configurationUnitVector((i-1)*9+5)
	multiplierArray(2,3) = configurationUnitVector((i-1)*9+6)
	multiplierArray(3,1) = configurationUnitVector((i-1)*9+7)
	multiplierArray(3,2) = configurationUnitVector((i-1)*9+8)
	multiplierArray(3,3) = configurationUnitVector((i-1)*9+9)
    If(present(distortionArrayIn))Then
	
	Else
	
	End If
!loop through configurations if 0 and calc volume for all configurations
    If(configurationID.eq.0)Then
	  Do i=1,configCount
	    xA = configLatticeParameters(i,1) * configurationUnitVector((i-1)*9+1)
        xB = configLatticeParameters(i,1) * configurationUnitVector((i-1)*9+2)
        xC = configLatticeParameters(i,1) * configurationUnitVector((i-1)*9+3)
        yA = configLatticeParameters(i,2) * configurationUnitVector((i-1)*9+4)
        yB = configLatticeParameters(i,2) * configurationUnitVector((i-1)*9+5)
        yC = configLatticeParameters(i,2) * configurationUnitVector((i-1)*9+6)
        zA = configLatticeParameters(i,3) * configurationUnitVector((i-1)*9+7)
        zB = configLatticeParameters(i,3) * configurationUnitVector((i-1)*9+8)
        zC = configLatticeParameters(i,3) * configurationUnitVector((i-1)*9+9)
	    crossProductI = (xB*yC-xC*yB)
	    crossProductJ = (xC*yA-xA*yC)
	    crossProductK = (xA*yB-xB*yA)
	    volume = 1.0D0*zA*crossProductI+zB*crossProductJ+zC*crossProductK	
        configurationVolume(i) = 1.0D0 * volume
	  End Do
	  volume = 0.0D0
	Else
	  i = configurationID
	  xA = configLatticeParameters(i,1) * configurationUnitVector((i-1)*9+1)
      xB = configLatticeParameters(i,1) * configurationUnitVector((i-1)*9+2)
      xC = configLatticeParameters(i,1) * configurationUnitVector((i-1)*9+3)
      yA = configLatticeParameters(i,2) * configurationUnitVector((i-1)*9+4)
      yB = configLatticeParameters(i,2) * configurationUnitVector((i-1)*9+5)
      yC = configLatticeParameters(i,2) * configurationUnitVector((i-1)*9+6)
      zA = configLatticeParameters(i,3) * configurationUnitVector((i-1)*9+7)
      zB = configLatticeParameters(i,3) * configurationUnitVector((i-1)*9+8)
      zC = configLatticeParameters(i,3) * configurationUnitVector((i-1)*9+9)
	  crossProductI = (xB*yC-xC*yB)
	  crossProductJ = (xC*yA-xA*yC)
	  crossProductK = (xA*yB-xB*yA)
	  volume = 1.0D0*zA*crossProductI+zB*crossProductJ+zC*crossProductK	
      configurationVolume(configurationID) = 1.0D0 * volume
	End If  
  End Function calcConfigurationVolumeB	
  
  !Subroutine calcUpdateConfigurationVolume()
  
  
  !End Function calcUpdateConfigurationVolume	
  
  
  Subroutine summariseConfigurations()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i
	Character(2) :: refForce
!If master process
	If(mpiProcessID.eq.0)Then
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
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
! Potential prep subroutines
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
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	
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
  
  Subroutine storeEAMToFile(eamKeyArray, eamDataArray, fileName, numberOfPointsIn)  	
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
	Character(*) :: fileName
!optional variables	
	Integer(kind=StandardInteger), optional :: numberOfPointsIn
	  If(Present(numberOfPointsIn))Then
	    numberOfPoints = numberOfPointsIn
	  Else
        numberOfPoints = 0	  
	  End If
!Store reduced potential to file
	If(mpiProcessID.eq.0)Then
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
	End If
!Wait for all processes to catch up	 
	Call synchMpiProcesses()
  End Subroutine storeEAMToFile  
  
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
  
  
  Subroutine prepDataStore()
!force declaration of all variables
	Implicit None	 
!Internal subroutine variables
    Integer(kind=StandardInteger) :: i, totalAtoms   
  End Subroutine prepDataStore 
   
  
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
		If(mpiProcessID.eq.0)Then
		!print *,k,x,yB,y,yA,curvature(k) 
		End If
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
      !eamKeyReduced, eamDataReduced
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
		!yArray = PointInterpolationArr(eamDataTrial,x,5,potStart,potLength,"N")
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