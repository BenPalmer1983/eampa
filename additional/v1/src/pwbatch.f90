Module pwbatch

!--------------------------------------------------------------!
! PWscf Batch File Subroutines                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 1st May 2014
!----------------------------------------

! Setup Modules
  Use kinds
  Use constants
  Use mpif
  Use general
  Use units
  Use strings		!string functions
  Use maths
  Use initialise
  Use input

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!declare global variables  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: pwbUnitVector, pwbUnitVectorWorking
  Real(kind=DoubleReal) :: pwbLatticeParameter
  Integer(kind=StandardInteger) :: pwbXCopy, pwbYCopy, pwbZCopy	
  Real(kind=DoubleReal) :: pwbXd, pwbYd, pwbZd	
  Character(len=8), Dimension(1:1024)   :: pwbAtomLabelsInput
  Real(kind=DoubleReal), Dimension(1:1024,1:3)   :: pwbAtomCoordsInput  
  Character(len=8), Dimension(1:4096)   :: pwbAtomLabels, pwbAtomLabelsWorking
  Real(kind=DoubleReal), Dimension(1:4096,1:3)   :: pwbAtomCoords, pwbAtomCoordsWorking  	
  Character(len=64), Dimension(1:128,1:2)   :: pwbAtomicSpeciesC  	
  Real(kind=DoubleReal), Dimension(1:128)   :: pwbAtomicSpeciesDP
!PWscf options
  Character(len=64) :: pwbRestartMode, pwbCalculation, pwbOutDir, pwbPseudoDir, pwbPrefix,&
                       pwbDiskIO, pwbOccupations, pwbSmearing, &
                       pwbDiagonalization, pwbMixingMode, pwbIonDynamics, pwbCellDynamics, &
                       pwbKpoints					   
  Character(len=6) :: pwbTprnfor, pwbTstress
  Real(kind=DoubleReal) :: pwbEtotConvThr, pwbForcConvThr, pwbDegauss, pwbMixingBeta, &
                           pwbConvThr, pwbPress, pwbCellFactor
  Integer(kind=StandardInteger) :: pwbNstep, pwbIbrav, pwbNat, pwbNtyp, pwbEcutwfc, pwbEcutrho
  Integer(kind=StandardInteger) :: pwbFixedAtoms
  
!Privacy of functions/subroutines/variables
  Private    
!Variables - Read in from pwbOptConf
  Public :: pwbUnitVector
  Public :: pwbLatticeParameter
  Public :: pwbXCopy, pwbYCopy, pwbZCopy
  Public :: pwbXd, pwbYd, pwbZd	
  Public :: pwbAtomLabels, pwbAtomLabelsInput, pwbAtomLabelsWorking
  Public :: pwbAtomCoords, pwbAtomCoordsInput, pwbAtomCoordsWorking
!Variables - PWscf options 
  Public :: pwbRestartMode, pwbCalculation, pwbOutDir, pwbPseudoDir, pwbPrefix, &
            pwbDiskIO, pwbOccupations, pwbSmearing, &
            pwbDiagonalization, pwbMixingMode, pwbIonDynamics, pwbCellDynamics, &
            pwbKpoints			
  Public :: pwbTprnfor, pwbTstress
  Public :: pwbEtotConvThr, pwbForcConvThr, pwbDegauss, pwbMixingBeta, &
            pwbConvThr, pwbPress, pwbCellFactor
  Public :: pwbNstep, pwbIbrav, pwbNat, pwbNtyp, pwbEcutwfc, pwbEcutrho
  Public :: pwbFixedAtoms
!Public subroutines
  Public :: runPWBatch
  
  
  
Contains
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Subroutines
!-------------------------------------------------------------------------  
! 


  Subroutine runPWBatch() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"PWscf batch files"
	  print *,ProgramTime(),trim(pwbOptConfFile)
	End If
!Read in pwbOptConf files
    Call readPWBatchConf()
!prep config
    Call prepPWBatchConf()
!make pwscf batch files
    If(pwbType.eq."BASIC")Then
      Call makePWBasic()
	End If
	If(pwbType.eq."SLAB")Then
      Call makePWSlab()
	End If
    If(pwbType.eq."FULL")Then
      Call makePWBatchFiles()
	End If
	
!write file
    !If(mpiProcessID.eq.0)Then	
	!  Call writePWscfFile("output/batch/pwscfbatch.in")
	!End If
	
  End Subroutine runPWBatch    
!------------------------------------------------------------------
! Read pwb opt file 
!------------------------------------------------------------------  
  Subroutine readPWBatchConf()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, atomCounter, atomicSpeciesCounter	
	Character(len=128) :: inputFile
	Character(len=255) :: fileRow, fileRowOrig	
	Character(len=64) :: bufferA, bufferB, bufferC, bufferD
!Default values
!--------------------------------------------------
    pwbUnitVector = 0.0D0	
	pwbAtomLabelsInput = "#BLANK##"
	pwbAtomCoordsInput = -2.1D20
	pwbAtomicSpeciesC = "#BLANK##"
	pwbAtomicSpeciesDP = -2.1D20
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
    pwbEtotConvThr = 1.0D-7
	pwbForcConvThr = 1.0D-7
	pwbDegauss = 0.05
	pwbMixingBeta = 0.7
	pwbConvThr = 1.0D-7
	pwbPress = 0.0
	pwbCellFactor = 2.0
!Default pwscf file values, integer
	pwbNstep = 50
	pwbIbrav = 14
	pwbNat = 0
	pwbNtyp = 0
	pwbEcutwfc = 40
	pwbEcutrho = 200
!Fixed atoms in coords
    pwbFixedAtoms = 0	
!--------------------------------------------------		
!re-read user input file for specific pwbatch entries to replace defaults	
!open & read in file	
  	Open(UNIT=1,FILE="temp/tmpinput.in") 
    Do i=1,maxFileRows 
!Read in line
	  Read(1,"(A255)",IOSTAT=ios) fileRow
	  If (ios /= 0) Then
	    EXIT 
	  End If
	  If(StrToUpper(fileRow(1:10)).eq."#PWBOUTDIR")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		pwbOutDir = trim(adjustl(fileRow))
	  End If
	  If(StrToUpper(fileRow(1:13)).eq."#PWBPSEUDODIR")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		pwbPseudoDir = trim(adjustl(fileRow))
	  End If
	  If(StrToUpper(fileRow(1:11)).eq."#PWBECUTWFC")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		Read(fileRow,*) pwbEcutwfc
	  End If
	  If(StrToUpper(fileRow(1:11)).eq."#PWBECUTRHO")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		Read(fileRow,*) pwbEcutrho
	  End If
	  If(StrToUpper(fileRow(1:14)).eq."#PWBMIXINGMODE")then
!read next line
	    Read(1,*,IOSTAT=ios) fileRow
		Read(fileRow,*) pwbMixingMode
	  End If
	  
	End Do
!--------------------------------------------------	
!open output file	
	inputFile = trim(currentWorkingDirectory)//"/"//pwbOptConfFile
	Open(unit=101,file=trim(inputFile))   
	atomCounter = 0
	atomicSpeciesCounter = 0
	Do i=1,maxFileRows 
!Read in line
	  Read(101,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
	  If (ios /= 0) Then
	    EXIT 
	  End If
	  fileRowOrig = fileRow
	  fileRow = StrToUpper(fileRow)
!Lattice parameter
      If(fileRow(1:3).eq."#LP")Then
	    Read(fileRow,*) bufferA, bufferB
		Read(bufferB,*) pwbLatticeParameter
	  End If
!Unit vectors - x
      If(fileRow(1:2).eq."#X")Then
	    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
		Read(bufferB,*) pwbUnitVector(1,1)
		Read(bufferC,*) pwbUnitVector(1,2)
		Read(bufferD,*) pwbUnitVector(1,3)
	  End If
!Unit vectors - y
      If(fileRow(1:2).eq."#Y")Then
	    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
		Read(bufferB,*) pwbUnitVector(2,1)
		Read(bufferC,*) pwbUnitVector(2,2)
		Read(bufferD,*) pwbUnitVector(2,3)
	  End If
!Unit vectors - z
      If(fileRow(1:2).eq."#Z")Then
	    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
		Read(bufferB,*) pwbUnitVector(3,1)
		Read(bufferC,*) pwbUnitVector(3,2)
		Read(bufferD,*) pwbUnitVector(3,3)
	  End If
!Copies
      If(fileRow(1:3).eq."#CC")Then
	    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
		Read(bufferB,*) pwbXCopy
		Read(bufferC,*) pwbYCopy
		Read(bufferD,*) pwbZCopy
	  End If
!Atom Speciec
      If(fileRow(1:3).eq."#AS")Then
	    atomicSpeciesCounter = atomicSpeciesCounter + 1
	    Read(fileRowOrig,*) bufferA, bufferB, bufferC, bufferD
		Read(bufferB,*) pwbAtomicSpeciesC(atomicSpeciesCounter,1)
		Read(bufferC,*) pwbAtomicSpeciesDP(atomicSpeciesCounter)
		Read(bufferD,*) pwbAtomicSpeciesC(atomicSpeciesCounter,2)
	  End If	  
!Co-ordinates
      If(fileRow(1:1).ne."#".and.fileRow(1:1).ne."!")Then
	    Read(fileRowOrig,*) bufferA, bufferB, bufferC, bufferD
		If(bufferA(1:1).ne." ".and.bufferB(1:1).ne." ".and.&
		bufferC(1:1).ne." ".and.bufferD(1:1).ne." ")Then
		  atomCounter = atomCounter + 1
		  !bufferA = StrToUpper(bufferA)
		  pwbAtomLabelsInput(atomCounter) = bufferA(1:8)
		  Read(bufferB,*) pwbAtomCoordsInput(atomCounter,1)
		  Read(bufferC,*) pwbAtomCoordsInput(atomCounter,2)
		  Read(bufferD,*) pwbAtomCoordsInput(atomCounter,3)
		End If
      End If
    End Do	
	pwbNtyp = atomicSpeciesCounter
!Calc other vars
    pwbXd = ((pwbXCopy*pwbUnitVector(1,1)*pwbLatticeParameter)**2+&
             (pwbYCopy*pwbUnitVector(1,2)*pwbLatticeParameter)**2+&
	         (pwbZCopy*pwbUnitVector(1,3)*pwbLatticeParameter)**2)**0.5
	pwbYd = ((pwbXCopy*pwbUnitVector(2,1)*pwbLatticeParameter)**2+&
             (pwbYCopy*pwbUnitVector(2,2)*pwbLatticeParameter)**2+&
	         (pwbZCopy*pwbUnitVector(2,3)*pwbLatticeParameter)**2)**0.5		
	pwbZd = ((pwbXCopy*pwbUnitVector(3,1)*pwbLatticeParameter)**2+&
             (pwbYCopy*pwbUnitVector(3,2)*pwbLatticeParameter)**2+&
	         (pwbZCopy*pwbUnitVector(3,3)*pwbLatticeParameter)**2)**0.5	 
  End Subroutine readPWBatchConf
!------------------------------------------------------------------
! Prep atom config 
!------------------------------------------------------------------  
  Subroutine prepPWBatchConf() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,m,n
!blank arrays
	pwbAtomLabels = "#BLANK##"
	pwbAtomCoords = -2.1D20
	pwbAtomLabelsWorking = "#BLANK##"
	pwbAtomCoordsWorking = -2.1D20
!Make array
    n = 0
    Do i=1,pwbXCopy
      Do j=1,pwbYCopy
        Do k=1,pwbZCopy
		  Do m=1,1024
		    If(pwbAtomLabelsInput(m).eq."#BLANK##".or.pwbAtomCoordsInput(m,1).lt.-2.1D20)Then
			  Exit
			End If
		    n = n + 1
		    pwbAtomLabels(n) = pwbAtomLabelsInput(m)
		    pwbAtomCoords(n,1) = ((i-1)+pwbAtomCoordsInput(m,1))/(1.0D0*pwbXCopy)
		    pwbAtomCoords(n,2) = ((j-1)+pwbAtomCoordsInput(m,2))/(1.0D0*pwbYCopy)
		    pwbAtomCoords(n,3) = ((k-1)+pwbAtomCoordsInput(m,3))/(1.0D0*pwbZCopy)
		  End Do
		End Do
	  End Do
    End Do	  
	pwbNat = n
!Store working labels and coords
    pwbAtomCoordsWorking = pwbAtomCoords
    pwbAtomLabelsWorking = pwbAtomLabels
  End Subroutine prepPWBatchConf  
!------------------------------------------------------------------
! [FULL] 
! 
! Make pwscf batch files 
!------------------------------------------------------------------    
  Subroutine makePWBatchFiles()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, n
	Integer(kind=StandardInteger) :: nodes, ppn, processes, numa, ompThreads
	Real(kind=DoubleReal), Dimension(1:60) :: hStrain
	Real(kind=DoubleReal), Dimension(1:10) :: oStrain
	Real(kind=DoubleReal), Dimension(1:10) :: mStrain
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainArray
	Real(kind=DoubleReal) :: x, xV, y, yV, z, zV, T
	Character(len=255) :: pwbFilesOutDir, mpiRun, apRun
	
	If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Preparing from "//trim(adjustl(IntToString(pwbStart)))//&
	  " to "//trim(adjustl(IntToString(pwbEnd)))
	End If
    If(mpiProcessID.eq.0)Then
	
!Prepare variables
    Read(pwbNodes,*) nodes
    Read(pwbPPN,*) ppn
	Read(pwbOpenMPThreads,*) ompThreads
	Read(pwbProcPerNuma,*) numa
	processes = nodes * ppn   
	
!----------------------------------------------------------	
! Init batch script files
!----------------------------------------------------------	
    pwbFilesOutDir = trim(outputDirectory)//"/batch/"//trim(pwbOutputDir)
	Call makeDir(pwbFilesOutDir)
!open files for writing
    open(unit=1021,file=trim(pwbFilesOutDir)//"/"//trim("bbBatch.sh"))	
    open(unit=1022,file=trim(pwbFilesOutDir)//"/"//trim("bbBatch.pbs"))	
    open(unit=1023,file=trim(pwbFilesOutDir)//"/"//trim("all.conf"))
	
!Blue Bear MSUB
	write(1021,"(A)") '#!/bin/bash'
	write(1021,"(A)") '#MOAB -l "nodes='//trim(pwbNodes)//':ppn='//&
	trim(pwbPPN)//',walltime='//trim(pwbWallTime)//'"'
	write(1021,"(A)") '#MOAB -j oe'
	write(1021,"(A)") '#MOAB -q '//trim(pwbMoabQueue)
	write(1021,"(A)") 'OMP_NUM_THREADS='//pwbOpenMPThreads
	write(1021,"(A)") 'export OMP_NUM_THREADS'
	write(1021,"(A)") 'cd "$HOME/pwscf/test\"'
	write(1021,"(A)") 'module load apps/intel/v2013.0.079'
	write(1021,"(A)") 'module load apps/QE/v5.0.2'
	
!Blue Bear PBS
	write(1022,"(A)") '#!/bin/bash --login'
	write(1022,"(A)") '#PBS -N pwscf_job'
	write(1022,"(A)") '#PBS -l mppwidth='//trim(pwbNodes)
	write(1022,"(A)") '#PBS -l mppnppn='//trim(pwbPPN)
	write(1022,"(A)") '#PBS -l walltime='//trim(pwbWallTime)
	write(1022,"(A)") '#PBS -A '//trim(pwbPbsAccount)
	write(1022,"(A)") 'OMP_NUM_THREADS='//pwbOpenMPThreads
	write(1022,"(A)") 'export OMP_NUM_THREADS'
	write(1022,"(A)") 'export PATH=$HOME/bin:$PATH'
	write(1022,"(A)") 'export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)'
	write(1022,"(A)") 'module load acml'
	write(1022,"(A)") 'cd $PBS_O_WORKDIR'
	
!Config file
	write(1023,"(A)") '! Configuration file'
	close(1023)
	
!Run commands
	mpiRun = 'mpirun '
	apRun = 'aprun -n '//&
	trim(intToString(processes))//&
	' -N '//trim(intToString(ppn))//' -d '//&
	trim(pwbOpenMPThreads)//' -S '//&
	trim(intToString(numa ))


!----------------------------------------------------------	
! Make 1101 to 1102		Homogeneous strain
!----------------------------------------------------------	
	
! 1101 vc-relax fixed positions of atoms, relax cell size only	
    pwbUnitVectorWorking = pwbUnitVector
	pwbCalculation = "vc-relax"
	pwbFixedAtoms = 1
	If(pwbStart.le.1101.and.pwbEnd.ge.1101)Then
	  Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	  TrimSpaces(intToString(1101))//".in")
	  write(1021,"(A)") mpiRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".out"
	  write(1022,"(A)") apRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".out"
	End If

! 1102 vc-relax, relax positions of atomsa and cell size		
	pwbUnitVectorWorking = pwbUnitVector
	pwbCalculation = "vc-relax"
	pwbFixedAtoms = 0
	If(pwbStart.le.1102.and.pwbEnd.ge.1102)Then
	  Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	  TrimSpaces(intToString(1102))//".in")
	  write(1021,"(A)") mpiRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1102))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1102))//".out"
	  write(1022,"(A)") apRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1102))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1102))//".out"
	End If

! 1103 vc-relax fixed positions of atoms, relax cell size only		
	pwbUnitVectorWorking = pwbUnitVector
	pwbCalculation = "scf"
	pwbFixedAtoms = 0
	If(pwbStart.le.1103.and.pwbEnd.ge.1103)Then
	  Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	  TrimSpaces(intToString(1103))//".in")
	  write(1021,"(A)") mpiRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1103))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1103))//".out"
	  write(1022,"(A)") apRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1103))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1103))//".out"
	End If

!----------------------------------------------------------	
! Make 1201 to 1260		Homogeneous strain
!----------------------------------------------------------	
    Do i=0,22 !1-23
	  hStrain(1+i) = 0.50D0+i*0.02D0
	End Do
    Do i=0,4  !24-28
	  hStrain(24+i) = 0.95D0+i*0.01D0
	End Do
    Do i=0,4 !29-33
	  hStrain(29+i) = 1.01D0+i*0.01D0
	End Do
    Do i=0,26 !34-60
	  hStrain(34+i) = 1.10D0+i*0.02D0
	End Do
!Loop through strains
    Do i=1,60 
      pwbUnitVectorWorking = pwbUnitVector
!Make strain array
	  strainArray = 0.0D0
	  strainArray(1,1) = hStrain(i)
	  strainArray(2,2) = hStrain(i)
	  strainArray(3,3) = hStrain(i)	  
!Make working vector 
      pwbUnitVectorWorking = matmul(strainArray,pwbUnitVectorWorking)
!Write pwscf input file
      If(pwbStart.le.(1200+i).and.pwbEnd.ge.(1200+i))Then
	    Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	    TrimSpaces(intToString(1200+i))//".in")
!Write to batch files
        write(1021,"(A)") mpiRun//" pw.x < "//&
	    "pwscfbatch"//TrimSpaces(intToString(1200+i))//".in > "//&
	    "pwscfbatch"//TrimSpaces(intToString(1200+i))//".out"
	    write(1022,"(A)") apRun//" pw.x < "//&
	    "pwscfbatch"//TrimSpaces(intToString(1200+i))//".in > "//&
	    "pwscfbatch"//TrimSpaces(intToString(1200+i))//".out"
	  End If
	

!mpirun pw.x < atom.in > atom.out
	End Do
!----------------------------------------------------------
! Make 1301 to 1320		Orthorhombic strain	
!----------------------------------------------------------
	Do i=0,9 
	  oStrain(1+i) = -0.1D0+i*0.01D0
	End Do
    Do i=0,9 
	  oStrain(11+i) = 0.01D0+i*0.01D0
	End Do	
	Do i=1,20 
      pwbUnitVectorWorking = pwbUnitVector
!Make strain array
	  strainArray = 0.0D0
	  strainArray(1,1) = 1.0D0+oStrain(i)
	  strainArray(2,2) = 1.0D0-oStrain(i)
	  strainArray(3,3) = 1.0D0+(oStrain(i)*oStrain(i))/&
	                     (1-oStrain(i)*oStrain(i))
!Make working vector
      pwbUnitVectorWorking = matmul(strainArray,pwbUnitVectorWorking)
	  
      If(pwbStart.le.(1300+i).and.pwbEnd.ge.(1300+i))Then
	    Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	    TrimSpaces(intToString(1300+i))//".in")
!Write to batch files
        write(1021,"(A)") "mpirun pw.x < "//&
	    "pwscfbatch"//TrimSpaces(intToString(1300+i))//".in > "//&
	    "pwscfbatch"//TrimSpaces(intToString(1300+i))//".out"
	  End If
	End Do
!----------------------------------------------------------	
! Make 1401 to 1410		Monoclinic strain
!----------------------------------------------------------	
	Do i=1,10 
	  mStrain(i) = 0.0D0+i*0.01D0
	End Do	
	Do i=1,10 
!load input unit vector
      pwbUnitVectorWorking = pwbUnitVector
!Make strain array
	  strainArray = 0.0D0
	  strainArray(1,1) = 1.0D0+mStrain(i)
	  strainArray(2,2) = 1.0D0+mStrain(i)
	  strainArray(3,3) = (1.0D0+mStrain(i))**(-2.0D0)
!Make working vector
      pwbUnitVectorWorking = matmul(strainArray,pwbUnitVectorWorking)	  
      If(pwbStart.le.(1400+i).and.pwbEnd.ge.(1400+i))Then
	    Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	    TrimSpaces(intToString(1400+i))//".in")
!Write to batch files
        write(1021,"(A)") "mpirun pw.x < "//&
	    "pwscfbatch"//TrimSpaces(intToString(1400+i))//".in > "//&
	    "pwscfbatch"//TrimSpaces(intToString(1400+i))//".out"
	  End If
	End Do
!----------------------------------------------------------	
! Make 1501 to 1544		Specific interstitials and vacancies
!----------------------------------------------------------		
	
	

!----------------------------------------------------------	
! Make 1601 to 1610		Random triclinic strains
!----------------------------------------------------------	
		

!----------------------------------------------------------	
! Make 1701 to 1714		Fixed positions, increasing pressures
!----------------------------------------------------------	    

	
	
!----------------------------------------------------------	
! Make 1801 to 1899		Optimum cell alat, randomly varied positions
!----------------------------------------------------------	  
!Loop through configurations
	Do i=1,99
!load input unit vector
      pwbUnitVectorWorking = pwbUnitVector
	  If(pwbStart.le.(1800+i).and.pwbEnd.ge.(1800+i))Then
!Set labels and co-ords
	    pwbAtomLabelsWorking = pwbAtomLabels
	    pwbAtomCoordsWorking = pwbAtomCoords
!loop through atoms/positions
	    Do n=1,4096
		  If(pwbAtomLabelsWorking(n).eq."#BLANK##".or.pwbAtomCoordsWorking(n,1).lt.-2.1D20)Then
			    Exit
		  End If
!Vary co-ordinate
          x = pwbAtomCoordsWorking(n,1)
          y = pwbAtomCoordsWorking(n,2)
          z = pwbAtomCoordsWorking(n,3)
		  xV = pwbVaryPoint(x, pwbXd, 0.3D0, 0.1D0)
		  !yV = pwbVaryPoint(y, pwbYd, 0.3D0, 0.1D0)
		  !zV = pwbVaryPoint(z, pwbZd, 0.3D0, 0.1D0)
		  
          xV = VaryPointRand(x,50.0D0,0.01D0/(1.0D0*pwbXCopy),0.0D0,0.0D0,0)
          yV = VaryPointRand(y,50.0D0,0.01D0/(1.0D0*pwbYCopy),0.0D0,0.0D0,0)
          zV = VaryPointRand(z,50.0D0,0.01D0/(1.0D0*pwbZCopy),0.0D0,0.0D0,0)
!Store variations
          pwbAtomCoordsWorking(n,1) = xV
          pwbAtomCoordsWorking(n,2) = yV
          pwbAtomCoordsWorking(n,3) = zV
	    End Do
!Make file
	    Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	    TrimSpaces(intToString(1800+i))//".in")
!Write to batch files
        write(1021,"(A)") "mpirun pw.x < "//&
	    "pwscfbatch"//TrimSpaces(intToString(1800+i))//".in > "//&
	    "pwscfbatch"//TrimSpaces(intToString(1800+i))//".out"
        write(1022,"(A)") "mpirun pw.x < "//&
	    "pwscfbatch"//TrimSpaces(intToString(1800+i))//".in > "//&
	    "pwscfbatch"//TrimSpaces(intToString(1800+i))//".out"
	  End If	
	End Do
	
!Close files
    Close(1021)
    Close(1022)
	
	End If
  End Subroutine makePWBatchFiles  
  
  
!------------------------------------------------------------------
! [SLAB] 
! 
! Make pwscf batch files 
!------------------------------------------------------------------    
  Subroutine makePWSlab()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, n
	Integer(kind=StandardInteger) :: nodes, ppn, processes, numa, ompThreads
	Real(kind=DoubleReal), Dimension(1:60) :: hStrain
	Real(kind=DoubleReal), Dimension(1:10) :: oStrain
	Real(kind=DoubleReal), Dimension(1:10) :: mStrain
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainArray
	Real(kind=DoubleReal) :: x, xV, y, yV, z, zV, T
	Character(len=255) :: pwbFilesOutDir, mpiRun, apRun	
!----------------------------------------------------------	
! Make directory
!----------------------------------------------------------	
    pwbFilesOutDir = trim(outputDirectory)//"/batch/"//trim(pwbOutputDir)
	Call makeDir(pwbFilesOutDir)
!----------------------------------------------------------	
! Make 1101 to 1102		Homogeneous strain
!----------------------------------------------------------	
! 1101 vc-relax fixed positions of atoms, relax cell size only	
    pwbUnitVectorWorking = pwbUnitVector
	pwbCalculation = "vc-relax"
	pwbFixedAtoms = 1
	If(pwbStart.le.1101.and.pwbEnd.ge.1101)Then
	  Call writePWscfFile(trim(pwbFilesOutDir),"pwscfbatch"//&
	  TrimSpaces(intToString(1101))//".in")
	  write(1021,"(A)") mpiRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".out"
	  write(1022,"(A)") apRun//" pw.x < "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".in > "//&
	  "pwscfbatch"//TrimSpaces(intToString(1101))//".out"
	End If
  End Subroutine makePWSlab
  
!------------------------------------------------------------------
! [BASIC] 
! 
! Make pwscf batch files 
!------------------------------------------------------------------    
  Subroutine makePWBasic()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k, n
	Integer(kind=StandardInteger) :: nodes, ppn, processes, numa, ompThreads
	Real(kind=DoubleReal), Dimension(1:60) :: hStrain
	Real(kind=DoubleReal), Dimension(1:10) :: oStrain
	Real(kind=DoubleReal), Dimension(1:10) :: mStrain
	Real(kind=DoubleReal), Dimension(1:3,1:3) :: strainArray
	Real(kind=DoubleReal) :: x, xV, y, yV, z, zV, T
	Character(len=255) :: pwbFilesOutDir, mpiRun, apRun	
!----------------------------------------------------------	
! Make directory
!----------------------------------------------------------	
    pwbFilesOutDir = trim(outputDirectory)//"/batch/"//trim(pwbOutputDir)
	Call makeDir(pwbFilesOutDir)
!----------------------------------------------------------	
! Make basic file
!----------------------------------------------------------	
    pwbUnitVectorWorking = pwbUnitVector
	pwbCalculation = "vc-relax"
	pwbFixedAtoms = 0
	Call writePWscfFile(trim(pwbFilesOutDir),"pwscf.in")
  End Subroutine makePWBasic
  
  
  
!------------------------------------------------------------------
! Write pwscf file 
!------------------------------------------------------------------  
  Subroutine writePWscfFile(outputFilePath, outputFileName)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8  
	Integer(kind=StandardInteger) :: i
	Real(kind=DoubleReal) :: a, b, c, cosBC, cosAC, cosAB
	Character(*) :: outputFileName, outputFilePath 
	Character(len=128) :: outputFileWrite, outputFileAllConfigs	
!calculate variables
    a = 1.0D0*sqrt(pwbUnitVectorWorking(1,1)**2+pwbUnitVectorWorking(1,2)**2+&
	    pwbUnitVectorWorking(1,3)**2)*&
		pwbXCopy*pwbLatticeParameter
	b = 1.0D0*(sqrt(pwbUnitVectorWorking(2,1)**2+pwbUnitVectorWorking(2,2)**2+&
	    pwbUnitVectorWorking(2,3)**2)*&
		pwbXCopy*pwbLatticeParameter)/a	
	c = 1.0D0*(sqrt(pwbUnitVectorWorking(3,1)**2+pwbUnitVectorWorking(3,2)**2+&
	    pwbUnitVectorWorking(3,3)**2)*&
		pwbXCopy*pwbLatticeParameter)/a		
	cosBC = 1.0D0*(pwbUnitVectorWorking(2,1)*pwbUnitVectorWorking(3,1)+&
	        pwbUnitVectorWorking(2,2)*pwbUnitVectorWorking(3,2)+&
	        pwbUnitVectorWorking(2,3)*pwbUnitVectorWorking(3,3))/(1.0D0*b*c)
	cosAC = 1.0D0*(pwbUnitVectorWorking(1,1)*pwbUnitVectorWorking(3,1)+&
	        pwbUnitVectorWorking(1,2)*pwbUnitVectorWorking(3,2)+&
	        pwbUnitVectorWorking(1,3)*pwbUnitVectorWorking(3,3))/(1.0D0*a*c)
	cosAB = 1.0D0*(pwbUnitVectorWorking(1,1)*pwbUnitVectorWorking(2,1)+&
	        pwbUnitVectorWorking(1,2)*pwbUnitVectorWorking(2,2)+&
	        pwbUnitVectorWorking(1,3)*pwbUnitVectorWorking(2,3))/(1.0D0*a*c)
!Convert from Angstrom to Bohr
    a = UnitConvert(a,"ANGS","BOHR")
  
!Open output files
	outputFileWrite = trim(outputFilePath)//"/"//trim(outputFileName)
	open(unit=103,file=trim(outputFileWrite))
	
!---------------------------
! Control
!---------------------------
    write(103,"(A8)") "&CONTROL"
    write(103,"(A)") "restart_mode = '"//TrimSpaces(pwbRestartMode)//"',"
    write(103,"(A)") "calculation = '"//TrimSpaces(pwbCalculation)//"',"
    write(103,"(A)") "outdir = '"//TrimSpaces(pwbOutDir)//"',"
    write(103,"(A)") "pseudo_dir = '"//TrimSpaces(pwbPseudoDir)//"',"
    write(103,"(A)") "prefix = '"//TrimSpaces(pwbPrefix)//"',"
    write(103,"(A)") "etot_conv_thr = "//TrimSpaces(dpToString(pwbEtotConvThr))//","
    write(103,"(A)") "forc_conv_thr = "//TrimSpaces(dpToString(pwbForcConvThr))//","
    write(103,"(A)") "nstep = "//TrimSpaces(intToString(pwbNstep))//","
    write(103,"(A)") "tprnfor = "//TrimSpaces(pwbTprnfor)//","
    write(103,"(A)") "tstress = "//TrimSpaces(pwbTstress)//","
    write(103,"(A)") "disk_io = '"//TrimSpaces(pwbDiskIO)//"',"
	write(103,"(A1)") "/" 
!---------------------------
! System
!---------------------------
    write(103,"(A7)") "&SYSTEM"  
    write(103,"(A)") "ibrav = "//TrimSpaces(intToString(pwbIbrav))//","
!----Unit cell 6 parameters----! 
    write(103,"(A)") "celldm(1) = "//TrimSpaces(dpToString(a))//","
	write(103,"(A)") "celldm(2) = "//TrimSpaces(dpToString(b))//","
	write(103,"(A)") "celldm(3) = "//TrimSpaces(dpToString(c))//","
	write(103,"(A)") "celldm(4) = "//TrimSpaces(dpToString(cosBC))//","
	write(103,"(A)") "celldm(5) = "//TrimSpaces(dpToString(cosAC))//","
	write(103,"(A)") "celldm(6) = "//TrimSpaces(dpToString(cosAB))//","

!----Unit cell 6 parameters----! 
    write(103,"(A)") "nat = "//TrimSpaces(intToString(pwbNat))//","
    write(103,"(A)") "ntyp = "//TrimSpaces(intToString(pwbNtyp))//","
    write(103,"(A)") "ecutwfc = "//TrimSpaces(intToString(pwbEcutwfc))//","
    write(103,"(A)") "ecutrho = "//TrimSpaces(intToString(pwbEcutrho))//","
    write(103,"(A)") "occupations = '"//TrimSpaces(pwbOccupations)//"',"
    write(103,"(A)") "smearing = '"//TrimSpaces(pwbSmearing)//"',"
    write(103,"(A)") "degauss = "//TrimSpaces(dpToString(pwbDegauss))//","
	write(103,"(A1)") "/" 
!---------------------------
! Electrons
!---------------------------
    write(103,"(A10)") "&ELECTRONS"  
    write(103,"(A)") "diagonalization = '"//TrimSpaces(pwbDiagonalization)//"',"
    write(103,"(A)") "mixing_mode = "//TrimSpaces(pwbMixingMode)//","
    write(103,"(A)") "mixing_beta = "//TrimSpaces(dpToString(pwbMixingBeta))//","
    write(103,"(A)") "conv_thr = "//TrimSpaces(dpToString(pwbConvThr))//","
	write(103,"(A1)") "/" 
!---------------------------
! Ions
!---------------------------
    write(103,"(A5)") "&IONS"  
    write(103,"(A)") "ion_dynamics = '"//TrimSpaces(pwbIonDynamics)//"',"
	write(103,"(A1)") "/" 
!---------------------------
! Cell
!---------------------------
    write(103,"(A5)") "&CELL"  
    write(103,"(A)") "cell_dynamics = '"//TrimSpaces(pwbCellDynamics)//"',"
    write(103,"(A)") "press = "//TrimSpaces(dpToString(pwbPress))//","
    write(103,"(A)") "cell_factor = "//TrimSpaces(dpToString(pwbCellFactor))//","
	write(103,"(A1)") "/" 
!---------------------------
! Atomic Species
!---------------------------
    write(103,"(A14)") "ATOMIC_SPECIES"  
    Do i=1,pwbNtyp
      write(103,"(A)") trim(pwbAtomicSpeciesC(i,1))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomicSpeciesDP(i)))//"  "//&
					   trim(pwbAtomicSpeciesC(i,2))
	End Do
!---------------------------
! Atomic Positions
!---------------------------
    If(pwbFixedAtoms.eq.0)Then	!don't fix atoms, set force to 1 (default)
      write(103,"(A24)") "ATOMIC_POSITIONS crystal"  
      Do i=1,pwbNat
        write(103,"(A)") trim(pwbAtomLabelsWorking(i))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomCoordsWorking(i,1)))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomCoordsWorking(i,2)))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomCoordsWorking(i,3)))//"  "
					   !TrimSpaces(dpToString(1.0D0))//"  "//&
					   !TrimSpaces(dpToString(1.0D0))//"  "//&
					   !TrimSpaces(dpToString(1.0D0))
	  End Do
	End If
    If(pwbFixedAtoms.eq.1)Then	!don't fix atoms, set force to 1 (default)
      write(103,"(A24)") "ATOMIC_POSITIONS crystal"  
      Do i=1,pwbNat
        write(103,"(A)") trim(pwbAtomLabelsWorking(i))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomCoordsWorking(i,1)))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomCoordsWorking(i,2)))//"  "//&
	                   TrimSpaces(dpToString(pwbAtomCoordsWorking(i,3)))//"  "
					   !TrimSpaces(dpToString(0.0D0))//"  "//&
					   !TrimSpaces(dpToString(0.0D0))//"  "//&
					   !TrimSpaces(dpToString(0.0D0))
	  End Do
	End If
!---------------------------
! K-Points
!---------------------------
    !write(103,"(A18)") "K_POINTS automatic" 
    !write(103,"(A)") TrimSpaces(pwbKpoints)
	write(103,"(A14)") "K_POINTS gamma" 
!Close file
    close(103)  
  	!pwbAtomLabelsInput = "#BLANK##"
	!pwbAtomCoordsInput = -2.1D20
	!pwbAtomicSpeciesC = "#BLANK##"
	!pwbAtomicSpeciesDP = -2.1D20  
  	!pwbRestartMode = "from_scratch"
	!pwbCalculation = "scf"
	!pwbDiskIO = "low"
	!pwbOccupations = "smearing"
	!pwbSmearing = "mv"
    !pwbEtotConvThr = 1.0D-7
	
  !Config file
    !open(unit=1023,file=trim(pwbOutDir)//"/"//trim("all.conf"))
	
	!outputFileAllConfigs = trim(outputFilePath)//"/all.conf"
	!open(unit=1023,file=trim(outputFileAllConfigs),status="old",position="append",action="write")
	!write(1023,"(A)") '#new'
	!write(1023,"(A)") '#LP '//TrimSpaces(dpToString(a))
	
	

!
!#LP 4.035
!#X 1 0 0
!#Y 0 1 0
!#Z 0 0 1
!#CC 4 4 4
!#RC 6.5                
!#F N         
!#BM 76 GPA
!#EC 3 114 61 31 GPA 
!#EV 4204 ANG3
!#CW 1 1 1 1
	
	write(1023,"(A)") '#end'
	write(1023,"(A)") '!------------------------------'
	close(1023)
  
	
  End Subroutine writePWscfFile
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 

  Function pwbVaryPoint(inputCoord, aLat, sigma, maxVariation) RESULT (outputCoord)
!force declaration of all variables
	Implicit None
!declare variables  
    Real(kind=DoubleReal) :: inputCoord, outputCoord
    Real(kind=DoubleReal) :: aLat, sigma, maxVariation, pointVariation, cellVariation
    Real(kind=DoubleReal) :: randNumber, addSubtract
!convert maximum variation from ang to lattice coords
	cellVariation = 1.0D0 * (maxVariation / aLat)
!Random float, Gaussian type distribution
	randNumber = RandomFloat(0.0D0,1.0D0,"G",sigma)
	pointVariation = randNumber * cellVariation  
    Call RANDOM_NUMBER(addSubtract)
	If(addSubtract.le.0.5D0)Then
	  outputCoord = inputCoord + pointVariation
	Else
	  outputCoord = inputCoord - pointVariation
	End If
  End Function pwbVaryPoint
  
  
End Module pwbatch  