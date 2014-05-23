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
  Use strings		!string functions
  Use maths
  Use initialise
  Use input

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!declare global variables  
  Real(kind=DoubleReal), Dimension(1:3,1:3) :: pwbUnitVector
  Real(kind=DoubleReal) :: pwbLatticeParameter
  Integer(kind=StandardInteger) :: pwbXCopy, pwbYCopy, pwbZCopy	
  Character(len=8), Dimension(1:1024)   :: pwbAtomLabelsInput
  Real(kind=DoubleReal), Dimension(1:1024,1:3)   :: pwbAtomCoordsInput  
  Character(len=8), Dimension(1:4096)   :: pwbAtomLabels, pwbAtomLabelsWorking
  Real(kind=DoubleReal), Dimension(1:4096,1:3)   :: pwbAtomCoords, pwbAtomCoordsWorking  	
  Character(len=64), Dimension(1:128,1:2)   :: pwbAtomicSpeciesC  	
  Real(kind=DoubleReal), Dimension(1:128)   :: pwbAtomicSpeciesDP
!PWscf options
  Character(len=16) :: pwbRestartMode, pwbCalculation, pwbDiskIO, pwbOccupations, pwbSmearing, &
                       pwbDiagonalization, pwbMixingMode, pwbIonDynamics, pwbCellDynamics  
  Character(len=6) :: pwbTprnfor, pwbTstress
  Real(kind=DoubleReal) :: pwbEtotConvThr, pwbForcConvThr, pwbDegauss, pwbMixingBeta, &
                           pwbConvThr, pwbPress, pwbCellFactor
  Integer(kind=StandardInteger) :: pwbNstep, pwbIbrav, pwbNat, pwbNtyp, pwbEcutwfc, pwbEcutrho
  
!Privacy of functions/subroutines/variables
  Private    
!Variables - Read in from pwbOptConf
  Public :: pwbUnitVector
  Public :: pwbLatticeParameter
  Public :: pwbXCopy, pwbYCopy, pwbZCopy	
  Public :: pwbAtomLabels, pwbAtomLabelsInput, pwbAtomLabelsWorking
  Public :: pwbAtomCoords, pwbAtomCoordsInput, pwbAtomCoordsWorking
!Variables - PWscf options 
  Public :: pwbRestartMode, pwbCalculation, pwbDiskIO, pwbOccupations, pwbSmearing, &
            pwbDiagonalization, pwbMixingMode, pwbIonDynamics, pwbCellDynamics 
  Public :: pwbTprnfor, pwbTstress
  Public :: pwbEtotConvThr, pwbForcConvThr, pwbDegauss, pwbMixingBeta, &
            pwbConvThr, pwbPress, pwbCellFactor
  Public :: pwbNstep, pwbIbrav, pwbNat, pwbNtyp, pwbEcutwfc, pwbEcutrho
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
	  print *,ProgramTime(),pwbOptConfFile
	End If
!Read in pwbOptConf files
    Call readPWBatchConf()
!prep config
    Call prepPWBatchConf()

	Call writePWscfFile("pwscfbatch.in")
	
  End Subroutine runPWBatch    
!------------------------------------------------------------------
! Read pwb opt file 
!------------------------------------------------------------------  
  Subroutine readPWBatchConf()
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8 
	Integer(kind=StandardInteger) :: ios, i, atomCounter	
	Character(len=128) :: inputFile
	Character(len=255) :: fileRow	
	Character(len=64) :: bufferA, bufferB, bufferC, bufferD
!Default values
!--------------------------------------------------
    pwbUnitVector = 0.0D0	
	pwbAtomLabelsInput = "#BLANK##"
	pwbAtomCoordsInput = -2.1D20
	pwbAtomicSpeciesC = "#BLANK##"
	pwbAtomicSpeciesD = -2.1D20
!Default pwscf file values, text
	pwbRestartMode = "from_scratch"
	pwbCalculation = "scf"
	pwbDiskIO = "low"
	pwbOccupations = "smearing"
	pwbSmearing = "mv"
	pwbDiagonalization = "david"
	pwbMixingMode = "plain"
	pwbIonDynamics = "bfgs"
	pwbCellDynamics = "bfgs"
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
	pwbNstep = 200
	pwbIbrav = 14
	pwbNat = 0
	pwbNtyp = 0
	pwbEcutwfc = 35
	pwbEcutrho = 208
!--------------------------------------------------	
!open output file	
	inputFile = trim(currentWorkingDirectory)//"/"//pwbOptConfFile
	Open(unit=101,file=trim(inputFile))   
	atomCounter = 0
	Do i=1,maxFileRows 
!Read in line
	  Read(101,"(A255)",IOSTAT=ios) fileRow
!break out if end of file
	  If (ios /= 0) Then
	    EXIT 
	  End If
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
!Co-ordinates
      If(fileRow(1:1).ne."#".and.fileRow(1:1).ne."!")Then
	    Read(fileRow,*) bufferA, bufferB, bufferC, bufferD
		If(bufferA(1:1).ne." ".and.bufferB(1:1).ne." ".and.&
		bufferC(1:1).ne." ".and.bufferD(1:1).ne." ")Then
		  atomCounter = atomCounter + 1
		  bufferA = StrToUpper(bufferA)
		  pwbAtomLabelsInput(atomCounter) = bufferA(1:8)
		  Read(bufferB,*) pwbAtomCoordsInput(atomCounter,1)
		  Read(bufferC,*) pwbAtomCoordsInput(atomCounter,2)
		  Read(bufferD,*) pwbAtomCoordsInput(atomCounter,3)
		End If
      End If
    End Do	
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
	
	!Do n=1,1024      
	!  If(pwbAtomLabels(n).ne."#BLANK##".and.mpiProcessID.eq.0)Then
	!    print *,pwbAtomLabels(n),pwbAtomCoords(n,1),pwbAtomCoords(n,2),pwbAtomCoords(n,3)
	!  End If
	!End Do
	
  End Subroutine prepPWBatchConf  
  
  
  
  
  
!------------------------------------------------------------------
! Write pwscf file 
!------------------------------------------------------------------  
  Subroutine writePWscfFile(outputFileName)
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger), Parameter :: maxFileRows = 1E8  
	Character(*) :: outputFileName
	Character(len=128) :: outputFilePath	
  
  
!Open output file	  
	outputFilePath = trim(currentWorkingDirectory)//"/"//trim(outputFileName)
	open(unit=103,file=trim(outputFilePath))
  
!---------------------------
! Control
!---------------------------
    write(103,"(A8)") "&CONTROL"
    write(103,"(A)") "restart_mode = '"//TrimSpaces(pwbRestartMode)//"',"
    write(103,"(A)") "calculation = '"//TrimSpaces(pwbCalculation)//"',"
    write(103,"(A)") "etot_conv_thr = '"//TrimSpaces(dpToString(pwbEtotConvThr))//"',"
    write(103,"(A)") "forc_conv_thr = '"//TrimSpaces(dpToString(pwbForcConvThr))//"',"
    write(103,"(A)") "nstep = '"//TrimSpaces(intToString(pwbNstep))//"',"
    write(103,"(A)") "tprnfor = '"//TrimSpaces(pwbTprnfor)//"',"
    write(103,"(A)") "tstress = '"//TrimSpaces(pwbTstress)//"',"
    write(103,"(A)") "disk_io = '"//TrimSpaces(pwbDiskIO)//"',"
	write(103,"(A1)") "/" 
!---------------------------
! System
!---------------------------
    write(103,"(A7)") "&SYSTEM"  
    write(103,"(A)") "ibrav = '"//TrimSpaces(intToString(pwbIbrav))//"',"
!----Unit cell 6 parameters----! 
    write(103,"(A)") "nat = '"//TrimSpaces(intToString(pwbNat))//"',"
    write(103,"(A)") "ntyp = '"//TrimSpaces(intToString(pwbNtyp))//"',"
    write(103,"(A)") "ecutwfc = '"//TrimSpaces(intToString(pwbEcutwfc))//"',"
    write(103,"(A)") "ecutrho = '"//TrimSpaces(intToString(pwbEcutrho))//"',"
    write(103,"(A)") "occupations = '"//TrimSpaces(pwbOccupations)//"',"
    write(103,"(A)") "smearing = '"//TrimSpaces(pwbSmearing)//"',"
    write(103,"(A)") "degauss = '"//TrimSpaces(dpToString(pwbDegauss))//"',"
	write(103,"(A1)") "/" 
!---------------------------
! Electrons
!---------------------------
    write(103,"(A10)") "&ELECTRONS"  
    write(103,"(A)") "diagonalization = '"//TrimSpaces(pwbDiagonalization)//"',"
    write(103,"(A)") "mixing_mode = '"//TrimSpaces(pwbMixingMode)//"',"
    write(103,"(A)") "mixing_beta = '"//TrimSpaces(dpToString(pwbMixingBeta))//"',"
    write(103,"(A)") "conv_thr = '"//TrimSpaces(dpToString(pwbConvThr))//"',"
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
    write(103,"(A)") "press = '"//TrimSpaces(dpToString(pwbPress))//"',"
    write(103,"(A)") "cell_factor = '"//TrimSpaces(dpToString(pwbCellFactor))//"',"
	write(103,"(A1)") "/" 
  
  
  
  	!pwbRestartMode = "from_scratch"
	!pwbCalculation = "scf"
	!pwbDiskIO = "low"
	!pwbOccupations = "smearing"
	!pwbSmearing = "mv"
    !pwbEtotConvThr = 1.0D-7
  
  
    close(103)
	
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

  
  
End Module pwbatch  