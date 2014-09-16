Module pwBatch

!--------------------------------------------------------------!
! General subroutines and functions                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

! Read user input file 

!----------------------------------------
! Updated: 12th Aug 2014
!----------------------------------------

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
!Include MPI header
  Include 'mpif.h'  

  
  
  
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: runPWBatch 
  
  
Contains
  Subroutine runPWBatch()
    Implicit None   ! Force declaration of all variables
! Private variables    

! Output to Terminal
    If(mpiProcessID.eq.0.and.printToTerminal.eq.1)Then
      Print *,"PWscf Batch Files ",trim(pwbConfigFilePath)
    End If      
! Load PWscf Batch Variables from User Input File
    Call loadPWBatchVars()
! Load PWscf Batch initial configuration
    Call loadPWBatchConfig()
! Make batch files
    Call makeBatchFiles()
    
! Prepare the temporary config file
    If(mpiProcessID.eq.0)Then
    
    End If
! Synch MPI processes    
    Call M_synchProcesses() 
  End Subroutine runPWBatch 
   
!----------------------------------------------------------
! Load batch variables from user input file
!----------------------------------------------------------
  Subroutine loadPWBatchVars()
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i
    Character(len=255) :: fileRow
!re-read user input file for specific pwbatch entries to replace defaults    
!open & read in file    
    Open(UNIT=1,FILE=Trim(inputFilePathT)) 
    Do i=1,maxFileRows 
!Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
      If (ios /= 0) Then
        EXIT 
      End If
      If(StrToUpper(fileRow(1:10)).eq."#PWBOUTDIR")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        pwbOutDir = trim(adjustl(fileRow))
      End If
      If(StrToUpper(fileRow(1:13)).eq."#PWBPSEUDODIR")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        pwbPseudoDir = trim(adjustl(fileRow))
      End If
      If(StrToUpper(fileRow(1:11)).eq."#PWBECUTWFC")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        Read(fileRow,*) pwbEcutwfc
      End If
      If(StrToUpper(fileRow(1:11)).eq."#PWBECUTRHO")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        Read(fileRow,*) pwbEcutrho
      End If
      If(StrToUpper(fileRow(1:14)).eq."#PWBMIXINGMODE")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        Read(fileRow,*) pwbMixingMode
      End If
      If(StrToUpper(fileRow(1:15)).eq."#PWBCALCULATION")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        Read(fileRow,*) pwbCalculation
      End If
      If(StrToUpper(fileRow(1:11)).eq."#PWBKPOINTS")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        pwbKpoints = StrToUpper(Trim(Adjustl(fileRow)))
      End If
      If(StrToUpper(fileRow(1:8)).eq."#PWBNBND")then
!read next line
        Read(1,"(A255)",IOSTAT=ios) fileRow
        Read(fileRow,*) pwbNbnd
      End If
      
      
      
    End Do  
  End Subroutine loadPWBatchVars 
  
  
   
!----------------------------------------------------------
! Load initial configuration from file
!----------------------------------------------------------  
  Subroutine loadPWBatchConfig()
    Implicit None   ! Force declaration of all variables
! Private variables   
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000 
    Integer(kind=StandardInteger) :: ios, i, j, k, m, n
    Integer(kind=StandardInteger) :: atomCounter, atomicSpeciesCounter
    Character(len=255) :: fileRow, fileRowOrig
    Character(len=64) :: bufferA, bufferB, bufferC, bufferD
! Open configuration file
    Open(unit=101,file=trim(pwbConfigFilePath))   
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
        Read(bufferB,*) pwbAtomicSpeciesL(atomicSpeciesCounter)
        Read(bufferC,*) pwbAtomicSpeciesDP(atomicSpeciesCounter)
        Read(bufferD,*) pwbAtomicSpeciesPP(atomicSpeciesCounter)
        pwbAtomicSpeciesL(atomicSpeciesCounter) = &
        trim(adjustl(pwbAtomicSpeciesL(atomicSpeciesCounter)))
        pwbAtomicSpeciesPP(atomicSpeciesCounter) = &
        trim(adjustl(pwbAtomicSpeciesPP(atomicSpeciesCounter)))
      End If      
!Co-ordinates
      If(fileRow(1:1).ne."#".and.fileRow(1:1).ne."!")Then
        Read(fileRowOrig,*) bufferA, bufferB, bufferC, bufferD
        If(bufferA(1:1).ne." ".and.bufferB(1:1).ne." ".and.&
        bufferC(1:1).ne." ".and.bufferD(1:1).ne." ")Then
          atomCounter = atomCounter + 1
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
! Close file
    Close(101)
! Expand co-ordinates
! Blank arrays
    pwbAtomLabels = "#BLANK##"
    pwbAtomCoords = -2.1D20
    pwbAtomLabelsWorking = "#BLANK##"
    pwbAtomCoordsWorking = -2.1D20
! Make array
    n = 0
    Do i=1,pwbXCopy
      Do j=1,pwbYCopy
        Do k=1,pwbZCopy
          Do m=1,1024
            If(pwbAtomLabelsInput(m).eq."#BLANK##".or.pwbAtomCoordsInput(m,1).lt.-2.1D20)Then
              Exit
            End If
            n = n + 1
            pwbAtomLabels(n) = Trim(adjustl(pwbAtomLabelsInput(m)))
            pwbAtomCoords(n,1) = ((i-1)+pwbAtomCoordsInput(m,1))/(1.0D0*pwbXCopy)
            pwbAtomCoords(n,2) = ((j-1)+pwbAtomCoordsInput(m,2))/(1.0D0*pwbYCopy)
            pwbAtomCoords(n,3) = ((k-1)+pwbAtomCoordsInput(m,3))/(1.0D0*pwbZCopy)
          End Do
        End Do
      End Do
    End Do      
    pwbNat = n
! Store working labels and coords
    pwbAtomCoordsWorking = pwbAtomCoords
    pwbAtomLabelsWorking = pwbAtomLabels
  End Subroutine loadPWBatchConfig 
  
!------------------------------------------------------------------
! Write pwscf file 
!------------------------------------------------------------------  
  Subroutine makeBatchFiles() 
!force declaration of all variables
    Implicit None    
!declare private variables
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: maxVariance
!-----------------------
! Single "batch" file
!-----------------------
    If(eampaRunType(1:4).eq."PWB1")Then ! Make single batch file
      pwbUnitVectorWorking = pwbUnitVector
      If(mpiProcessID.eq.0)Then
        Call writePWscfFile("pwscf.in")
      End If
! Synch MPI processes    
      Call M_synchProcesses() 
    End If
!-----------------------
! Single "batch" file, randomised coords
!-----------------------
    If(eampaRunType(1:4).eq."PWB2")Then ! Make single batch file
! Set unit vector
      pwbUnitVectorWorking = pwbUnitVector
! Calculate primary lattice parameter     
      maxVariance = pwbVarianceMax/(1.0D0*sqrt(pwbUnitVectorWorking(1,1)**2+&
      pwbUnitVectorWorking(1,2)**2+&
      pwbUnitVectorWorking(1,3)**2)*pwbXCopy*pwbLatticeParameter)
      Do j=1,3
        Do i=1,size(pwbAtomCoordsWorking,1)
          If(pwbAtomLabelsWorking(j).eq."#BLANK##")Then
            Exit
          End If
          pwbAtomCoordsWorking(i,j) = &
          RandomVaryPoint(pwbAtomCoordsWorking(i,j), maxVariance, pwbVarianceSigma)
        End Do
      End Do
      If(mpiProcessID.eq.0)Then
        Call writePWscfFile("pwscf.in")
      End If
! Synch MPI processes    
      Call M_synchProcesses() 
    End If    

    
  
  
  
  End Subroutine makeBatchFiles 
  
    
!------------------------------------------------------------------
! Write pwscf file 
!------------------------------------------------------------------  
  Subroutine writePWscfFile(outputFileName)
!force declaration of all variables
    Implicit None    
!declare private variables
    Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000  
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: a, b, c, cosBC, cosAC, cosAB, zero
    Character(*) :: outputFileName 
    Character(len=255) :: outputFileWrite    
! Initialise variables
    If(pwbBatchDir(1:1).eq." ")Then      
      outputFileWrite = trim(outputDirectory)//"/"//trim(outputFileName)
      Call makeDir(trim(outputDirectory))
    Else      
      outputFileWrite = trim(outputDirectory)//"/"//trim(pwbBatchDir)//&
      "/"//trim(outputFileName)
      Call makeDir(trim(outputDirectory)//"/"//trim(pwbBatchDir))
    End If    
    zero = 0.0D0
! Calculate variables    
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
!Convert from Angstrom to Bohr (Bohr is the unit used by PWscf)
    a = UnitConvert(a,"ANGS","BOHR")
!Open output file
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
    write(103,"(A12,F12.7,A1)") "celldm(1) = ",a,","
    write(103,"(A12,F12.7,A1)") "celldm(2) = ",b,","
    write(103,"(A12,F12.7,A1)") "celldm(3) = ",c,","
    write(103,"(A12,F12.7,A1)") "celldm(4) = ",cosBC,","
    write(103,"(A12,F12.7,A1)") "celldm(5) = ",cosAC,","
    write(103,"(A12,F12.7,A1)") "celldm(6) = ",cosAB,","
!----Unit cell 6 parameters----! 
    write(103,"(A)") "nat = "//TrimSpaces(intToString(pwbNat))//","
    write(103,"(A)") "ntyp = "//TrimSpaces(intToString(pwbNtyp))//","
    write(103,"(A)") "nbnd = "//TrimSpaces(intToString(pwbNbnd))//","
    write(103,"(A)") "ecutwfc = "//TrimSpaces(intToString(pwbEcutwfc))//","
    write(103,"(A)") "ecutrho = "//TrimSpaces(intToString(pwbEcutrho))//","
    write(103,"(A)") "occupations = '"//TrimSpaces(pwbOccupations)//"',"
    write(103,"(A)") "smearing = '"//TrimSpaces(pwbSmearing)//"',"
    write(103,"(A10,F5.3,A1)") "degauss = ",pwbDegauss,","
    write(103,"(A1)") "/" 
!---------------------------
! Electrons
!---------------------------
    write(103,"(A10)") "&ELECTRONS"  
    write(103,"(A)") "diagonalization = '"//TrimSpaces(pwbDiagonalization)//"',"
    write(103,"(A)") "mixing_mode = '"//TrimSpaces(pwbMixingMode)//"',"
    write(103,"(A14,F5.3,A1)") "mixing_beta = ",pwbMixingBeta,","
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
      write(103,"(A8,A2,F12.7,A2,A)") SpacesRight(pwbAtomicSpeciesL(i)),"  ",&
                       pwbAtomicSpeciesDP(i),"  ",&
                       trim(pwbAtomicSpeciesPP(i))
    End Do
!---------------------------
! Atomic Positions
!---------------------------
    write(103,"(A24)") "ATOMIC_POSITIONS crystal"  
    If(pwbFixedAtoms.eq.0)Then    !don't fix atoms, set force to 1 (default)
      Do i=1,pwbNat
        write(103,"(A8,A1,F12.7,A1,F12.7,A1,F12.7,A1)") &
            SpacesRight(pwbAtomLabelsWorking(i))," ",&
            pwbAtomCoordsWorking(i,1)," ",&
            pwbAtomCoordsWorking(i,2)," ",&
            pwbAtomCoordsWorking(i,3)," "
      End Do
    End If
    If(pwbFixedAtoms.eq.1)Then    !don't fix atoms, set force to 1 (default)
      Do i=1,pwbNat
        write(103,"(A8,A1,F12.7,A1,F12.7,A1,F12.7,A1,F12.7,A1,F12.7,A1,F12.7,A1)") &
            trim(pwbAtomLabelsWorking(i))," ",&
            pwbAtomCoordsWorking(i,1)," ",&
            pwbAtomCoordsWorking(i,2)," ",&
            pwbAtomCoordsWorking(i,3)," ",zero," ",zero," ",zero," "
      End Do
    End If
!---------------------------
! K-Points
!---------------------------
    If(pwbKpoints(1:5).eq."GAMMA")Then
      write(103,"(A14)") "K_POINTS gamma" 
    Else
      write(103,"(A18)") "K_POINTS automatic" 
      write(103,"(A)") TrimSpaces(pwbKpoints)
    End If
!Close file
    close(103)  

  
    
  End Subroutine writePWscfFile
  
  
  
  
  
  
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!   



  
  
End Module pwBatch  