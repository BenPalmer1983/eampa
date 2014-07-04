Module preppotfit

!--------------------------------------------------------------!
! Prep EAM File Subroutines                        
! Ben Palmer, University of Birmingham   
!--------------------------------------------------------------!

!----------------------------------------
! Updated: 1st May 2014
!----------------------------------------

! Setup Modules
  Use kinds
  Use constants
  Use general
  Use mpif
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use output
  Use prep

!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
!Privacy of functions/subroutines/variables
  Private    
!Public subroutines
  Public :: runPrepPotfit
  
  
  
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


  Subroutine runPrepPotfit() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
    Character(len=255) :: filePath
!create potfit folder
	filePath = trim(outputDirectory)//"/potfit"
	Call makeDir(filePath) 
!log	
    If(printToTerminal.eq.1.and.mpiProcessID.eq.0)Then
	  print *,ProgramTime(),"Prepare Potfit input files"
	  print *,ProgramTime(),"Saved in "//trim(filePath)
	End If	
!Prepare potential starting file
    Call potfitRunFile(filePath)
!Prepare pot file
	Call potfitEAMFile(filePath)
	
	
  End Subroutine runPrepPotfit
  
  
  Subroutine potfitRunFile(filePath) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
    Character(len=255) :: filePath
    Character(len=255) :: fileName
!If on master process
	If(mpiProcessID.eq.0)Then
!Open file
      fileName = trim(filePath)//"/potfit.run"	
      open(unit=10,file=trim(fileName))
!Write lines
      write(10,"(A)") "# Potfit File created by EAMPA"
	  write(10,"(A)") "ntypes      "//trim(intToString(size(elements,1)))
	  write(10,"(A)") "# Input Files"
	  write(10,"(A)") "config          conf.dat"
	  write(10,"(A)") "startpot        start.pot"
      write(10,"(A)") "tempfile        temp.pot"
	  write(10,"(A)") "flagfile        STOP"
	  write(10,"(A)") "# Final pot"
	  write(10,"(A)") "endpot          opt.pot"
	  write(10,"(A)") "plotfile        TestEAM.plot"
	  write(10,"(A)") "write_pair      1"
	  write(10,"(A)") "# Options"
	  write(10,"(A)") "plotpointfile   TestSampled.plot"
	  write(10,"(A)") "write_lammps    1"
	  write(10,"(A)") "# Output style"
	  write(10,"(A)") "imdpotsteps     1000 "
	  write(10,"(A)") "extend          2"
	  write(10,"(A)") "output_prefix   pf_testEAM"
	  write(10,"(A)") "# Minimization Options"
	  write(10,"(A)") "opt             1 "
	  write(10,"(A)") "anneal_temp     100.0"
	  write(10,"(A)") "eng_weight      100.0 "
	  write(10,"(A)") "stress_weight   20.0"
	  write(10,"(A)") "seed	         123"
	  write(10,"(A)") "d_eps           0.001"
!Close file	
	  close(10) 
	End If  
  End Subroutine potfitRunFile
  
  
  Subroutine potfitEAMFile(filePath) 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
    Character(len=255) :: filePath
    Character(len=255) :: fileName
	
	If(mpiProcessID.eq.0)Then
!Open file
      fileName = trim(filePath)//"/potfitstart.pot"	
      open(unit=10,file=trim(fileName))
	
      write(10,"(A)") "#F 4 "//trim(intToString(size(eamKey,1)))
	
	
	
	
!Close file	
	  close(10)	
	End If
  End Subroutine potfitEAMFile
  
  
  
  
  
  
  
  
  
  Subroutine potentialStartFile() 
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i	
  
  
  
  
  
  End Subroutine potentialStartFile
  
  
    
!------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!  
  
! List of Functions
!-------------------------------------------------------------------------  
! 

  
  
End Module preppotfit  