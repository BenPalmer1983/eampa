Module output

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use prep


!force declaration of all variables
  Implicit None
!Include MPI header
  Include 'mpif.h'  
  
!declare global variables  

!Privacy of functions/subroutines/variables
  Private
!Variables  
!Subroutines  
  Public :: calcOutput
  Public :: outputPrepareEAM	    
!Functions
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!------------------------------------------------------------------------!
!                                                                        !
! OUTPUT TO TERMINAL                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!







!------------------------------------------------------------------------!
!                                                                        !
! OUTPUT TO OUTPUT.DAT                                                   !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
!------------------------------------------------------------------------!
! ***calcOutput*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!  
  Subroutine calcOutput () 
!Outputs configuration calculation details
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k 
	Integer(kind=StandardInteger) :: configCounter, atomCounter, ecComponentCounter
!Write energies to output file, if master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	 
!Write to file	  
	  write(999,"(A44)") "--------------------------------------------"
	  write(999,"(A19,I4)") "Calculation count: ",globalCounter(1)
	  write(999,"(A44)") "--------------------------------------------"
	  Do i=1,configCount
!Configuration Energy
	    write(999,"(I4,A1,F12.6,A4,I8,A10,I8)") i," ",&
	      (1.0*configurationEnergy(i))," eV ",&
	    configAtoms(i)," atoms     ",configAtomsMap(i,3)
		If(configAtomsMap(i,5).gt.0)Then
	      write(999,"(I4,F12.6,A14,F12.6,A8)") i,&
		  (configurationEnergy(i)/configAtoms(i))," eV per atom (",&
		  configurationRefEnergy(i)," eV Ref)"
		Else
		  write(999,"(I4,F12.6,A21)") i,&
		  (configurationEnergy(i)/configAtoms(i))," eV per atom (no ref)"
		End If
!Calculated stress
		If(calcStressOnOff.eq.1)Then
		  write(999,"(I4,A22)") i,&
		  "  Calculated Stresses:"
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Sxx:",configurationStress((i-1)*9+1)," Sxy:",&
		  configurationStress((i-1)*9+2)," Sxz:",configurationStress((i-1)*9+3)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Syx:",configurationStress((i-1)*9+4)," Syy:",&
		  configurationStress((i-1)*9+5)," Syz:",configurationStress((i-1)*9+6)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Szx:",configurationStress((i-1)*9+7)," Szy:",&
		  configurationStress((i-1)*9+8)," Szz:",configurationStress((i-1)*9+9)
		End If
		If(configurationRefStress((i-1)*9+1).gt.-2.0D20)Then
		  write(999,"(I4,A21)") i,&
		  "  Reference Stresses:"
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Sxx:",configurationRefStress((i-1)*9+1)," Sxy:",&
		  configurationRefStress((i-1)*9+2)," Sxz:",configurationRefStress((i-1)*9+3)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Syx:",configurationRefStress((i-1)*9+4)," Syy:",&
		  configurationRefStress((i-1)*9+5)," Syz:",configurationRefStress((i-1)*9+6)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Szx:",configurationRefStress((i-1)*9+7)," Szy:",&
		  configurationRefStress((i-1)*9+8)," Szz:",configurationRefStress((i-1)*9+9)
		End If
!Calculated equilibrium volume
		If(configurationEquVolume(i).gt.-2.0D20)Then
		    write(999,"(I4,A22,F18.10,A15,F18.10)") i,&
		    "  Equilibrium volume: ",configurationEquVolume(i),&
			" Input volume: ",configurationVolume(i)
		End If
!Calculated equilibrium lattice parameter multiplier
		If(configurationEquVolume(i).gt.-2.0D20)Then
		    write(999,"(I4,A44,F18.10)") i,&
		    "  Equilibrium lattice parameter multiplier: ",configurationEquLat(i)
		End If
!Calculated bulk modulus
		If(configurationBM(i).gt.-2.0D20)Then
		  If(configurationRefBM(i).gt.-2.0D20)Then
		    write(999,"(I4,A20,F18.10,A9,F18.10)") i,&
		    "  Bulk Modulus/GPa: ",configurationBM(i)," Ref BM: ",&
		    configurationRefBM(i)
		  Else
		    write(999,"(I4,A20,F18.10,A13)") i,&
		    "  Bulk Modulus/GPa: ",configurationBM(i)," (No Ref BM) "
		  End If
		End If
!Calculated bulk modulus
		If(configurationEC(i,1).gt.-2.0D20)Then
		  ecComponentCounter = 0
		  Do j=1,6
		    If(configurationEC(i,j).gt.-2.0D20)Then
			  ecComponentCounter = ecComponentCounter + 1
			End If
		  End Do
		  If(ecComponentCounter.eq.3)Then		!Cubic
		    write(999,"(I4,A37)") i,&
		    "  Cubic elastic constants C11 C12 C44"
			If(configurationRefEC(i,1).gt.-2.0D20)Then
			  write(999,"(I4,A23,F18.10,A2,F18.10,A2,F18.10)") i,&
		      "  C11,C12,C44/GPa:     ",&
			  configurationEC(i,1),", ",configurationEC(i,2),", ",configurationEC(i,3)
			  write(999,"(I4,A23,F18.10,A2,F18.10,A2,F18.10)") i,&
			  "  Ref C11,C12,C44/GPa: ",&
		      configurationRefEC(i,1),", ",configurationRefEC(i,2),", ",configurationRefEC(i,3)
			Else
			  write(999,"(I4,A23,F18.10,A2,F18.10,A2,F18.10)") i,&
		      "  C11,C12,C44/GPa:     ",&
			  configurationEC(i,1),", ",configurationEC(i,2),", ",configurationEC(i,3)
			End If
		  End If	
		End If		
!RSS between calculated and 
		If(configurationRSS(i,1).gt.0.0D0)Then
		  write(999,"(I4,A14,F18.10)") i,&
		  "  Energy RSS: ",configurationRSS(i,1)
		End If
		If(configurationRSS(i,2).gt.0.0D0)Then
		  write(999,"(I4,A13,F18.10)") i,&
		  "  Stress RSS: ",configurationRSS(i,2)
		End If
		If(configurationRSS(i,3).gt.0.0D0)Then
		  write(999,"(I4,A13,F18.10)") i,&
		  "  Force RSS: ",configurationRSS(i,3)
		End If
		If(configurationRSS(i,4).gt.0.0D0)Then
		  write(999,"(I4,A13,F18.10)") i,&
		  "  Eq Vol RSS: ",configurationRSS(i,4)
		End If
	    write(999,"(A1)") " "
	  End Do 
	  If(trialResidualSquareSum.gt.-2.0D20)Then
		write(999,"(A24,E20.10)") "RSS all configurations: ",trialResidualSquareSum
	  End If
	  !If(eamFunctionWobbliness(1).gt.0.0D0)Then	    
		!write(999,"(A24)") "EAM Functions Wobbliness: "
		!Do i=1,size(eamFunctionWobbliness,1)
		!  write(999,"(A4,I4,E20.10)") "Pot ",i,eamFunctionWobbliness(i)
		!End Do
	  !End If
	  write(999,"(A1)") " "
!Close output file
      close(999)
	End If  
  End Subroutine calcOutput

  
  
!------------------------------------------------------------------------!
! ***outputPrepareEAM*** 
! Calculates the energy (and forces, if required) of all the configurations
! splitting the work over the MPI processes
!------------------------------------------------------------------------!  
  Subroutine outputPrepareEAM () 
!Outputs configuration calculation details
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i, j, k 
	Integer(kind=StandardInteger) :: configCounter, atomCounter, ecComponentCounter
!Write energies to output file, if master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!open output file	
	  outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	  open(unit=999,file=trim(outputFile),status="old",position="append",action="write")	 
!Write to file	  
	  write(999,"(A59)") "EAM potential functions read in and output in EAMPA format."
	  write(999,"(A64)") trim(eamPreparedFile)
!Close output file
      close(999)
	End If  
  End Subroutine outputPrepareEAM

!------------------------------------------------------------------------!
!                                                                        !
! OUTPUT TO OTHER FILE                                                   !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!



  
  

End Module output