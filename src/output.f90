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
  Public :: outputForces
  Public :: outputBM
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
		If(configurationStresses(i,1).gt.-2.0D20)Then
		  write(999,"(I4,A22)") i,&
		  "  Calculated Stresses:"
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Sxx:",configurationStresses(i,1)," Sxy:",&
		  configurationStresses(i,2)," Sxz:",configurationStresses(i,3)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Syx:",configurationStresses(i,4)," Syy:",&
		  configurationStresses(i,5)," Syz:",configurationStresses(i,6)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Szx:",configurationStresses(i,7)," Szy:",&
		  configurationStresses(i,8)," Szz:",configurationStresses(i,9)
		End If
		If(configurationRefStresses(i,1).gt.-2.0D20)Then
		  write(999,"(I4,A21)") i,&
		  "  Reference Stresses:"
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Sxx:",configurationRefStresses(i,1)," Sxy:",&
		  configurationRefStresses(i,2)," Sxz:",configurationRefStresses(i,3)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Syx:",configurationRefStresses(i,4)," Syy:",&
		  configurationRefStresses(i,5)," Syz:",configurationRefStresses(i,6)
		  write(999,"(I4,A6,F18.10,A5,F18.10,A5,F18.10)") i,&
		  "  Szx:",configurationRefStresses(i,7)," Szy:",&
		  configurationRefStresses(i,8)," Szz:",configurationRefStresses(i,9)
		End If
!Calculated equilibrium volume
		If(configurationEquVolume(i).gt.-2.0D20)Then
		    write(999,"(I4,A22,F18.10)") i," Input volume: ",&
			configurationVolume(i)
		    write(999,"(I4,A22,F18.10)") i," Equilibrium volume: ",&
			configurationEquVolume(i)
		    write(999,"(I4,A22,F18.10)") i," Ref Equilibrium volume: ",&
			configurationRefEquVolume(i)
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
		  write(999,"(I4,A20,F18.10)") i,&
		  "  Energy RSS (1):   ",configurationRSS(i,1)
		End If
		If(configurationRSS(i,2).gt.0.0D0)Then
		  write(999,"(I4,A20,F18.10)") i,&
		  "  Stress RSS (2):   ",configurationRSS(i,2)
		End If
		If(configurationRSS(i,3).gt.0.0D0)Then
		  write(999,"(I4,A20,F18.10)") i,&
		  "  Force RSS (3):    ",configurationRSS(i,3)
		End If
		If(configurationRSS(i,4).gt.0.0D0)Then
		  write(999,"(I4,A20,F18.10)") i,&
		  "  Eq Vol RSS (4):   ",configurationRSS(i,4)
		End If
		If(configurationRSS(i,5).gt.0.0D0)Then
		  write(999,"(I4,A20,F18.10)") i,&
		  "  BM RSS (5):       ",configurationRSS(i,5)
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

  Subroutine outputForces() 
!Outputs configuration forces to file
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n
	Integer(kind=StandardInteger) :: configAtomsCount
!If root/master process
    If(mpiProcessID.eq.0)Then
!-------------------------------------
! Output Forces to force file
!-------------------------------------
!Write forces to output file, if master process
!open forces output file	
	  open(unit=5,file=trim(trim(outputDirectory)//"/forces.dat"),action="write")
	  !open(unit=989,file=trim(outputFileForces),status="old",position="append",action="write")	   
!Loop through configurations
	  Do i=1,configCount	    
	    configAtomsCount = configAtoms(i)
		write(5,"(A15,I4)") "Configuration: ",i
	    Do j=1,configAtomsCount
	      n = configAtomsStart(i) + j		  
		  If(configurationRefForce(n,1).gt.-2.0D20)Then
		    write(5,"(I4,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10)") &
			j,&
			ForceZero(configurationForce(n,1),1.0D-15),&
			ForceZero(configurationForce(n,2),1.0D-15),&
			ForceZero(configurationForce(n,3),1.0D-15),&
			ForceZero(configurationRefForce(n,1),1.0D-15),&
			ForceZero(configurationRefForce(n,2),1.0D-15),&
			ForceZero(configurationRefForce(n,3),1.0D-15)
		  Else	
			write(5,"(I4,E20.10,E20.10,E20.10)") &
			j,&
			ForceZero(configurationForce(n,1),1.0D-15),&
			ForceZero(configurationForce(n,2),1.0D-15),&
			ForceZero(configurationForce(n,3),1.0D-15)
		  End If
	    End Do
	  End Do	 
	  close(5)	
!----------------------------------
  	End If  
  End Subroutine outputForces
  
  
  Subroutine outputBM(configurationID,bmData)
!Outputs configuration forces to file
!force declaration of all variables
	Implicit None	
!declare private variables
	Integer(kind=StandardInteger) :: i,j,k,n,configurationID
	Integer(kind=StandardInteger) :: configAtomsCount
	Real(kind=DoubleReal), Dimension(1:50,1:3) :: bmData
!If root/master process
    If(mpiProcessID.eq.0)Then  
!----------------------------------
!open forces output file	
	  open(unit=5,file=trim(trim(outputDirectory)//"/bm.dat"),action="write")	
	  write(5,"(A15,I4)") "Configuration: ",configurationID
	  Do i=1,size(bmData,1)
	    If(bmData(i,1).gt.-2.0D20)Then
		  write(5,"(I8,E20.10,E20.10,E20.10)") i,bmData(i,1),bmData(i,2),bmData(i,3)
		End If
	  End Do
	  close(5)
!----------------------------------
  	End If    
  End Subroutine outputBM
  
  

  
  

End Module output