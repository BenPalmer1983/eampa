Module output

! Setup Modules
  Use kinds
  Use constants
  Use strings		!string functions
  Use maths
  Use initialise
  Use input
  Use neighbour
  Use calc


!force declaration of all variables
  Implicit None
  
!declare global variables  

!Privacy of functions/subroutines/variables
  Private
  Public :: runOutput		        !Subroutine
  
  
!------------------------------------------------------------------------!
!                                                                        !
! MODULE SUBROUTINES                                                     !
!                                                                        !
!                                                                        !
!------------------------------------------------------------------------!
  
contains 

!Run all the input subroutines

  Subroutine runOutput()
	
!Internal subroutine variables
	Integer(kind=StandardInteger) :: i, j, k
!open output file	
	outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
	open(unit=999,file=trim(outputFile),status="old",position="append",action="write")		
!write to output file
    write(999,"(F8.4)") ProgramTime()
    write(999,"(A1)") " "
    write(999,"(A1)") " "
    write(999,"(A48)") "================================================"	
    write(999,"(A48)") "|            Calculation Summary               |"
    write(999,"(A48)") "================================================"	
    write(999,"(A1)") " "
	write(999,"(A10,F8.4)") "Run time: ",ProgramTime()
    write(999,"(A1)") " "
	do i=1,configCount
!loop through configurations
	  write(999,"(A7,I4,A1)") "Config ",i,":"
	  
!if energy run
	  if(calcRunType(1:6).eq."ENERGY")then
	    write(999,"(A21,I8,A1,I8)") "Atoms:          ",configAtomsUnitCell(i),&
		"/",configAtoms(i) 
	    write(999,"(A21,F20.14,A8)") "Volume:              ", &
		configurationVolume(i)," Ang^3  "
	    write(999,"(A21,F20.14,A8)") "Energy:              ", &
		configurationEnergy(i)," eV   "
	  endif
	  
!if bulk modulus run
	  if(calcRunType(1:11).eq."BULKMODULUS")then	    
		write(999,"(A21,I8,A1,I8)") "Atoms:          ",configAtomsUnitCell(i),&
		"/",configAtoms(i) 
	    write(999,"(A21,F20.14,A8)") "Opt volume:          ", &
		configurationOptVolume(i)," Ang^3  "
	    write(999,"(A21,F20.14,A8)") "Opt energy:          ", &
		configurationOptEnergy(i)," eV     "
	    write(999,"(A21,F20.14,A8)") "Bulk modulus:        ", &
		configurationBM(i)," GPa    "
	  endif
	  
	  
	  
      write(999,"(A1)") " "
	enddo
	
    write(999,"(A48)") "================================================"	
	
!close output file
	close(999)

  End Subroutine runOutput
  
  
    
  

End Module output