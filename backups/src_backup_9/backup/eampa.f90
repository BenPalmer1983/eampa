Program eampa

! University of Birmingham
! Ben Palmer
!


! Setup Modules
Use kinds				!data kinds
Use constants			!physical constants module
Use units				!unit conversion and normalisation 
Use strings		        !string functions
Use maths				!maths functions
Use initialise			! initialise program
Use input				! input
Use neighbour			! neighbour list
Use calc			    ! calc

!run initialisation module
Call runInitialise()

!open output file	


!read data from input files
outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
write(999,"(A21,F8.4)") "###Start runInput:   ",ProgramTime()
close(999)
Call runInput()
outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
write(999,"(A19,F8.4)") "###End runInput:   ",ProgramTime()
close(999)


!Build neighbour list
outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
write(999,"(A25,F8.4)") "###Start runNeighbour:   ",ProgramTime()
close(999)
Call runNeighbour()
outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
write(999,"(A23,F8.4)") "###End runNeighbour:   ",ProgramTime()
close(999)

!run calculation
outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
write(999,"(A20,F8.4)") "###Start runCalc:   ",ProgramTime()
close(999)
Call runCalc()
outputFile = trim(currentWorkingDirectory)//"/"//"output.dat"
open(unit=999,file=trim(outputFile),status="old",position="append",action="write")
write(999,"(A18,F8.4)") "###End runCalc:   ",ProgramTime()
close(999)




End