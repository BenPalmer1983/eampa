Module output

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
  Use globals
  Use initialise
  Use loadData  
! Force declaration of all variables
  Implicit None
!Privacy of variables/functions/subroutines
  Private    
!Public Subroutines
  Public :: saveEamFile
  Public :: outputForcesFile
  Public :: outputEvaluate
  Public :: outputTimeTaken
  
Contains


 
  Subroutine saveEamFile(fileName)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(len=64) :: fileName
    Character(len=255) :: filePath
    Integer(kind=StandardInteger) :: i, j, k, functionCounter    
    fileName = Trim(Adjustl(fileName))
    functionCounter = 0
    If(fileName(1:1).ne." ")Then
      filePath = Trim(outputDirectory)//"/"//Trim(fileName)    
      Open(UNIT=1,FILE=Trim(filePath)) 
! Loop through EAM Functions
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
          If(eamKey(i,2).gt.0)Then
            write(1,"(A4,A1,A2,A1,A2)") eamFunctionTypes(eamKey(i,3))," ",&
            elements(eamKey(i,1))," ",elements(eamKey(i,2))
          Else
            write(1,"(A4,A1,A2)") eamFunctionTypes(eamKey(i,3))," ",&
            elements(eamKey(i,1))       
          End If        
          k = 0
          Do j=eamKey(i,4),eamKey(i,6)
            k = k + 1
            write(1,"(E17.10,A1,E17.10,A1,E17.10,A1,E17.10,A1,I5,A1,I5,A1,I5)") &
            eamData(j,1)," ",eamData(j,2)," ",eamData(j,3)," ",&
            eamData(j,4)," ",i," ",k," ",j
          End Do
        End If
        If(functionCounter.eq.eamFunctionCount)Then
          Exit
        End If  
      End Do
    End If
! Close file
    Close(1)
  End Subroutine saveEamFile 
    

  Subroutine outputForcesFile()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, printOut, startKey, endKey   
! Only on master process
    If(mpiProcessID.eq.0)Then
    Open(UNIT=1,FILE=Trim(outputDirectory)//"/forces.dat") 
    Do i=1,configCount
      write(1,"(A15,I8)") "Configuration: ",i
      printOut = 0
      startKey = configurationCoordsKeyG(i,1)
      endKey = configurationCoordsKeyG(i,3)
      If(configRefForces(startKey,1).gt.-2.0D20)Then
        printOut = 1
      End If
      If(configCalcForces(startKey,1).gt.-2.0D20)Then
        printOut = 2
      End If
      If(configRefForces(startKey,1).gt.-2.0D20.and.configCalcForces(startKey,1).gt.-2.0D20)Then
        printOut = 3
      End If      
      If(printOut.gt.0)Then
        Do j=configurationCoordsKeyG(i,1),configurationCoordsKeyG(i,3)
          If(printOut.eq.1)Then
            write(1,"(I8,A3,E17.10,A1,E17.10,A1,E17.10)") &
            j," R ",configRefForces(j,1)," ",configRefForces(j,2)," ",configRefForces(j,3)
          End If
          If(printOut.eq.2)Then
            write(1,"(I8,A3,E17.10,A1,E17.10,A1,E17.10)") &
            j," C ",configCalcForces(j,1)," ",configCalcForces(j,2)," ",configCalcForces(j,3)
          End If
          If(printOut.eq.3)Then
             write(1,"(I8,A3,E17.10,A1,E17.10,A1,E17.10,A5,E17.10,A1,E17.10,A1,E17.10)") &
            j," R ",configRefForces(j,1)," ",configRefForces(j,2)," ",configRefForces(j,3),"   C ",&
            configCalcForces(j,1)," ",configCalcForces(j,2)," ",configCalcForces(j,3)
          End If
        End Do
      End If
    End Do
    Close(1)
    End If    
  End Subroutine outputForcesFile 
  
  
  

  
  
!----------------------------------------
! Save to output file
!----------------------------------------  
  
  Subroutine outputEvaluate()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, totalAtoms
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A1)") " "
      write(999,"(F8.4,A2,A26)") ProgramTime(),"  ","Configuration Evaluations:"
      write(999,"(A5,A5,A7,A13,A13,A13)") "Cfg  ","Proc ","Atoms  ","Ref Energy   ",&
      "Calc Energy  ","RSS          "
      totalAtoms = 0
      Do configID=1,configCount
        write(999,"(I4,A1,I4,A1,I6,A1,F12.4,A1,F12.4,A1,F12.4,A1)") &
        configID," ",processMap(configID,1)," ",&
        configurationCoordsKeyG(configID,2)," ",&
        (configRefEnergies(configID)*configurationCoordsKeyG(configID,2))," ",&
        configCalcEnergies(configID)," ",&
        configRSS(configID,size(configRSS,2))," "    
        totalAtoms = totalAtoms + configurationCoordsKeyG(configID,2)       
      End Do
      write(999,"(A36,I8)")  "Total atoms:                        ",totalAtoms
      write(999,"(A36,F12.4)") "Total RSS all configurations:       ",totalRSS
      write(999,"(A1)") " "
  
  
      Close(999)
    End If  
  End Subroutine outputEvaluate 
  
  
  Subroutine outputTimeTaken(textOut,duration)
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: duration
    Character(*) :: textOut
    Character(Len=48) :: textPrint    
! Only on root process
    If(mpiProcessID.eq.0)Then
      Do i=1,48
        textPrint(i:i) = "."
      End Do  
      Do i=1,Len(textOut)
        textPrint(i:i) = textOut(i:i)
      End Do    
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A32,F14.6,A1)") textPrint, duration, "s"
      Close(999)
    End If  
  End Subroutine outputTimeTaken 
  
  
  
  
  
  
  
  
  
  
End Module output  