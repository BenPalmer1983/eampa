Module output

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

! ----------------------------------------
! Updated: 12th Aug 2014
! ----------------------------------------

! Setup Modules
  Use kinds
  Use types
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
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines - Output to specific file
  Public :: saveEamFile
  Public :: saveEamNodes
  Public :: outputForcesFile
  Public :: outputAtomEnergiesFile
  Public :: outputNLSeparationFile
  Public :: outputNLFile
  Public :: outputSplineNodes
  Public :: outputInputFiles
  Public :: outputBpData
! Public Subroutines - Output to output file
  Public :: outputNLSummary
  Public :: outputEvaluate
  Public :: outputTimeTaken
  Public :: outputProcessMap
  Public :: outputCpuTimes
  Public :: outputEquilibriumPoints
  Public :: outputZBL
  Public :: outputNLMinMax
  Public :: outputEmbeRescale
  Public :: outputALatTest
  Public :: outputALat
  Public :: outputTestingSummary
  Public :: outputCleanupList
  Public :: outputConfigPoints
  Public :: outputConfigBPPoints
  Public :: outputOptLine
! Public Subroutines - Output to terminal
  Public :: outputConfigSummaryT
  Public :: outputNLSummaryT
  Public :: outputEndT
  Public :: outputNLMinMaxT
  Public :: outputALatT
  Public :: outputTestingSummaryT
  Public :: outputProcessMapT
  Public :: outputPreCalcSummaryT
  Public :: outputPreCalcSummaryBpT
  Public :: outputEAMFunctionsT
  Public :: outputEAMAnParams
  Public :: outputEnergyT
  Public :: outputBpT
  Public :: outputEAMSummaryT
  Public :: outputRssT
  
  
  Contains

! ---------------------------------------------------------------------------------------------------
! Save to specific file
! ---------------------------------------------------------------------------------------------------

  Subroutine saveEamFile(fileName)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(*) :: fileName
    Character(len=255) :: filePath
    Integer(kind=StandardInteger) :: i, j, k, functionCounter
! Only on master process
    If(mpiProcessID.eq.0)Then
      !fileName = Trim(Adjustl(fileName))
      FunctionCounter = 0
      If(fileName(1:1).ne." ")Then
        filePath = Trim(outputDirectory)//"/"//fileName
        Open(UNIT=118,FILE=Trim(filePath))
! Loop through EAM Functions
        Do i=1,size(eamKey,1)
          If(eamKey(i,1).gt.0)Then
            FunctionCounter = functionCounter + 1
            If(eamKey(i,2).gt.0)Then
              write(118,"(A4,A1,A2,A1,A2)") eamFunctionTypes(eamKey(i,3))," ",&
              elements(eamKey(i,1))," ",elements(eamKey(i,2))
            Else
              write(118,"(A4,A1,A2)") eamFunctionTypes(eamKey(i,3))," ",&
              elements(eamKey(i,1))
            End If
            k = 0
            Do j=eamKey(i,4),eamKey(i,6)
              k = k + 1
              write(118,"(E17.10,A1,E17.10,A1,E17.10,A1,E17.10,A1,I5,A1,I5,A1,I5)") &
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
      Close(118)
    End If
  End Subroutine saveEamFile
  
  Subroutine saveEamNodes(fileName)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(*) :: fileName
    Character(len=255) :: filePath
    Integer(kind=StandardInteger) :: i, j, k, functionCounter
! Only on master process
    If(mpiProcessID.eq.0)Then
      !fileName = Trim(Adjustl(fileName))
      FunctionCounter = 0
      If(fileName(1:1).ne." ")Then
        filePath = Trim(outputDirectory)//"/"//fileName
        Open(UNIT=118,FILE=Trim(filePath))
! Loop through EAM Functions
        Do i=1,size(splineNodesKey,1)
          If(splineNodesKey(i,1).gt.0)Then
            FunctionCounter = functionCounter + 1
            If(splineNodesKey(i,2).gt.0)Then
              write(118,"(A4,A1,A2,A1,A2)") eamFunctionTypes(splineNodesKey(i,3))," ",&
              elements(splineNodesKey(i,1))," ",elements(splineNodesKey(i,2))
            Else
              write(118,"(A4,A1,A2)") eamFunctionTypes(splineNodesKey(i,3))," ",&
              elements(splineNodesKey(i,1))
            End If
            k = 0
            Do j=splineNodesKey(i,4),splineNodesKey(i,6)
              k = k + 1
              write(118,"(E17.10,A1,E17.10,A1,E17.10,A1,E17.10,A1,I5,A1,I5,A1,I5)") &
              splineNodesData(j,1)," ",splineNodesData(j,2)," ",splineNodesData(j,3)," ",&
              splineNodesData(j,4)," ",i," ",k," ",j
            End Do
          End If
          If(functionCounter.eq.eamFunctionCount)Then
            Exit
          End If
        End Do
      End If
! Close file
      Close(118)
    End If
! Output to Terminal
    If(TerminalPrint())Then
      Print *,"Saved to: ",trim(filePath)
    End If
  End Subroutine saveEamNodes

  Subroutine outputForcesFile()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, printOut, startKey, endKey
! Only on master process
    If(mpiProcessID.eq.0)Then
      Open(UNIT=1,FILE=Trim(outputDirectory)//"/"//"forces.dat",&
      status="old",position="append",action="write")
      Do i=1,configCount
        write(1,"(A14,I8)") "Configuration ",i
        write(1,"(A40)") "----------------------------------------"
        write(1,"(A20)") "Crystal unit vectors:"
        write(1,"(F14.7,F14.7,F14.7)") crystalUnitCell(i,1), &
        crystalUnitCell(i,2), crystalUnitCell(i,3)
        write(1,"(F14.7,F14.7,F14.7)") crystalUnitCell(i,4), &
        crystalUnitCell(i,5), crystalUnitCell(i,6)
        write(1,"(F14.7,F14.7,F14.7)") crystalUnitCell(i,7), &
        crystalUnitCell(i,8), crystalUnitCell(i,9)
        write(1,"(A20)") "Stress matrix:"
        write(1,"(F14.7,F14.7,F14.7)") configCalcStresses(i,1),&
        configCalcStresses(i,2),configCalcStresses(i,3)
        write(1,"(F14.7,F14.7,F14.7)") configCalcStresses(i,4),&
        configCalcStresses(i,5),configCalcStresses(i,6)
        write(1,"(F14.7,F14.7,F14.7)") configCalcStresses(i,7),&
        configCalcStresses(i,8),configCalcStresses(i,9)
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
        write(1,"(A20)") "Forces:"
        If(printOut.gt.0)Then
          Do j=configurationCoordsKeyG(i,1),configurationCoordsKeyG(i,3)
            If(printOut.eq.1)Then
              write(1,"(I8,A5,A2,A1,E17.10,A1,E17.10,A1,E17.10)") &
              j," Ref ",&
              elements(configurationCoordsIG(j,1))," ",&
              configRefForces(j,1)," ",configRefForces(j,2)," ",configRefForces(j,3)
            End If
            If(printOut.eq.2)Then
              write(1,"(I8,A6,A2,A1,E17.10,A1,E17.10,A1,E17.10)") &
              j," Calc ",&
              elements(configurationCoordsIG(j,1))," ",&
              configCalcForces(j,1)," ",configCalcForces(j,2)," ",configCalcForces(j,3)
            End If
            If(printOut.eq.3)Then
              write(1,"(I8,A5,A2,A1,E17.10,A1,E17.10,A1,E17.10,A5,E17.10,A1,E17.10,A1,E17.10)") &
              j," Ref ",&
              elements(configurationCoordsIG(j,1))," ",&
              configRefForces(j,1)," ",configRefForces(j,2)," ",configRefForces(j,3),"   C ",&
              configCalcForces(j,1)," ",configCalcForces(j,2)," ",configCalcForces(j,3)
            End If
          End Do
        End If
      End Do
      Close(1)
    End If
  End Subroutine outputForcesFile

  Subroutine outputAtomEnergiesFile()
! Output atom energies to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i,j,n
! Only on master process
    If(mpiProcessID.eq.0)Then
      Open(UNIT=173,FILE=Trim(outputDirectory)//"/"//"AtomEnergies1.dat",&
      status="old",position="append",action="write")
      Do i=1,configCount
        write(173,"(A14,I8)") "Configuration ",i
        write(173,"(A40)") "----------------------------------------"
        write(173,"(A20)") "Atom Energies:"
        n = 0
        write(173,"(A12,A16,A16,A16)") &
        "            ","Pair","Embe","Total"
        Do j=configurationCoordsKeyG(i,1),configurationCoordsKeyG(i,3)
          n = n + 1
          write(173,"(I8,A4,F16.8,F16.8,F16.8)") &
          n,elements(configurationCoordsIG(j,1)),configAtomEnergy(j,1),&
          configAtomEnergy(j,2),(configAtomEnergy(j,1)+configAtomEnergy(j,2))
        End Do
      End Do
      write(173,"(A1)") " "
      write(173,"(A1)") " "
      Close(173)
    End If
  End Subroutine outputAtomEnergiesFile

  Subroutine outputNLSeparationFile()
! Output neighbour list to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Only on root process
    If(mpiProcessID.eq.0)Then
! Save tally to file
      Open(UNIT=112,FILE=Trim(outputDirectory)//"/"//"nlSeparation.dat",&
      status="old",position="append",action="write")
      Write(112,"(A8,A1,A6,A1,A8)") "ID      "," ","R/ang "," ","Count   "
      Do i=1,size(atomSeparationSpread,1)
        Write(112,"(I8,A1,F6.3,A1,I8)") i," ",1.0D0*(i/100.0D0)," ",atomSeparationSpread(i)
      End Do
      Close(112)
    End If
  End Subroutine outputNLSeparationFile

  Subroutine outputNLFile()
! Output neighbour list to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, nlStart, nlLength, nlEnd
    Integer(kind=StandardInteger) :: idA, idB, keyA, keyB, keyAB
! Only on root process
    If(mpiProcessID.eq.0)Then
! Save tally to file
      Open(UNIT=114,FILE=Trim(outputDirectory)//"/"//"nlFile.dat")
      Do i=1,maxConfigs
        If(neighbourListKey(i,1).gt.0)Then
          Write(114,"(A20,I8)") "Configuration:      ",i
          nlStart = neighbourListKey(i,1)
          nlLength = neighbourListKey(i,2)
          nlEnd = neighbourListKey(i,3)
          Do j=nlStart,nlEnd
            idA = neighbourListI(j,1)
            idB = neighbourListI(j,2)
            keyA = neighbourListI(j,3)
            keyB = neighbourListI(j,4)
            keyAB = neighbourListI(j,5)
            Write(114,"(I5,A1,A2,A1,A2,A1,I5,A1,I5,A1,I5,A3,F8.4,A3,F8.4,A1,F8.4,A1,F8.4,A3,F8.4,A1,F8.4,A1,F8.4)") &
            j," ",elements(idA)," ",elements(idB)," ",&
            keyA," ",keyB," ",keyAB," | ",&
            neighbourListR(j)," | ",&
            neighbourListCoords(j,7)," ",neighbourListCoords(j,8)," ",neighbourListCoords(j,9)," | ",&
            neighbourListCoords(j,10)," ",neighbourListCoords(j,11)," ",neighbourListCoords(j,12)
          End Do
        End If
      End Do
      Close(114)
    End If
  End Subroutine outputNLFile

  Subroutine outputSplineNodes(fileName, newFileIn)
! Output neighbour list to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, functionCounter
    Integer(kind=StandardInteger), Optional :: newFileIn
    Integer(kind=StandardInteger) :: newFile
    Character(*) :: fileName
! Optional variables
    newFile = 1
    If(Present(newFileIn))Then
      newFile = newFileIn
    End If
! Only on root process
    If(mpiProcessID.eq.0)Then
! Save tally to file
      If(newFile.eq.1)Then
        Open(UNIT=115,FILE=Trim(outputDirectory)//"/"//Trim(fileName))
      Else
        Open(unit=115,file=trim(trim(outputDirectory)//"/"//Trim(fileName)),&
        status="old",position="append",action="write")
      End If
! Loop through EAM Functions
      FunctionCounter = 0
      Do i=1,size(splineNodesKey,1)
        If(splineNodesKey(i,1).gt.0)Then
          FunctionCounter = functionCounter + 1
          write(115,"(A6,I4)") "Node: ",functionCounter
          Do j=splineNodesKey(i,4),splineNodesKey(i,6)
            write(115,"(I4,A1,E17.10,A1,E17.10,A1,E17.10,A1,E17.10,A1,I4,A1,I4)") &
            j," ",splineNodesData(j,1)," ",splineNodesData(j,2)," ",&
            splineNodesData(j,3)," ",splineNodesData(j,4)," ",&
            RoundDP(splineNodesData(j,5))," ",&
            RoundDP(splineNodesData(j,6))
          End Do
          write(1,"(A1)") " "
        End If
        If(functionCounter.eq.eamFunctionCount)Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
      Close(115)
    End If
  End Subroutine outputSplineNodes

  Subroutine outputData(fileName, dataPoints)
! Output neighbour list to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: dataPoints
    Character(*) :: fileName
! Only on root process
    If(mpiProcessID.eq.0)Then
      Open(UNIT=201,FILE=Trim(outputDirectory)//"/"//Trim(fileName))
      Do i=1,size(dataPoints,1)
        If(dataPoints(i,1).lt.-2.0D20)Then
          Exit
        End If
        If(size(dataPoints,1).eq.2)Then
          write(201,"(E17.10,E17.10)") dataPoints(i,1),dataPoints(i,2)
        End If
        If(size(dataPoints,1).eq.3)Then
          write(201,"(E17.10,E17.10,E17.10)") dataPoints(i,1),&
          dataPoints(i,2),dataPoints(i,3)
        End If
        If(size(dataPoints,1).eq.4)Then
          write(201,"(E17.10,E17.10,E17.10,E17.10)") dataPoints(i,1),&
          dataPoints(i,2),dataPoints(i,3),dataPoints(i,4)
        End If
      End Do
      Close(201)
    End If
  End Subroutine outputData

  Subroutine outputBpData(configID,datapoints)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, configID
    Real(kind=DoubleReal), Dimension(:,:) :: datapoints
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"bpData.dat"),&
      status="old",position="append",action="write")
      write(999,"(A32)") "--------------------------------"
      write(999,"(A12,I8)") "Config ID:  ",configID
      write(999,"(A32)") "--------------------------------"
      Do i=1,size(datapoints,1)
        write(999,"(I8,F16.8,F16.8)") &
        i,datapoints(i,1),datapoints(i,2)
      End Do
      write(999,"(A1)") " "
      
      write(999,"(A42)") "     Reference       Calculated        RSS"
      write(999,"(A5,F16.8,A1,F16.8,A1,F16.8)") "aLat ",&
      bpInArr(configID)%alat," ",calcBulkProperties(configID)%alat,&
      " ",(bpInArr(configID)%alat-calcBulkProperties(configID)%alat)**2      
      write(999,"(A5,F16.8,A1,F16.8,A1,F16.8)") "v0   ",&
      bpInArr(configID)%v0," ",calcBulkProperties(configID)%v0,&
      " ",(bpInArr(configID)%v0-calcBulkProperties(configID)%v0)**2      
      write(999,"(A5,F16.8,A1,F16.8,A1,F16.8)") "e0   ",&
      bpInArr(configID)%e0," ",calcBulkProperties(configID)%e0,&
      " ",(bpInArr(configID)%e0-calcBulkProperties(configID)%e0)**2      
      write(999,"(A5,F16.8,A1,F16.8,A1,F16.8)") "b0   ",&
      bpInArr(configID)%b0," ",calcBulkProperties(configID)%b0,&
      " ",(bpInArr(configID)%b0-calcBulkProperties(configID)%b0)**2      
      write(999,"(A5,F16.8,A1,F16.8,A1,F16.8)") "bp0   ",&
      bpInArr(configID)%bp0," ",calcBulkProperties(configID)%bp0,&
      " ",(bpInArr(configID)%bp0-calcBulkProperties(configID)%bp0)**2      
      Close(999)
    End If
    
    
  End Subroutine outputBpData

! ---------------------------------------------------------------------------------------------------
! Save to output file
! ---------------------------------------------------------------------------------------------------

  Subroutine outputInputFiles()
! Output neighbour list summary to output file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
    Character(len=255) :: fileRow
! Only on root process
    If(mpiProcessID.eq.0)Then
! User Input File
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"input.in"),&
      action="write")
      Do i=1,maxFileRows
        fileRow = userInputData(i)
        If(fileRow(1:1).eq." ")Then
          EXIT
        End If
        write(999,"(A)") trim(fileRow)
      End Do
      Close(999)
! EAM File
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"input.pot"),&
      action="write")
      Do i=1,maxFileRows
        fileRow = eamInputData(i)
        If(fileRow(1:1).eq." ")Then
          EXIT
        End If
        write(999,"(A)") trim(fileRow)
      End Do
      Close(999)
! Config File
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"input.config"),&
      action="write")
      Do i=1,maxFileRows
        fileRow = configInputData(i)
        If(fileRow(1:1).eq." ")Then
          EXIT
        End If
        write(999,"(A)") trim(fileRow)
      End Do
      Close(999)
    End If
  End Subroutine outputInputFiles

  Subroutine outputNLSummary()
! Output neighbour list summary to output file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A1)") ""
      write(999,"(A22)") "Neighbour List Summary"
      write(999,"(A32)") "Start   End     Length  Rcutoff "
      Do i=1,maxConfigs
        If(neighbourListKey(i,1).gt.0)Then
          write(999,"(I8,I8,I8,F12.6)") neighbourListKey(i,1),&
          neighbourListKey(i,3),neighbourListKey(i,2),neighbourListKeyR(i,1)
        End If
      End Do
      write(999,"(A1)") ""
      Close(999)
    End If
  End Subroutine outputNLSummary

  Subroutine outputEvaluate()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
  End Subroutine outputEvaluate

  Subroutine outputTimeTaken(textOut,duration)
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
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
      duration = duration
! Temporarily disabled
! open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
! status="old",position="append",action="write")
! write(999,"(A32,F14.6,A1)") textPrint, duration, "s"
! Close(999)
    End If
  End Subroutine outputTimeTaken

  Subroutine outputProcessMap()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(A64)") "                         Process Map                            "
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(A64)") "MapID  Process                                                  "
      Do i=1,maxConfigs
        If(neighbourListKey(i,1).gt.0)Then
          write(999,"(I6,I6,I6)") i,processMap(i)
        End If
      End Do
      write(999,"(A64)") "MapID  Process  (BP)                                            "
      Do i=1,maxConfigsBP
        If(neighbourListKeyBP(i,1).gt.0)Then
          write(999,"(I6,I6,I6)") i,processMapBP(i)
        End If
      End Do
      Close(999)
    End If
  End Subroutine outputProcessMap

  Subroutine outputCpuTimes()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(A64)") "                        CPU Times/s                             "
      write(999,"(A64)") "----------------------------------------------------------------"
      Do i=1,size(cpuTime,1)
        If(cpuTime(i).gt.0.0D0)Then
          write(999,"(I3,A1,A48,F10.5)") i," ",cpuTimeLabels(i),cpuTime(i)
        End If
      End Do
      Close(999)
    End If
  End Subroutine outputCpuTimes

  Subroutine outputEquilibriumPoints(dataPointsL, dataPointsV)
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPointsL
    Real(kind=DoubleReal), Dimension(1:7,1:2) :: dataPointsV
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(A64)") "                 Equilibrium Data Points                        "
      write(999,"(A64)") "----------------------------------------------------------------"
      Do i=1,7
        write(999,"(I3,A1,F16.8,A1,F16.8,A5,F16.8,A1,F16.8)") &
        i," ",dataPointsL(i,1)," ",dataPointsL(i,2),"     ",dataPointsV(i,1)," ",dataPointsV(i,2)
      End Do
      Close(999)
    End If
  End Subroutine outputEquilibriumPoints

  Subroutine outputZBL(zA, zB, pointA, pointB, splineCoeffs)
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: zA, zB
    Real(kind=DoubleReal), Dimension(1:4) :: pointA, pointB
    Real(kind=DoubleReal), Dimension(1:6) :: splineCoeffs
! Only on root process
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A10)") "ZBL Spline"
      write(999,"(A4,I3,A6,I3)") "ZA: ",zA,"  ZB: ",zB
      write(999,"(A9,F14.8,A1,F14.9,A1,F14.8,A1,F14.8)") &
      "Point A: ",pointA(1)," ",pointA(2)," ",pointA(3)," ",pointA(4)
      write(999,"(A9,F14.8,A1,F14.9,A1,F14.8,A1,F14.8)") &
      "Point B: ",pointB(1)," ",pointB(2)," ",pointB(3)," ",pointB(4)
      write(999,"(E10.4,A1,E10.4,A2,E10.4,A4,E10.4,A4,E10.4,A4,E10.4,A3)") &
      splineCoeffs(1),"+",splineCoeffs(2),"x+",splineCoeffs(3),"x^2+",&
      splineCoeffs(4),"x^3+",splineCoeffs(5),"x^4+",splineCoeffs(6),"x^5"
      Close(999)
    End If
  End Subroutine outputZBL

  Subroutine outputNLMinMax(rMin, rMax)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: rMin, rMax
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A39)") "Neighbour List min/max atom separation:"
      write(999,"(A16,F12.6)") "rMin/Angstrom:  ",rMin
      write(999,"(A16,F12.6)") "rMax/Angstrom:  ",rMax
      Close(999)
    End If
  End Subroutine outputNLMinMax

  Subroutine outputEmbeRescale(i,j,rhoMin,rhoMax,embeMin,embeMax)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j
    Real(kind=DoubleReal) :: rhoMin,rhoMax,embeMin,embeMax
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A14)") "Embed Rescale:"
      write(999,"(A21,I8)") "Density Function:    ",i
      write(999,"(A21,I8)") "Embedding Function:  ",j
      write(999,"(A18,F12.6,A2,F12.6)") "rhoMin/rhoMax:    ",rhoMin,"  ",rhoMax
      write(999,"(A18,F12.6,A2,F12.6)") "embeMin/embeMax:  ",embeMin,"  ",embeMax
      Close(999)
    End If
  End Subroutine outputEmbeRescale

! -------- Analysis/Testing Output

  Subroutine outputALatTest(datapoints)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: datapoints
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A32)") "--------------------------------"
      write(999,"(A32)") "Lattice Parameter vs Energy     "
      write(999,"(A32)") "--------------------------------"
      Do i=1,size(datapoints,1)
        If(datapoints(i,1).gt.-2.0D20)Then
          write(999,"(F16.8,F16.8,F16.8,F16.8,F16.8)") &
          datapoints(i,1),datapoints(i,2),datapoints(i,3),datapoints(i,4),datapoints(i,5)
        End If
      End Do
      Close(999)
    End If
  End Subroutine outputALatTest

  Subroutine outputALat(textIn, aLatResult, minEnergyResult, minVolResult, bmResult)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: aLatResult, minEnergyResult, minVolResult, bmResult
    Character(*) :: textIn
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A32)") textIn
      write(999,"(A32,F12.6,A2,F12.6)") "Lattice Parameter/Min Energy:   ", aLatResult, "  ", minEnergyResult
      write(999,"(A32,F12.6,A2,F12.6)") "Opt Vol/Min Energy:             ", minVolResult, "  ", minEnergyResult
      write(999,"(A32,F12.6)")          "Bulk Modulus:                   ", bmResult
      Close(999)
    End If
  End Subroutine outputALat

  Subroutine outputTestingSummary()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: aLatResult, minEnergyResult
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(A64)") "                     EAM Testing Results                        "
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(A64)") "                                                                "
      write(999,"(A64)") "FCC:                                                            "
      write(999,"(A20,F12.6,A3,F12.6,A1)") "Alat:               ",fccALatMurn,"  (",fccALat,")"
      write(999,"(A20,F12.6,A3,F12.6,A1)") "Min Energy:         ",fccEMinMurn,"  (",fccEMin,")"
      write(999,"(A20,F12.6,A3,F12.6,A1)") "Opt Volume:         ",fccVolMinMurn,"  (",fccVolMin,")"
      write(999,"(A20,F12.6,A3,F12.6,A1)") "Bulk Modulus:       ",fccBMMurn,"  (",fccBM,")"
      write(999,"(A20,F12.6)") "C11:                ",fccEC(1)
      write(999,"(A20,F12.6)") "C12:                ",fccEC(2)
      write(999,"(A20,F12.6)") "C44:                ",fccEC(3)
      write(999,"(A64)") "                                                                "
      write(999,"(A64)") "BCC:                                                            "
      write(999,"(A20,F12.6)") "Alat:               ",bccALat
      write(999,"(A20,F12.6)") "Min Energy:         ",bccEMin
      write(999,"(A20,F12.6)") "Opt Volume:         ",bccVolMin
      write(999,"(A20,F12.6,A3,F12.6,A1)") "Bulk Modulus:       ",bccBMMurn,"  (",bccBM,")"
      write(999,"(A20,F12.6)") "C11:                ",bccEC(1)
      write(999,"(A20,F12.6)") "C12:                ",bccEC(2)
      write(999,"(A20,F12.6)") "C44:                ",bccEC(3)
      write(999,"(A64)") "                                                                "
      write(999,"(A64)") "----------------------------------------------------------------"
      write(999,"(F12.6,A2,F12.6)") aLatResult, "  ", minEnergyResult
      Close(999)
    End If
  End Subroutine outputTestingSummary

  Subroutine outputCleanupList()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
    Character(len=512) :: testLine
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      write(999,"(A32)") "Cleaning temporary files:       "
      Do i=1,100
        testLine = fileCleanupList(i)
        If(testLine(1:1).eq." ")Then
          Exit
        Else
          write(999,"(A)") adjustl(trim(testLine))
        End If
      End Do
      Close(999)
    End If
  End Subroutine outputCleanupList

  Subroutine outputConfigPoints(printTypeIn)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! In/Out    
    Integer(kind=StandardInteger), Optional :: printTypeIn
! Private variables
    Integer(kind=StandardInteger) :: printType
    Integer(kind=StandardInteger) :: configID,i
    Integer(kind=StandardInteger) :: coordStart, coordLength, coordEnd
! Optional
    printType = 0
    If(Present(printTypeIn))Then
      printType = printTypeIn
    End If
! Root process    
    If(mpiProcessID.eq.0)Then
      If(printType.eq.0)Then
        open(unit=999,file=trim(trim(outputDirectory)//"/"//"configPoints.dat"))
        write(999,"(A22,I8)") "Configurations:       ", configCount
        Do configID=1,configCount
          coordStart = configurationCoordsKeyG(configID,1)
          coordLength = configurationCoordsKeyG(configID,2)
          coordEnd = configurationCoordsKeyG(configID,3)
          write(999,"(A7,I4,A3,I8,I8,A1,I8,A1)") "Config ",configID,"   ",coordStart,coordEnd,"(",coordLength,")"
          Do i=coordStart,coordEnd
            write(999,"(A8,A2,F14.6,A2,F14.6,A2,F14.6,A4,F14.6,A2,F14.6,A2,F14.6)") &
            adjustl(elements(configurationCoordsIG(i,1))),"  ",&
            configurationCoordsRG(i,1),"  ",&
            configurationCoordsRG(i,2),"  ",&
            configurationCoordsRG(i,3),"    ",&
            configurationCoordsRG(i,4),"  ",&
            configurationCoordsRG(i,5),"  ",&
            configurationCoordsRG(i,6)
          End Do
        End Do  
        Close(999)
      End If
      If(printType.eq.1)Then
        open(unit=999,file=trim(trim(outputDirectory)//"/"//"configPointsAngs.dat"))
        write(999,"(A22,I8)") "Configurations:       ", configCount
        Do configID=1,configCount
          coordStart = configurationCoordsKeyG(configID,1)
          coordLength = configurationCoordsKeyG(configID,2)
          coordEnd = configurationCoordsKeyG(configID,3)
          write(999,"(A7,I4,A3,I8,I8,A1,I8,A1)") "Config ",configID,"   ",coordStart,coordEnd,"(",coordLength,")"
          Do i=coordStart,coordEnd
            write(999,"(A8,A2,F14.6,A2,F14.6,A2,F14.6)") &
            adjustl(elements(configurationCoordsIG(i,1))),"  ",&
            configurationCoordsRG(i,1),"  ",&
            configurationCoordsRG(i,2),"  ",&
            configurationCoordsRG(i,3)
          End Do
        End Do
        Close(999)
      End If
      If(printType.eq.2)Then
        open(unit=999,file=trim(trim(outputDirectory)//"/"//"configPointsFrac.dat"))
        write(999,"(A22,I8)") "Configurations:       ", configCount
        Do configID=1,configCount
          coordStart = configurationCoordsKeyG(configID,1)
          coordLength = configurationCoordsKeyG(configID,2)
          coordEnd = configurationCoordsKeyG(configID,3)
          write(999,"(A7,I4,A3,I8,I8,A1,I8,A1)") "Config ",configID,"   ",coordStart,coordEnd,"(",coordLength,")"
          Do i=coordStart,coordEnd
            write(999,"(A8,A2,F14.6,A2,F14.6,A2,F14.6)") &
            adjustl(elements(configurationCoordsIG(i,1))),"  ",&
            configurationCoordsRG(i,4),"  ",&
            configurationCoordsRG(i,5),"  ",&
            configurationCoordsRG(i,6)
          End Do
        End Do
        Close(999)
      End If
    End If
  End Subroutine outputConfigPoints

  Subroutine outputConfigBPPoints()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID,i
    Integer(kind=StandardInteger) :: coordStart, coordLength, coordEnd
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"configPointsBP.dat"))
      write(999,"(A22,I8)") "Configurations:       ", configCountBP
      Do configID=1,configCountBP
        coordStart = configurationCoordsKeyBP(configID,1)
        coordLength = configurationCoordsKeyBP(configID,2)
        coordEnd = configurationCoordsKeyBP(configID,3)
        write(999,"(A7,I4,A3,I8,I8,A1,I8,A1)") "Config ",configID,"   ",coordStart,coordEnd,"(",coordLength,")"
        Do i=coordStart,coordEnd
          write(999,"(A8,A4,A2,F14.6,A2,F14.6,A2,F14.6,A4,F14.6,A2,F14.6,A2,F14.6)") &
          "        ",&
          elements(configurationCoordsIBP(i,1)),"  ",&
          configurationCoordsRBP(i,1),"  ",&
          configurationCoordsRBP(i,2),"  ",&
          configurationCoordsRBP(i,3)
        End Do
      End Do
      Close(999)
    End If
  End Subroutine outputConfigBPPoints
  
  
  
  Subroutine outputOptLine()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    If(mpiProcessID.eq.0)Then
      optLogCounter = optLogCounter + 1
      open(unit=999,file=trim(outputDirectory)//"/"//"OptLog.dat",&
      status="old",position="append",action="write")
      write(999,"(I8,A20)") optLogCounter,"Input Configurations"
      Do configID=1,configCount
        write(999,"(I4,A1,F16.8,A1,F16.8)")&
        configID," ",&
        configCalcEnergies(configID)," ",rssConfigsArr(configID)%total
      End Do
      write(999,"(I8,A20)") optLogCounter,"BP Configurations"
      Do configID=1,configCountBP
       write(999,"(I4,A1,F16.8,A1,F16.8,A1,F16.8,A1,F16.8)") &
       configID," ",&
       calcBulkProperties(configID)%alat," ",calcBulkProperties(configID)%e0," ",&
       calcBulkProperties(configID)%b0," ",rssBPArr(configID)%total
      End Do   
      write(999,"(A12,F16.8)") "Total RSS:  ",totalRSS
      write(999,"(A1)") " "
      Close(999)
    End If  
  End Subroutine outputOptLine
! ---------------------------------------------------------------------------------------------------
! Output to terminal
! ---------------------------------------------------------------------------------------------------

  Subroutine outputProcessMapT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Only on root process
    If(TerminalPrint())Then
      print *,""
      print *,"                             Process Map                                 "
      Call printBR() ! -----------------------------------------------------------
      print *,"Config  Process                                                          "
      Do i=1,maxConfigs
        If(processMap(i).ge.0)Then
          print "(I6,A2,I6)",i,"  ",processMap(i)
        Else
          exit
        End If
      End Do
      print *,"Config  Process  (BP)                                           "
      Do i=1,maxConfigsBP
        If(processMapBP(i).ge.0)Then
          print "(I6,A2,I6)",i,"  ",processMapBP(i)
        Else
          exit
        End If
      End Do      
      print *,""
    End If
  End Subroutine outputProcessMapT

  Subroutine outputPreCalcSummaryT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
! Only on root process
    If(TerminalPrint())Then
      print *,""
      print *,"                           Config Summary                                "
      Call printBR() ! -----------------------------------------------------------
      print *,"Conf    Atoms    Vol        Rcut   Vcut   Rmin   Rmax   Nl Len  ConfW       "
      Do configID=1,configCount
        print "(I4,A3,I6,A3,F10.4,A1,F6.3,A1,F6.3,A1,F6.3,A1,F6.3,A1,I6,A2,F6.3)",&
        configID,"   ",&
        configurationCoordsKeyG(configID,2),"   ",&
        configVolume(configID)," ",&
        neighbourListKeyR(configID,1)," ",&
        neighbourListKeyR(configID,6)," ",&
        neighbourListKeyR(configID,2)," ",&
        neighbourListKeyR(configID,3)," ",&
        neighbourListKey(configID,2),"  ",&
        configWeighting(configID)        
      End Do
      print *,"Total atoms: ",configsAtomTotal
      print *,""
    End If
  End Subroutine outputPreCalcSummaryT
  
  
  Subroutine outputEAMFunctionsT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, functionCounter
! Only on root process
    If(TerminalPrint())Then
      print *,"                    Read EAM Potential Functions                         "
      Call printBR() ! -----------------------------------------------------------
      If(potentialType.eq.1)Then 
        print *,"Spline/Tabulated Potential"
      End If
      If(potentialType.eq.2)Then 
        print *,"Analytic Potential ",apfData(1)%functionCount 
      End If
 
      print *,"EAM File: ",trim(eamFilePath)
      print *,"EAM Type: ",eamType
      If(eamType.eq.1)Then
        print *,"Pair: ",eamPairCount,", Dens: ",eamDensCount,&
        ", Dens: ",eamEmbeCount,", Total: ",eamFunctionCount
      End If
      print *,"EAM Potential Functions:"
      functionCounter = 0
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
          If(eamKey(i,2).gt.0)Then
            print "(I4,A1,A4,A1,I2,A2,A2,A1,I2,A2,A2,A1,I2,A2,I7,A1,I7,A1,I7)",&
            functionCounter," ",eamFunctionTypes(eamKey(i,3)),"(",eamKey(i,3),") ",&
            elements(eamKey(i,1)),"(",eamKey(i,1),") ",&
            elements(eamKey(i,2)),"(",eamKey(i,2),") ",&
            eamKey(i,4)," ",eamKey(i,5)," ",eamKey(i,6)
          Else
            print "(I4,A1,A4,A1,I2,A2,A2,A1,I2,A2,A7,I7,A1,I7,A1,I7)",&
            functionCounter," ",eamFunctionTypes(eamKey(i,3)),"(",eamKey(i,3),") ",&
            elements(eamKey(i,1)),"(",eamKey(i,1),") ",&
            "       ",&
            eamKey(i,4)," ",eamKey(i,5)," ",eamKey(i,6)
          End If
        End If
        If(functionCounter.eq.eamFunctionCount )Then
          Exit  ! Exit, all functions cycled through
        End If
      End Do
      print *,""
    End If
  End Subroutine outputEAMFunctionsT  
  
  Subroutine outputEAMAnParams(apfData_Fin)
! Output params
    Implicit None   ! Force declaration of all variables
! In/Out    
    Type(analyticFunctions), Dimension(:) :: apfData_Fin
! Private variables
    Integer(kind=StandardInteger) :: funcID, paramID

    Do funcID=1,apfData_Fin(1)%functionCount
      Print *, "Function ID: ",funcID
      Print *, "Parameters: ",apfData_Fin(funcID)%functionParameterCount
      Do paramID=1,apfData_Fin(funcID)%functionParameterCount
        Print *, "   ",paramID,apfData_Fin(funcID)%functionParameters(paramID)
      End Do
    End Do


  End Subroutine outputEAMAnParams  
  
  
  Subroutine outputPreCalcSummaryBpT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
! Only on root process
    If(TerminalPrint())Then
      print *,""
      print *,"                      Bulk Property Config Summary                       "
      Call printBR() ! -----------------------------------------------------------
      print *,"Conf    Atoms    Vol        Rcut   Vcut   Rmin   Rmax   Nl Len          "
      Do configID=1,configCountBP
        print "(I4,A3,I6,A3,F10.4,A1,F6.3,A1,F6.3,A1,F6.3,A1,F6.3,A1,I6)",&
        configID,"   ",&
        configurationCoordsKeyBP(configID,2),"   ",&
        0.0D0," ",&
        neighbourListKeyRBP(configID,1)," ",&
        neighbourListKeyRBP(configID,6)," ",&
        neighbourListKeyRBP(configID,2)," ",&
        neighbourListKeyRBP(configID,3)," ",&
        neighbourListKeyBP(configID,2)
      End Do
      print *,"Total atoms: ",configsAtomTotalBP
      print *,""
    End If
  End Subroutine outputPreCalcSummaryBpT

  Subroutine outputEnergyT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
! Only on root process
    If(TerminalPrint())Then
      print *,""
      print *,"                         Config Energies                                 "
      Call printBR() ! -----------------------------------------------------------
      print *,"Conf   Atoms    Energy      Energy per Atom                              "
      Do configID=1,configCount
        print "(I4,A3,I8,A1,F14.7,A1,F14.7)",&
        configID,"   ",configurationCoordsKeyG(configID,2)," ",&
        (configCalcEnergies(configID)*1.0D0*&
        configurationCoordsKeyG(configID,2))," ",configCalcEnergies(configID)
      End Do
      print *,""
    End If
  End Subroutine outputEnergyT
  
  
  Subroutine outputBpT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
! Only on root process
    If(TerminalPrint())Then
      Print *,""
      Print *,"Configurations:"
      print "(A7,A9,A15,A15,A15,A15)",&
      " Config",&
      " Atoms",&
      "  Calc Total E ",&
      "  Calc EPA ",&
      "  Ref EPA ",&
      "  RSS "
      
      Do configID=1,configCount
        print "(I4,A3,I8,A1,ES14.5E3,A1,ES14.5E3,A2,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        configID,"   ",&
        configurationCoordsKeyG(configID,2)," ",&
        (configCalcEnergies(configID)*1.0D0*configurationCoordsKeyG(configID,2)),&
        " ",configCalcEnergies(configID),&
        " (",configRefEnergies(configID),")",&
        " [",rssConfigsArr(configID)%energy,"]"
      End Do
      print *,"Bulk Protperties:"
      Do configID=1,configCountBP
        print *,"Config ",configID      
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "aLat: ",calcBulkProperties(configID)%aLat,&
        " (ref: ",bpInArr(configID)%alat,")"," [",rssBPArr(configID)%aLat,"]"
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "v0: ",calcBulkProperties(configID)%v0,&
        " (ref: ",bpInArr(configID)%v0,")"," [",rssBPArr(configID)%v0,"]"
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "e0: ",calcBulkProperties(configID)%e0,&
        " (ref: ",bpInArr(configID)%e0,")"," [",rssBPArr(configID)%e0,"]"
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "b0: ",calcBulkProperties(configID)%b0,&
        " (ref: ",bpInArr(configID)%b0,")"," [",rssBPArr(configID)%b0,"]"    
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "C11: ",calcBulkProperties(configID)%c11,&
        " (ref: ",bpInArr(configID)%c11,")"," [",rssBPArr(configID)%c11,"]"     
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "C12: ",calcBulkProperties(configID)%c12,&
        " (ref: ",bpInArr(configID)%c12,")"," [",rssBPArr(configID)%c12,"]"    
        print "(A6,ES14.5E3,A7,ES14.5E3,A1,A2,ES14.5E3,A1)",&
        "C44: ",calcBulkProperties(configID)%c44,&
        " (ref: ",bpInArr(configID)%c44,")"," [",rssBPArr(configID)%c44,"]"      
        print "(A9,ES14.5E3)","EoS rss: ",rssBPArr(configID)%eos
      End Do
      print "(A12,ES14.5E3)", "Config RSS: ",rssConfigsArrTotal%total
      print "(A12,ES14.5E3)", "BP RSS:     ",rssBPArrTotal%total
      print "(A12,ES14.5E3)", "Total RSS:  ",totalRSS
    End If
  End Subroutine outputBpT 
  
  
  Subroutine outputRssT()
! Output ref and calc forces to file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal) :: configRSS, bpConfigRSS    
! Init values
    configRSS = 0.0D0  
    bpConfigRSS = 0.0D0     
! Config RSS
    print *,"Config RSS"
    print *,"Config       Energy ref    Energy calc   Energy rss    ",&
    "Forces rss     Stress rss     "
    Do i=1,configCount
      print "(A1,I3,A6,F14.7,F14.7,F14.7,F14.7,F14.7)",&
      " ",i,"     ",&
      configRefEnergies(i),configCalcEnergies(i),&
      rssConfigsArr(i)%energy,rssConfigsArr(i)%force,&
      rssConfigsArr(i)%stress
! sum RSS
      configRSS = configRSS + rssConfigsArr(i)%energy + &
      rssConfigsArr(i)%force + rssConfigsArr(i)%stress      
    End Do
    print *,"Config rss: ",configRSS
   
    print *,"BP RSS"    
    print *,"    Alat           V0             E0             B0             EoS   "
    Do i=1,configCountBP
      print "(F14.7,A1,F14.7,A1,F14.7,A1,F14.7,A1,F14.7)",&
      rssBPArr(i)%alat," ",rssBPArr(i)%v0," ",rssBPArr(i)%e0," ",&
      rssBPArr(i)%b0," ",rssBPArr(i)%eos
      bpConfigRSS = bpConfigRSS + rssBPArr(i)%v0 + rssBPArr(i)%alat + rssBPArr(i)%e0 +&
      rssBPArr(i)%b0 + rssBPArr(i)%eos
    End Do
    print *,"BP Config rss: ",bpConfigRSS
    print *,"Total rss: ",(configRSS+bpConfigRSS)," (",totalRSS,")"
  
  End Subroutine outputRssT 
  
  
  Subroutine outputEAMSummaryT()
! Output summary of the current EAM after calculations have been run
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, configID, functionCounter
    Type(tableObj) :: table
    Character(Len=16), Dimension(1:1) :: rowA_Str
    Real(kind=DoubleReal), Dimension(1:1) :: rowA
    Real(kind=DoubleReal), Dimension(1:3) :: row
    Real(kind=DoubleReal) :: rssVal
    Character(Len=512) :: rowHeadings
    Character(Len=4) :: configStr
! Only on root process
    If(TerminalPrint())Then
      Call printBR() ! ----------------
      Print *,"                         Summary of current EAM Potential                              "
      Call printBR() ! ----------------
      Print *,""
      Print *,"EAM Potential details:"
      Print *,"Type: ",eamType
      Do i=1,size(eamKey,1)
        If(eamKey(i,1).gt.0)Then
          functionCounter = functionCounter + 1
          If(eamKey(i,2).gt.0)Then
            Print *,i,eamFunctionTypes(eamKey(i,3))," ",&
            elements(eamKey(i,1))," ",elements(eamKey(i,2)),eamKey(i,4)," to ",eamKey(i,6)
          Else
            Print *,i,eamFunctionTypes(eamKey(i,3)),"    ",&
            elements(eamKey(i,1)),eamKey(i,4)," to ",eamKey(i,6)
          End If
        End If
        If(functionCounter.eq.eamFunctionCount)Then
          Exit
        End If
      End Do       
      
      Print *,""
      print *,"                                Bulk Protperties                                    "
      Call printBR() ! -----------------------------------------------------------
      Do configID=1,configCountBP
        Call printTableInit(table)
        table%colWidth = 20        
        table%padR = 1
        table%padL = 1
        configStr = IntToStr(configID)
        Call printTableAddHeadersRC(table,trim("Config "//adjustl(configStr)))
        rowHeadings = BlankString(rowHeadings)
        rowHeadings = "Structure,Atoms in Calc,NL Length,"
        rowHeadings = trim(rowHeadings)//"Alat (ang),v0 (ang3),e0 (eV),b0 (ev/ang3),c11 (ev/ang3),"
        rowHeadings = trim(rowHeadings)//"c12 (ev/ang3),c44 (ev/ang3),b0 (GPa),"
        rowHeadings = trim(rowHeadings)//"c11 (GPa),c12 (GPa),c44 (GPa)"
        Call printTableAddHeadersR(table,"Calculated,Reference,RSS")    
        Call printTableAddHeadersC(table,trim(rowHeadings))  
        
        ! atom structure 
        rowA_Str(1) = Trim(bpInArr(configID)%structure)
        Call printTableAddRow(table,rowA_Str)    
        ! atoms in calc 
        rowA(1) = configurationCoordsKeyBP(configID,2)
        Call printTableAddRow(table,rowA)    
        ! neighbour list length 
        rowA(1) = neighbourListKeyBP(configID,2)
        Call printTableAddRow(table,rowA)       
        ! alat
        row(1) = calcBulkProperties(configID)%aLat
        row(2) = bpInArr(configID)%alat
        row(3) = rssBPArr(configID)%aLat
        Call printTableAddRow(table,row)
        ! v0
        row(1) = calcBulkProperties(configID)%v0
        row(2) = bpInArr(configID)%v0
        row(3) = rssBPArr(configID)%v0
        Call printTableAddRow(table,row)    
        ! e0
        row(1) = calcBulkProperties(configID)%e0
        row(2) = bpInArr(configID)%e0
        row(3) = rssBPArr(configID)%e0
        Call printTableAddRow(table,row)   
        ! b0
        row(1) = calcBulkProperties(configID)%b0
        row(2) = bpInArr(configID)%b0
        row(3) = rssBPArr(configID)%b0
        Call printTableAddRow(table,row)   
        ! c11
        row(1) = calcBulkProperties(configID)%c11
        row(2) = bpInArr(configID)%c11
        row(3) = rssBPArr(configID)%c11
        Call printTableAddRow(table,row)   
        ! c12
        row(1) = calcBulkProperties(configID)%c12
        row(2) = bpInArr(configID)%c12
        row(3) = rssBPArr(configID)%c12
        Call printTableAddRow(table,row)  
        ! c44
        row(1) = calcBulkProperties(configID)%c44
        row(2) = bpInArr(configID)%c44
        row(3) = rssBPArr(configID)%c44
        Call printTableAddRow(table,row)  
        ! b0
        row(1) = UnitConvert(calcBulkProperties(configID)%b0,"EVAN3","GPA")
        row(2) = UnitConvert(bpInArr(configID)%b0,"EVAN3","GPA")
        rssVal = RSSCalc(row(1),row(2),rssWeighting(6))
        row(3) = rssVal
        Call printTableAddRow(table,row)  
        ! c11
        row(1) = UnitConvert(calcBulkProperties(configID)%c11,"EVAN3","GPA")
        row(2) = UnitConvert(bpInArr(configID)%c11,"EVAN3","GPA")
        rssVal = RSSCalc(row(1),row(2),rssWeighting(7))
        row(3) = rssVal
        Call printTableAddRow(table,row)   
        ! c12
        row(1) = UnitConvert(calcBulkProperties(configID)%c12,"EVAN3","GPA")
        row(2) = UnitConvert(bpInArr(configID)%c12,"EVAN3","GPA")
        rssVal = RSSCalc(row(1),row(2),rssWeighting(7))
        row(3) = rssVal
        Call printTableAddRow(table,row)  
        ! c44
        row(1) = UnitConvert(calcBulkProperties(configID)%c44,"EVAN3","GPA")
        row(2) = UnitConvert(bpInArr(configID)%c44,"EVAN3","GPA")
        rssVal = RSSCalc(row(1),row(2),rssWeighting(7))
        row(3) = rssVal
        Call printTableAddRow(table,row)  
    
        Call printTableMake(table)       
        
      End Do
    End If
  End Subroutine outputEAMSummaryT 
  
  
  
  
  
  
  
  

  Subroutine outputConfigSummaryT()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Print out
    If(TerminalPrint())Then
      Print *, "Configuration Summary"
      Print *, "Count: ",configCount,configCountT
      Do i=1,1024
        If(configurationCoordsKeyG(i,1).gt.0)Then
          Print "(A9,I4,A1,I8,A1,I8,A1,I8,A1,F8.5,A1,I4,A1,I4,A1,I4,A1,F12.6)",&
          "  Config ",i," ",configurationCoordsKeyG(i,1),&
          " ",configurationCoordsKeyG(i,3)," ",configurationCoordsKeyG(i,2),&
          " ",configurationsR(i,1)," ",configurationsI(i,1)," ",configurationsI(i,2),&
          " ",configurationsI(i,3)," ",configVolume(i)
        End If
      End Do
    End If
  End Subroutine outputConfigSummaryT

  Subroutine outputNLSummaryT()
! Output neighbour list summary to output file
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i
! Only on root process
    If(TerminalPrint())Then
      print *,"Neighbour List Summary"
      print *,"ID    Start   End     Length  Rcutoff "
      Do i=1,1024
        If(neighbourListKey(i,1).gt.0)Then
          print "(I5,A1,I8,I8,I8,F12.6)",i," ",neighbourListKey(i,1),&
          neighbourListKey(i,3),neighbourListKey(i,2),neighbourListKeyR(i,1)
        End If
      End Do
    End If
  End Subroutine outputNLSummaryT

  Subroutine outputEndT()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
! Print out
    If(TerminalPrint())Then
      Print *,"----------------------------------------------------------------"
      Print *, "The program is about to terminate."
      Print *, "Program time: ",ProgramTime()
      Print *,"----------------------------------------------------------------"
    End If
  End Subroutine outputEndT

  Subroutine outputNLMinMaxT(rMin, rMax)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: rMin, rMax
! Print out
    If(TerminalPrint())Then
      Print *,"Neighbour List min/max atom separation:"
      Print *,"rMin/Angstrom:  ",rMin
      Print *,"rMax/Angstrom:  ",rMax
    End If
  End Subroutine outputNLMinMaxT

  Subroutine outputALatT(textIn, aLatResult, minEnergyResult, minVolResult, bmResult)
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: aLatResult, minEnergyResult, minVolResult, bmResult
    Character(*) :: textIn
! Print out
    If(TerminalPrint())Then
      print *," "
      Print *, textIn
      Print *, "Lattice Parameter(ang)/Min Energy(eV):   ", aLatResult, "  ", minEnergyResult
      Print *, "Opt Vol (ang3)/Min Energy(eV):           ", minVolResult, "  ", minEnergyResult
      Print *, "Bulk Modulus(GPa):                       ", bmResult
      print *," "
    End If
  End Subroutine outputALatT

  Subroutine outputTestingSummaryT()
! Saves the eam file to the output directory
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: aLatResult, minEnergyResult
! Print out
    If(mpiProcessID.eq.0)Then
      open(unit=999,file=trim(trim(outputDirectory)//"/"//"output.dat"),&
      status="old",position="append",action="write")
      print "(A64)","----------------------------------------------------------------"
      print "(A64)","                     EAM Testing Results                        "
      print "(A64)","----------------------------------------------------------------"
      print "(A64)","                                                                "
      print "(A64)","FCC:                                                            "
      print "(A20,F12.6,A3,F12.6,A1)","Alat:               ",fccALatMurn,"  (",fccALat,")"
      print "(A20,F12.6,A3,F12.6,A1)","Min Energy:         ",fccEMinMurn,"  (",fccEMin,")"
      print "(A20,F12.6,A3,F12.6,A1)","Opt Volume:         ",fccVolMinMurn,"  (",fccVolMin,")"
      print "(A20,F12.6,A3,F12.6,A1)","Bulk Modulus:       ",fccBMMurn,"  (",fccBM,")"
      print "(A20,F12.6)","C11:                ",fccEC(1)
      print "(A20,F12.6)","C12:                ",fccEC(2)
      print "(A20,F12.6)","C44:                ",fccEC(3)
      print "(A64)","                                                                "
      print "(A64)","BCC:                                                            "
      print "(A20,F12.6)","Alat:               ",bccALat
      print "(A20,F12.6)","Min Energy:         ",bccEMin
      print "(A20,F12.6)","Opt Volume:         ",bccVolMin
      print "(A20,F12.6,A3,F12.6,A1)","Bulk Modulus:       ",bccBMMurn,"  (",bccBM,")"
      print "(A20,F12.6)","C11:                ",bccEC(1)
      print "(A20,F12.6)","C12:                ",bccEC(2)
      print "(A20,F12.6)","C44:                ",bccEC(3)
      print "(A64)","                                                                "
      print "(A20,F12.6)","Testing RSS:        ",testingRSS
      print "(A64)","                                                                "
      print "(A64)","----------------------------------------------------------------"
      write(999,"(F12.6,A2,F12.6)") aLatResult, "  ", minEnergyResult
      Close(999)
    End If
  End Subroutine outputTestingSummaryT

End Module output
