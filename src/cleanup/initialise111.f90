

    
    
! Loop through DFT files
    Open(UNIT=1,FILE=Trim(configFilePathTC))
    Do i=1,maxFileRows
! Read in line
      Read(1,"(A255)",IOSTAT=ios) fileRow
! Break out If end of file
      If(ios /= 0)Then
        EXIT
      End If
      If(fileRow(1:7).eq."#NEWDFT")Then  ! Clean variables
        dftFilePath = BlankString(dftFilePath)
        dftType = BlankString(dftType)
        dftReplaceLabel = BlankString2DArray(dftReplaceLabel)
        eqVol = -2.1D20
      End If
      If(fileRow(1:5).eq."#PATH")Then  ! Path to DFT file
        Read(fileRow,*) bufferA, dftFilePath
      End If
      If(fileRow(1:5).eq."#TYPE")Then  ! Ab init file type to read in
        Read(fileRow,*) bufferA, dftType
      End If
      If(fileRow(1:3).eq."#EV")Then  ! Equilibrium volume
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) eqVol
        eqVol = UnitConvert(eqVol, bufferC, "ANG3")
      End If
      If(fileRow(1:3).eq."#CW")Then  ! Config Weighting
        Read(fileRow,*) bufferA, bufferB
        Read(bufferB,*) confWeight
      End If
      If(fileRow(1:3).eq."#RC")Then  ! Radius cutoff
        Read(fileRow,*) bufferA, bufferB, bufferC
        Read(bufferB,*) radiusCutoff
        radiusCutoff = UnitConvert(radiusCutoff, bufferC, "ANGS")
      End If
      If(fileRow(1:5).eq."#REPL")Then  ! Replace label
        Read(fileRow,*) bufferA, bufferB, bufferC
        Do j=1,size(dftReplaceLabel,1)
          If(dftReplaceLabel(j,1).eq."        ")Then
            dftReplaceLabel(j,1) = trim(bufferB)
            dftReplaceLabel(j,2) = trim(bufferC)
            Exit
          End If
        End Do
      End If
      If(fileRow(1:7).eq."#ENDDFT")Then  ! Clean variables
        If(trim(dftType).eq."PWSCF")Then
          Call readPWSCFFile(dftFilePath, configFilePathT, &
          eqVol, radiusCutoff, confWeight)
        End If
      End If
    End Do
    Close(1)
! Remove unnecessary files
    Call rmFile(configFilePathTA)
    Call rmFile(configFilePathTB)
    Call rmFile(configFilePathTC)
    Call rmFile(configFilePathTDFT)
! Add files to clean
! Call fileToClean(configFilePathT)
! Call fileToClean(configFilePathTA)
! Call fileToClean(configFilePathTB)
! Call fileToClean(configFilePathTC)
! Call fileToClean(configFilePathTDFT)
! Call rmFile(configFilePathT)
! Call rmFile(configFilePathTA)
! Call rmFile(configFilePathTB)
! Call rmFile(configFilePathTC)
! Call rmFile(configFilePathTDFT)