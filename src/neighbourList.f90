Module neighbourList

! --------------------------------------------------------------!
! General subroutines and functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! Read user input file

! ----------------------------------------
! Updated: 16th Sept 2014
! ----------------------------------------

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
  Use output
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: makeNeighbourList
  Public :: clearNeighbourList
  Public :: saveConfigNL
  Public :: loadConfigNL
  Public :: isotropicDistortion
  Public :: tetragonalDistortion
  Public :: orthogonalDistortion
  Public :: monoclinicDistortion
  Public :: applyDistortionNL

  Contains
  Subroutine makeNeighbourList(configIDStartIn,configIDEndIn,forcePrintIn)
! Make neighbour list of all atom pairs for all configs where separation is .le. rcutoff
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, coordStart, coordLength, coordEnd
    Integer(kind=StandardInteger) :: atomA, atomB, nlKey, neighbourListCount, asKey, lastNLC
    Integer(kind=StandardInteger) :: configStart, configLength
    Integer(kind=StandardInteger) :: xCopy, yCopy, zCopy
    Integer(kind=StandardInteger) :: i, l, m, n
    Real(kind=DoubleReal) :: rCutoffSq
    Real(kind=DoubleReal) :: aLat, xShift, yShift, zShift
    Real(kind=DoubleReal) :: xA, xB, yA, yB, zA, zB, xdSq, ydSq, zdSq, rdSq
    Real(kind=DoubleReal) :: rMin, rMax
    Integer(kind=StandardInteger), optional :: configIDStartIn, configIDEndIn, forcePrintIn
    Integer(kind=StandardInteger) :: configIDStart, configIDEnd, forcePrint
    configIDStart = 1                   ! Starting configID
    If(present(configIDStartIn))Then
      configIDStart = configIDStartIn
    End If
    configIDEnd = 1024                   ! Ending configID
    If(present(configIDEndIn))Then
      configIDEnd = configIDEndIn
    End If
    forcePrint = 0                   ! Ending configID
    If(present(forcePrintIn))Then
      forcePrint = forcePrintIn
    End If
! Start time
    Call cpu_time(timeStart)
! Init variables
    neighbourListCount = 0
    configStart = 1
    rMin = 2.0D21
    rMax = -2.0D21
! Prepare NL key array - clear all above configIDStart
    Do i=configIDStart,1024
      neighbourListKey(i,1) = 0
      neighbourListKey(i,2) = 0
      neighbourListKey(i,3) = 0
    End Do
    If(configIDStart.gt.1)Then
      lastNLC = 0
      Do i=1,configIDStart
        If(lastNLC.lt.neighbourListKey(i,3))Then
          lastNLC = neighbourListKey(i,3)
        End If
      End Do
      neighbourListCount = lastNLC
      configStart = lastNLC+1
    End If
! print *,"NL",configStart
! Loop through configurations
    Do configID=configIDStart,configIDEnd
! Check config is there
      If(configurationCoordsKeyG(configID,1).gt.0)Then
! Init looping variables
        atomA = 0
        atomB = 0
        nlKey = 0
        configLength = 0
        coordStart = configurationCoordsKeyG(configID,1)
        coordLength = configurationCoordsKeyG(configID,2)
        coordEnd = configurationCoordsKeyG(configID,3)
        If(nlCutoff.lt.0.0D0)Then
          rCutoffSq = configurationsR(configID,11)**2
        Else
          rCutoffSq = nlCutoff**2
        End If
        xCopy = configurationsI(configID,1)
        yCopy = configurationsI(configID,2)
        zCopy = configurationsI(configID,3)
        aLat = configurationsR(configID,1)
! loop through Atom B 3x3x3
        Do l=-1,1
          Do m=-1,1
            Do n=-1,1
! Set co-ordinate shift
              xShift = aLat * xCopy * l
              yShift = aLat * yCopy * m
              zShift = aLat * zCopy * n
! Reset unique key list
              nlUniqueKeys = 0
! Loop through atom pairs
              Do atomA=1,coordLength
                Do atomB=1,coordLength
                  If(l.eq.0.and.m.eq.0.and.n.eq.0.and.atomA.eq.atomB)Then  ! Don't self count atom
                  Else
                    If(atomA.lt.atomB)Then
                      nlKey = (atomB-1)*(atomB-2)/2+atomA
                    Else
                      nlKey = (atomA-1)*(atomA-2)/2+atomB
                    End If
                    If(nlUniqueKeys(nlKey).eq.0)Then
                      nlUniqueKeys(nlKey) = 1
! check range one co-ord at a time then distance squared
                      xA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,1)
                      xB = 1.0D0*(xshift + configurationCoordsRG(coordStart+atomB-1,1))
                      yA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,2)
                      yB = 1.0D0*(yshift + configurationCoordsRG(coordStart+atomB-1,2))
                      zA = 1.0D0*configurationCoordsRG(coordStart+atomA-1,3)
                      zB = 1.0D0*(zshift + configurationCoordsRG(coordStart+atomB-1,3))
                      xdSq = (xA-xB)**2
                      If(xdSq.le.rCutoffSq)Then
                        ydSq = (yA-yB)**2
                        If(ydSq.le.rCutoffSq)Then
                          zdSq = (zA-zB)**2
                          If(zdSq.le.rCutoffSq)Then
                            rdSq = xdSq + ydSq + zdSq
                            If(rdSq.le.rCutoffSq)Then
                              neighbourListCount = neighbourListCount + 1
                              configLength = configLength + 1
! Store atom type/id data
                              neighbourListI(neighbourListCount,1) = &
                              configurationCoordsIG(coordStart+atomA-1,1)  !Atom A type
                              neighbourListI(neighbourListCount,2) = &
                              configurationCoordsIG(coordStart+atomB-1,1)  !Atom B type
                              neighbourListI(neighbourListCount,3) = atomA  !Atom A id
                              neighbourListI(neighbourListCount,4) = atomB  !Atom B id
                              neighbourListI(neighbourListCount,5) = nlKey  !Atom A-B Key
                              If(l.eq.0.and.m.eq.0.and.n.eq.0)Then
                                neighbourListI(neighbourListCount,6) = 1
                              Else
                                neighbourListI(neighbourListCount,6) = 0
                              End If
! Store atom separation
                              neighbourListR(neighbourListCount) = rdSq**0.5
! Atom coordinate data
                              neighbourListCoords(neighbourListCount,1) = xA
                              neighbourListCoords(neighbourListCount,2) = yA
                              neighbourListCoords(neighbourListCount,3) = zA
                              neighbourListCoords(neighbourListCount,4) = xB
                              neighbourListCoords(neighbourListCount,5) = yB
                              neighbourListCoords(neighbourListCount,6) = zB
                              neighbourListCoords(neighbourListCount,7) = xB-xA
                              neighbourListCoords(neighbourListCount,8) = yB-yA
                              neighbourListCoords(neighbourListCount,9) = zB-zA
                              neighbourListCoords(neighbourListCount,10) = &
                              (xB-xA)/neighbourListR(neighbourListCount)
                              neighbourListCoords(neighbourListCount,11) = &
                              (yB-yA)/neighbourListR(neighbourListCount)
                              neighbourListCoords(neighbourListCount,12) = &
                              (zB-zA)/neighbourListR(neighbourListCount)
! Tally atom separation
                              asKey = Ceiling(neighbourListR(neighbourListCount)*100)
                              If(asKey.lt.1)Then
                                asKey = 1
                              End If
                              atomSeparationSpread(asKey) = atomSeparationSpread(asKey) + 1
                              If(neighbourListR(neighbourListCount).lt.rMin)Then
                                rMin = neighbourListR(neighbourListCount)
                              End If
                              If(neighbourListR(neighbourListCount).gt.rMax)Then
                                rMax = neighbourListR(neighbourListCount)
                              End If
                            End If
                          End If
                        End If
                      End If
                    End If
                  End If
                End Do
              End Do
            End Do
          End Do
        End Do
! Store nl key
        neighbourListKey(configID,1) = configStart
        neighbourListKey(configID,2) = configLength
        neighbourListKey(configID,3) = configStart+configLength-1
        neighbourListKeyR(configID,1) = rCutoffSq**0.5D0
        configStart = configStart + configLength
      End If
    End Do  ! End loop configs
    If(configIDStart.eq.1.or.forcePrint.eq.1)Then
! Save summary to output file
      Call outputNLSummary()
      Call outputNLSummaryT()
! End If
! If(configIDStart.eq.1)Then
! Output NL separation to file
      Call outputNLSeparationFile()
! Output entire neighbour list to file
      If(saveNLToFile.eq.1)Then
        Call outputNLFile()
      End If
! Output min max atom separation
      Call outputNLMinMaxT(rMin,rMax)  ! To terminal
      Call outputNLMinMax(rMin,rMax)  ! To output file
    End If
! End time
    Call cpu_time(timeEnd)
! Record time taken to make neighbour list
    Call outputTimeTaken("Neighbour List",timeEnd-timeStart)
    Call storeTime(10,timeEnd-timeStart)
  End Subroutine makeNeighbourList
! -------------------------------------------
  Subroutine clearNeighbourList()
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: i, j, k
    If(configCount.gt.0)Then
! Clear data
      Do i=1,configCount
        Do j=neighbourListKey(i,1),neighbourListKey(i,3)
          If(neighbourListKey(i,1).gt.0.and.neighbourListKey(i,3).gt.0)Then
            neighbourListR(j) = 0.0D0
            Do k=1,6
              neighbourListI(j,k) = 0
            End Do
            Do k=1,12
              neighbourListCoords(j,k) = 0.0D0
            End Do
          End If
        End Do
      End Do
! Clear keys
      Do i=1,configCount
        Do j=1,3
          neighbourListKey(i,j) = 0
        End Do
      End Do
! Reset counter
      neighbourListCount = 0
    End If
  End Subroutine clearNeighbourList
! -------------------------------------------
  Subroutine storeCNL()
    Implicit None   ! Force declaration of all variables
! Private variables
    neighbourListCountInput = neighbourListCount
    neighbourListKeyInput = neighbourListKey
    neighbourListKeyRInput = neighbourListKeyR
    neighbourListIInput = neighbourListI
    neighbourListRInput = neighbourListR
    neighbourListCoordsInput = neighbourListCoords
    configCountInput = configCount
    configurationsIInput = configurationsI
    configurationsRInput = configurationsR
    coordCountInput = coordCount
    configurationCoordsKeyInput = configurationCoordsKey
    configurationCoordsIInput = configurationCoordsI
    configurationCoordsRInput = configurationCoordsR
    configurationForcesRInput = configurationForcesR
    coordCountGInput = coordCountG
    configurationCoordsKeyGInput = configurationCoordsKeyG
    configurationCoordsIGInput = configurationCoordsIG
    configurationCoordsRGInput = configurationCoordsRG
    configVolumeInput = configVolume
    configRefInput = configRef
    configRefForcesInput = configRefForces
    configRefStressesInput = configRefStresses
    configRefEnergiesInput = configRefEnergies
    configRefEVInput = configRefEV
    configRefBMInput = configRefBM
  End Subroutine storeCNL

! ------------------------------------------------
! Neighbour List - Save & Load
! ------------------------------------------------

  Subroutine saveConfigNL(configID)
! Save NL for one config to the temp memory array
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, keyS, keyE, i, j
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
! copy data
    Do j=1,size(neighbourListI,2)
      Do i=keyS,keyE
        neighbourListIT(i,j) = neighbourListI(i,j)
      End Do
    End Do
    Do i=keyS,keyE
      neighbourListRT(i) = neighbourListR(i)
    End Do
    Do j=1,size(neighbourListCoords,2)
      Do i=keyS,keyE
        neighbourListCoordsT(i,j) = neighbourListCoords(i,j)
      End Do
    End Do
  End Subroutine saveConfigNL

  Subroutine loadConfigNL(configID)
! Load NL for one config from the temp memory array
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, keyS, keyE, i, j
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
! copy data
    Do j=1,size(neighbourListI,2)
      Do i=keyS,keyE
        neighbourListI(i,j) = neighbourListIT(i,j)
      End Do
    End Do
    Do i=keyS,keyE
      neighbourListR(i) = neighbourListRT(i)
    End Do
    Do j=1,size(neighbourListCoords,2)
      Do i=keyS,keyE
        neighbourListCoords(i,j) = neighbourListCoordsT(i,j)
      End Do
    End Do
  End Subroutine loadConfigNL

! ------------------------------------------------
! Distortion + NL subroutines
! ------------------------------------------------

  Subroutine isotropicDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn
    End If
! set distortion matrix
    distortion(1,1) = factor*(1.0D0+distortionAmount)
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = factor*(1.0D0+distortionAmount)
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*(1.0D0+distortionAmount)
    Call applyDistortionNL(configID, distortion, 1)
  End Subroutine isotropicDistortion

  Subroutine tetragonalDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn
    End If
! set distortion matrix
    distortion(1,1) = factor*(1.0D0+distortionAmount)
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = factor*(1.0D0+distortionAmount)
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*((1.0D0+distortionAmount)**(-2.0D0))
    Call applyDistortionNL(configID, distortion, 1)
  End Subroutine tetragonalDistortion

  Subroutine orthogonalDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn
    End If
! set distortion matrix
    distortion(1,1) = factor*(1.0D0+distortionAmount)
    distortion(1,2) = 0.0D0
    distortion(1,3) = 0.0D0
    distortion(2,1) = 0.0D0
    distortion(2,2) = factor*(1.0D0-distortionAmount)
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*(1.0D0+(distortionAmount**2)/(1.0D0-distortionAmount**(2.0D0)))
    Call applyDistortionNL(configID, distortion, 1)
  End Subroutine orthogonalDistortion

  Subroutine monoclinicDistortion(configID, distortionAmount, factorIn)
! Apply an isotropic distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: distortionAmount, factor
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Real(kind=DoubleReal), Optional :: factorIn
! Optional input
    factor = 1.0D0
    If(Present(factorIn))Then
      factor = factorIn
    End If
! set distortion matrix
    distortion(1,1) = factor*1.0D0
    distortion(1,2) = factor*0.5D0*distortionAmount
    distortion(1,3) = 0.0D0
    distortion(2,1) = factor*0.5D0*distortionAmount
    distortion(2,2) = factor*1.0D0
    distortion(2,3) = 0.0D0
    distortion(3,1) = 0.0D0
    distortion(3,2) = 0.0D0
    distortion(3,3) = factor*(1.0D0+(distortionAmount**2)/(4.0D0-distortionAmount**(2.0D0)))
    Call applyDistortionNL(configID, distortion, 1)
  End Subroutine monoclinicDistortion

  Subroutine applyDistortionNL(configID, distortion, typeCalc)
! Apply a distortion to the neighbour list
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, typeCalc, distOption
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: distortion
    Integer(kind=StandardInteger) :: i, j, arrayCheck, keyS, keyE
    Real(kind=DoubleReal) :: xCA, yCA, zCA, xCB, yCB, zCB
! Init variables
    keyS = neighbourListKey(configID,1)
    keyE = neighbourListKey(configID,3)
    distOption = 2 ! Full nl modification for distortion
! Check type of distortion (just energy, isotropic distortion)
    If(typeCalc.eq.1)Then  ! Energy only calculation
! Check distortion array
      If(distortion(1,1).eq.distortion(2,2).and.distortion(1,1).eq.distortion(3,3))Then
        arrayCheck = 0
        Do i=1,3
          Do j=1,3
            If(j.ne.i)Then
              If(distortion(i,j).ne.0.0D0)Then
                arrayCheck = 1
              End If
            End If
          End Do
        End Do
        If(arrayCheck.eq.0)Then
          distOption = 1
        End If
      End If
    End If
! Distort neighbour list
    If(distOption.eq.1)Then
! Isotropic distortion, energy calculation only - just change nl separation
      Do i=keyS,keyE
        neighbourListR(i) = distortion(1,1)*neighbourListR(i)
      End Do
    Else
! Non isotropic/energy-force-stress calculations
      Do i=keyS,keyE
! Apply distortion vector - Point A
        xCA = neighbourListCoords(i,1)*distortion(1,1)+&
        neighbourListCoords(i,2)*distortion(1,2)+&
        neighbourListCoords(i,3)*distortion(1,3)
        yCA = neighbourListCoords(i,1)*distortion(2,1)+&
        neighbourListCoords(i,2)*distortion(2,2)+&
        neighbourListCoords(i,3)*distortion(2,3)
        zCA = neighbourListCoords(i,1)*distortion(3,1)+&
        neighbourListCoords(i,2)*distortion(3,2)+&
        neighbourListCoords(i,3)*distortion(3,3)
        neighbourListCoords(i,1) = xCA
        neighbourListCoords(i,2) = yCA
        neighbourListCoords(i,3) = zCA
! Apply distortion vector - Point A
        xCB = neighbourListCoords(i,4)*distortion(1,1)+&
        neighbourListCoords(i,5)*distortion(1,2)+&
        neighbourListCoords(i,6)*distortion(1,3)
        yCB = neighbourListCoords(i,4)*distortion(2,1)+&
        neighbourListCoords(i,5)*distortion(2,2)+&
        neighbourListCoords(i,6)*distortion(2,3)
        zCB = neighbourListCoords(i,4)*distortion(3,1)+&
        neighbourListCoords(i,5)*distortion(3,2)+&
        neighbourListCoords(i,6)*distortion(3,3)
        neighbourListCoords(i,4) = xCB
        neighbourListCoords(i,5) = yCB
        neighbourListCoords(i,6) = zCB
! Store differences
        neighbourListCoords(i,7) = neighbourListCoords(i,4) - neighbourListCoords(i,1)
        neighbourListCoords(i,8) = neighbourListCoords(i,5) - neighbourListCoords(i,2)
        neighbourListCoords(i,9) = neighbourListCoords(i,6) - neighbourListCoords(i,3)
! Store separation
        neighbourListR(i) = (neighbourListCoords(i,7)**2 + &
        neighbourListCoords(i,8)**2 + &
        neighbourListCoords(i,9)**2)**0.5
! Store directon magnitudes
        neighbourListCoords(i,10) = neighbourListCoords(i,7)/neighbourListR(i)
        neighbourListCoords(i,11) = neighbourListCoords(i,8)/neighbourListR(i)
        neighbourListCoords(i,12) = neighbourListCoords(i,9)/neighbourListR(i)
      End Do
    End If
  End Subroutine applyDistortionNL

End Module neighbourList
