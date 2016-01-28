Module evalBP
! --------------------------------------------------------------!
! Evaluate configurations/EAM
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Calls the calcEAM subroutines and evaluates the results
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
  Use plotTypes
  Use plot
  Use globals
  Use initialise
  Use loadData
  Use output
  Use eamGen
  Use neighbourListBP
  Use bpCalcEAM
! Force declaration of all variables
  Implicit None
! Privacy of variables/functions/subroutines
  Private
! Public Subroutines
  Public :: evalBulkProperties
  Public :: evalBP_PrintRSS
  Contains
! ---------------------------------------------------------------------------------------------------  
  Subroutine evalBulkProperties(printRssIn,makeChartIn)
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID
    Logical, optional :: printRssIn
    Logical, optional :: makeChartIn
    Logical :: printRss
    Logical :: makeChart
! Optional arguments
    printRss = .false.
    makeChart = .false.
    If(Present(printRssIn))Then
      printRss = printRssIn
    End If
    If(Present(makeChartIn))Then
      makeChart = makeChartIn
    End If
! Store the neighbour list to recall after distortion
    Call bpStoreUndistortedNL()    
! Calculate alat (and fit EoS)
    Call bpEoS(makeChart) !evalBP.f90 
! Calculate elastic constants    
    Call bpEC(makeChart) !evalBP.f90       
! Loop through configs and calculate RSS
    If(mpiProcessID.eq.0)Then
      rssBPArrTotal%total = 0.0D0
      rssBPArrTotal%aLat = 0.0D0
      rssBPArrTotal%v0 = 0.0D0
      rssBPArrTotal%e0 = 0.0D0
      rssBPArrTotal%b0 = 0.0D0
      rssBPArrTotal%bp0 = 0.0D0
      rssBPArrTotal%eos = 0.0D0
      rssBPArrTotal%c11 = 0.0D0
      rssBPArrTotal%c12 = 0.0D0
      rssBPArrTotal%c44 = 0.0D0
      Do configID=1,configCountBP
        Call evalBP_RSS(configID)
! Totals      
        rssBPArrTotal%total = rssBPArrTotal%total+rssBPArr(configID)%total
        rssBPArrTotal%aLat = rssBPArrTotal%aLat+rssBPArr(configID)%aLat
        rssBPArrTotal%v0 = rssBPArrTotal%v0+rssBPArr(configID)%v0
        rssBPArrTotal%e0 = rssBPArrTotal%e0+rssBPArr(configID)%e0
        rssBPArrTotal%b0 = rssBPArrTotal%b0+rssBPArr(configID)%b0
        rssBPArrTotal%bp0 = rssBPArrTotal%bp0+rssBPArr(configID)%bp0
        rssBPArrTotal%eos = rssBPArrTotal%eos+rssBPArr(configID)%eos
        rssBPArrTotal%c11 = rssBPArrTotal%c11+rssBPArr(configID)%c11
        rssBPArrTotal%c12 = rssBPArrTotal%c12+rssBPArr(configID)%c12
        rssBPArrTotal%c44 = rssBPArrTotal%c44+rssBPArr(configID)%c44
! Total RSS
        totalRSS = totalRSS + rssBPArr(configID)%total
      End Do
    End If
! Distribute value
    Call M_distDouble(totalRSS)
! Print   
    If(printRss)Then
      Call evalBP_PrintRSS()
    End If  
  End Subroutine evalBulkProperties
! ---------------------------------------------------------------------------------------------------
  Subroutine bpStoreUndistortedNL()
! Store neighbour list for each bulk property config 
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Do configID=1,configCountBP    
      Call bpStoreNL(configID)
    End Do
  End Subroutine bpStoreUndistortedNL
! ---------------------------------------------------------------------------------------------------
  Subroutine bpEoS(makeChart)
! Calculates the lattice parameter - also gives the optimum energy, bulk modulus and derivative 
! using Birch Murnaghan fit  
    Implicit None   ! Force declaration of all variables
! Private variables
    Logical :: makeChart
    Real(kind=DoubleReal) :: dLatInc
    Real(kind=DoubleReal) :: bpTimeStart, bpTimeEnd
! Start time    
    Call cpu_time(bpTimeStart)
! One attempt
    dLatInc = 0.02D0
    Call bpEoSProcess(dLatInc, makeChart)
! End time    
    Call cpu_time(bpTimeEnd)
    evalTimeBP = evalTimeBP + bpTimeEnd - bpTimeStart
! Results stored in calcBulkProperties
  End Subroutine bpEoS
! ---------------------------------------------------------------------------------------------------
  Subroutine bpEoSProcess(dLatInc, makeChart)
! Calculates the lattice parameter - also gives the optimum energy, bulk modulus and derivative 
! using Birch Murnaghan fit  
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, i, j, k 
    Integer(kind=StandardInteger) :: unitCopies
    Real(kind=DoubleReal) :: alat, dAlat, aLatCalc, dLatInc
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: dMatrix
    Real(kind=DoubleReal), Dimension(1:aLatSamples,1:2) :: energyVolume, energyAlat
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBM
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: chartPoints    
    Type(plotData) :: bpPlotData
    Logical :: makeChart
    Character(len=32) :: fileName
    Character(len=32) :: tempStr
    Character(len=8) :: tempName  
! Init variables
    dMatrix = 0.0D0
    energyVolume = 0.0D0       
    k = 0
    Do i=LowerInc(aLatSamples),UpperInc(aLatSamples)
      k = k + 1
! Change in aLat
      dAlat = 1.0D0*i*dLatInc
! Volume per atom
      Do configID=1,configCountBP  
        alat = bpInArr(configID)%alat
        unitCopies = bpInArr(configID)%size
        configVolBP(configID,k) = (alat * unitCopies * (1.0D0+dAlat))**3/&
        configurationCoordsKeyBP(configID,2)  
      End Do        
! Load config nl points
      If(i.gt.LowerInc(aLatSamples))Then
        Do configID=1,configCountBP  
          Call bpLoadNL(configID)
        End Do  
      End If  
! Prepare distortion matrix
      dMatrix = 0.0D0
      Do j=1,3
        dMatrix(j,j) = 1.0D0+dAlat
      End Do
! Apply distortion matrix to each config
      Do configID=1,configCountBP  
        Call applyDistortion(configID, dMatrix)
      End Do
! Calc energies of all configurations with distortion applied      
! Need to sort out MPI optimisation better here
      Call calcEnergiesBP()  ! bpCalcEAM.f90
! Save undistorted energy for use later     
      If(i.eq.0)Then
        Do configID=1,configCountBP  
          undistortedCellEnergies(configID) = configCalcEnergiesBP(configID)
        End Do  
      End If  
! Store Data - root process only
      If(mpiProcessID.eq.0)Then
        Do configID=1,configCountBP  
          alatEnergies(configID,k) = configCalcEnergiesBP(configID)
        End Do
      End If  
    End Do
! Fit energy data with bm fit (root process only)
    If(mpiProcessID.eq.0)Then
      Do configID=1,configCountBP 
        energyVolume = 0.0D0    
        Do i=1,aLatSamples
          energyVolume(i,1) = configVolBP(configID,i)
          energyVolume(i,2) = alatEnergies(configID,i)
          energyAlat(i,1) = (bpInArr(configID)%atomsPerUnit*configVolBP(configID,i))**(1.0D0/3.0D0)
          energyAlat(i,2) = alatEnergies(configID,i)
        End Do  
! Fit points to birch-murn           
        coefficientsBM = BirchMurnFit(energyVolume, 1.0D0, 9.0D0)         
! Make chart if required
        If(makeChart)Then
! Chart settings   
          Call plotInit(bpPlotData)
          tempName = RandName()
          write(tempStr,"(I4)") configID          
          fileName = "EoS_"//trim(adjustl(tempStr))//"_"//tempName
          bpPlotData%tempDirectory = trim(tempDirectory)
          bpPlotData%outputDirectory = trim(outputDirectory)
          bpPlotData%outputName = trim(fileName)
          bpPlotData%title = "Equation of State"    
          bpPlotData%xAxis = "Unit Cell Alat/Ang"    
          bpPlotData%yAxis = "Energy per Atom/eV" 
          bpPlotData%cleanPyFile = .true.
          bpPlotData%dataFile = .true.
! Add data
          Call plotAdd(bpPlotData, energyAlat,"","BM1,BM2")
          Call plotStyle(bpPlotData,"o","--")
! Make
          Call plotMake(bpPlotData)             
        End If  
        aLatCalc = (bpInArr(configID)%atomsPerUnit*coefficientsBM(2))**(1.0D0/3.0D0)     
! ---- Store results
        calcBulkProperties(configID)%alat = aLatCalc
        calcBulkProperties(configID)%v0 = coefficientsBM(2)
        calcBulkProperties(configID)%e0 = coefficientsBM(1)
        calcBulkProperties(configID)%b0 = coefficientsBM(3)
        calcBulkProperties(configID)%bp0 = coefficientsBM(4)
! ---- Output data
        If(bpPrintData)Then   
          Call outputBpData(configID,energyVolume)
        End If  
      End Do
    End If    
  End Subroutine bpEoSProcess
! ---------------------------------------------------------------------------------------------------
  Subroutine bpEC(makeChart)
! Elastic Constants 
! Method of M Mehl 1990
! calcBulkProperties(configID)%alat, calcBulkProperties(configID)%b0 
! using Birch Murnaghan fit  
    Implicit None   ! Force declaration of all variables
! Private variables
    Logical :: makeChart
    Real(kind=DoubleReal) :: dLatInc
    Real(kind=DoubleReal) :: bpTimeStart, bpTimeEnd
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: dMatrix
    Integer(kind=StandardInteger) :: configID
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal) :: dStrain, strain
    Real(kind=DoubleReal), Dimension(1:ecSamples,1:2) :: dataPoints
    Real(kind=DoubleReal), Dimension(1:(2*ecSamples+1),1:2) :: dataPointsT
    Real(kind=DoubleReal), Dimension(1:maxConfigsBP) :: tetragonal, orthorhombic, monoclinic
    Real(kind=DoubleReal), Dimension(1:3) :: coefficients
    Real(kind=DoubleReal) :: volume, orthVal, tetraVal, monoclinicVal, bulkVal, avgOrthTetra
    Type(plotData) :: ecPlotData
    Character(len=32) :: fileName
    Character(len=32) :: tempStr
! Start time    
    Call cpu_time(bpTimeStart)
! Elastic constants are calculated for Cubic Crystals only - C11 C12 C44
!
! Init energy-strain arrays
    ecEnergiesBP = 0.0D0
    ecStrainBP = 0.0D0
    Do configID=1,configCountBP
      ecEnergiesBP(configID,1) = undistortedCellEnergies(configID)
    End Do    
    dStrain = 0.01D0
! ---------------------------------------------------
! Tetragonal
! ---------------------------------------------------
! c11-c12 = c2/(3V)  
    k = 0  
    Do i=(-1*ecSamples),ecSamples
      k = k + 1
      If(i.eq.0)Then
        Do configID=1,configCountBP  
          ecEnergiesBP_T(configID,k) = undistortedCellEnergies(configID)
          ecStrainBP_T(configID,k) = 0.0D0
        End Do  
      Else
        strain = (i)*dStrain
! Load       
        Do configID=1,configCountBP  
          Call bpLoadNL(configID)
        End Do      
! Prepare distortion matrix
        dMatrix = 0.0D0
        Do j=1,3
          dMatrix(j,j) = 1.0D0
        End Do
        dMatrix(1,1) = dMatrix(1,1) + strain
        dMatrix(2,2) = dMatrix(2,2) + strain
        dMatrix(3,3) = dMatrix(3,3) + (1.0D0+strain)**(-2) - 1.0D0      
! Apply distortion matrix to each config
        Do configID=1,configCountBP  
          Call applyDistortion(configID, dMatrix)
        End Do
! calculate energies      
        Call calcEnergiesBP()  ! bpCalcEAM.f90      
        Do configID=1,configCountBP  
          ecEnergiesBP_T(configID,k) = configCalcEnergiesBP(configID)
          ecStrainBP_T(configID,k) = strain
        End Do      
      End If
    End Do  
! Store data points and fit (make chart if requested)    
    Do configID=1,configCountBP
      Do i=1,(2*ecSamples+1)
        dataPointsT(i,1) = ecStrainBP_T(configID,i)
        dataPointsT(i,2) = ecEnergiesBP_T(configID,i)
      End Do    
      coefficients = PolyFit(dataPointsT,2)      
      tetragonal(configID) = coefficients(3)
! Make chart if required
      If(makeChart)Then
! Chart settings   
        Call plotInit(ecPlotData)
        write(tempStr,"(I4)") configID          
        fileName = "EC_Tetragonal_"//trim(adjustl(tempStr))
        ecPlotData%tempDirectory = trim(tempDirectory)
        ecPlotData%outputDirectory = trim(outputDirectory)
        ecPlotData%outputName = trim(fileName)
        ecPlotData%title = "Equation of State"    
        ecPlotData%xAxis = "Unit Cell Alat/Ang"    
        ecPlotData%yAxis = "Energy per Atom/eV" 
        ecPlotData%cleanPyFile = .true.
        ecPlotData%dataFile = .true.
! Add data
        Call plotAdd(ecPlotData, dataPointsT,"","POLY2")
        Call plotStyle(ecPlotData,"o","--")
! Make
        Call plotMake(ecPlotData)             
      End If 
    End Do
! ---------------------------------------------------
! Orthorhombic
! ---------------------------------------------------
! c11-c12 = c2/V    
    Do i=2,ecSamples
      strain = (i-1)*dStrain
! Load       
      Do configID=1,configCountBP  
        Call bpLoadNL(configID)
      End Do      
! Prepare distortion matrix
      dMatrix = 0.0D0
      Do j=1,3
        dMatrix(j,j) = 1.0D0
      End Do
      dMatrix(1,1) = dMatrix(1,1) + strain
      dMatrix(2,2) = dMatrix(2,2) - strain
      dMatrix(3,3) = dMatrix(3,3) + strain**2/(1-strain**2)
! Apply distortion matrix to each config
      Do configID=1,configCountBP  
        Call applyDistortion(configID, dMatrix)
      End Do
! calculate energies      
      Call calcEnergiesBP()  ! bpCalcEAM.f90      
      Do configID=1,configCountBP  
        ecEnergiesBP(configID,i) = configCalcEnergiesBP(configID)
        ecStrainBP(configID,i) = strain
      End Do      
    End Do  
! Store data points and fit (make chart if requested)    
    Do configID=1,configCountBP
      Do i=1,ecSamples
        dataPoints(i,1) = ecStrainBP(configID,i)
        dataPoints(i,2) = ecEnergiesBP(configID,i)
      End Do    
      coefficients = PolyFit(dataPoints,2)      
      orthorhombic(configID) = coefficients(3)
! Make chart if required
      If(makeChart)Then
! Chart settings   
        Call plotInit(ecPlotData)
        write(tempStr,"(I4)") configID          
        fileName = "EC_Orthorhombic_"//trim(adjustl(tempStr))
        ecPlotData%tempDirectory = trim(tempDirectory)
        ecPlotData%outputDirectory = trim(outputDirectory)
        ecPlotData%outputName = trim(fileName)
        ecPlotData%title = "Equation of State"    
        ecPlotData%xAxis = "Unit Cell Alat/Ang"    
        ecPlotData%yAxis = "Energy per Atom/eV" 
        ecPlotData%cleanPyFile = .true.
        ecPlotData%dataFile = .true.
! Add data
        Call plotAdd(ecPlotData, dataPoints,"","POLY2")
        Call plotStyle(ecPlotData,"o","--")
! Make
        Call plotMake(ecPlotData)             
      End If 
    End Do
! ---------------------------------------------------
! Monoclinic
! ---------------------------------------------------
! c44 = 2c2/V    
    Do i=2,ecSamples
      strain = (i-1)*dStrain
! Load       
      Do configID=1,configCountBP  
        Call bpLoadNL(configID)
      End Do      
! Prepare distortion matrix
      dMatrix = 0.0D0
      Do j=1,3
        dMatrix(j,j) = 1.0D0
      End Do
      dMatrix(2,1) = dMatrix(2,1) + 0.5D0*strain
      dMatrix(1,2) = dMatrix(1,2) + 0.5D0*strain
      dMatrix(3,3) = dMatrix(3,3) + strain**2/(4.0D0-strain**2)
! Apply distortion matrix to each config
      Do configID=1,configCountBP  
        Call applyDistortion(configID, dMatrix)
      End Do
! calculate energies      
      Call calcEnergiesBP()  ! bpCalcEAM.f90      
      Do configID=1,configCountBP  
        ecEnergiesBP(configID,i) = configCalcEnergiesBP(configID)
        ecStrainBP(configID,i) = strain
      End Do      
    End Do  
! Store data points and fit (make chart if requested)    
    Do configID=1,configCountBP
      Do i=1,ecSamples
        dataPoints(i,1) = ecStrainBP(configID,i)
        dataPoints(i,2) = ecEnergiesBP(configID,i)
      End Do    
      coefficients = PolyFit(dataPoints,2)      
      monoclinic(configID) = coefficients(3)
! Make chart if required
      If(makeChart)Then
! Chart settings   
        Call plotInit(ecPlotData)
        write(tempStr,"(I4)") configID          
        fileName = "EC_Monoclinic_"//trim(adjustl(tempStr))
        ecPlotData%tempDirectory = trim(tempDirectory)
        ecPlotData%outputDirectory = trim(outputDirectory)
        ecPlotData%outputName = trim(fileName)
        ecPlotData%title = "Equation of State"    
        ecPlotData%xAxis = "Unit Cell Alat/Ang"    
        ecPlotData%yAxis = "Energy per Atom/eV" 
        ecPlotData%cleanPyFile = .true.
        ecPlotData%dataFile = .true.
! Add data
        Call plotAdd(ecPlotData, dataPoints,"","POLY2")
        Call plotStyle(ecPlotData,"o","--")
! Make
        Call plotMake(ecPlotData)             
      End If 
    End Do    
! Calculate elastic constants
    Do configID=1,configCountBP
      volume = calcBulkProperties(configID)%v0
      bulkVal = calcBulkProperties(configID)%b0
      tetraVal = tetragonal(configID)/(3.0D0*volume)
      orthVal = orthorhombic(configID)/(1.0D0*volume)
      avgOrthTetra = 0.5D0*(tetraVal+orthVal)
      monoclinicVal = (2.0D0*monoclinic(configID))/(1.0D0*volume)
      calcBulkProperties(configID)%c12 = (3.0D0*bulkVal-avgOrthTetra)/3.0D0
      calcBulkProperties(configID)%c11 = avgOrthTetra+calcBulkProperties(configID)%c12
      calcBulkProperties(configID)%c44 = monoclinicVal
      calcBulkProperties(configID)%shearConstant = (calcBulkProperties(configID)%c11-&
                                                   calcBulkProperties(configID)%c12)/2.0D0
      
    End Do  
    
! End time    
    Call cpu_time(bpTimeEnd)
    evalTimeBP = evalTimeBP + bpTimeEnd - bpTimeStart
! Results stored in calcBulkProperties
  End Subroutine bpEC  
  
! ---------------------------------------------------------------------------------------------------
  Subroutine evalBP_RSS(configID)
! calculate rss from bulk property calculated values
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: aLatRSS, v0RSS, e0RSS, b0RSS, bp0RSS
    Real(kind=DoubleReal) :: c11RSS, c12RSS, c44RSS, eosRSS 
    Type(eos) :: refEoS, calcEoS
! Init
    aLatRSS = 0.0D0
    v0RSS = 0.0D0
    e0RSS = 0.0D0
    b0RSS = 0.0D0
    bp0RSS = 0.0D0
    c11RSS = 0.0D0
    c12RSS = 0.0D0
    c44RSS = 0.0D0
! calculate each rss
    rssBPArr(configID)%alat = RSSCalc(bpInArr(configID)%alat,calcBulkProperties(configID)%alat,rssWeighting(4))
    rssBPArr(configID)%v0 = RSSCalc(bpInArr(configID)%v0,calcBulkProperties(configID)%v0,rssWeighting(4))
    rssBPArr(configID)%e0 = RSSCalc(bpInArr(configID)%e0,calcBulkProperties(configID)%e0,rssWeighting(5))
    rssBPArr(configID)%b0 = RSSCalc(bpInArr(configID)%b0,calcBulkProperties(configID)%b0,rssWeighting(6))
    rssBPArr(configID)%bp0 = RSSCalc(bpInArr(configID)%bp0,calcBulkProperties(configID)%bp0,1.0D0)
    rssBPArr(configID)%c11 = RSSCalc(bpInArr(configID)%c11,calcBulkProperties(configID)%c11,rssWeighting(7))
    rssBPArr(configID)%c12 = RSSCalc(bpInArr(configID)%c12,calcBulkProperties(configID)%c12,rssWeighting(7))
    rssBPArr(configID)%c44 = RSSCalc(bpInArr(configID)%c44,calcBulkProperties(configID)%c44,rssWeighting(7))
    
    
! Zero B'0 - input is just a guess at the moment
    rssBPArr(configID)%bp0 = 0.0D0
! Check NaN
    If(ISNAN(calcBulkProperties(configID)%alat))Then
      rssBPArr(configID)%alat = 1.0D20
    End If
    If(ISNAN(calcBulkProperties(configID)%v0))Then
      rssBPArr(configID)%v0 = 1.0D20
    End If
    If(ISNAN(calcBulkProperties(configID)%e0))Then
      rssBPArr(configID)%e0 = 1.0D20
    End If
    If(ISNAN(calcBulkProperties(configID)%b0))Then
      rssBPArr(configID)%b0 = 1.0D20
    End If
! EoS RSS    
! Set up data types to compare                 
    refEoS%e0 = bpInArr(configID)%e0
    refEoS%v0 = bpInArr(configID)%v0   
    refEoS%b0 = bpInArr(configID)%b0
    refEoS%bp0 = bpInArr(configID)%bp0
    calcEoS%e0 = calcBulkProperties(configID)%e0
    calcEoS%v0 = calcBulkProperties(configID)%v0
    calcEoS%b0 = calcBulkProperties(configID)%b0
    calcEoS%bp0 = calcBulkProperties(configID)%bp0
! calculate and store EoS RSS
    Call evalBP_EosRSS(refEoS, calcEoS, eosRSS)  
    rssBPArr(configID)%eos = eosRSS     
! Sum
    rssBPArr(configID)%total = rssBPArr(configID)%alat+rssBPArr(configID)%v0+&
    rssBPArr(configID)%e0+rssBPArr(configID)%b0+rssBPArr(configID)%bp0+&
    rssBPArr(configID)%c11+rssBPArr(configID)%c12+rssBPArr(configID)%c44+&
    rssBPArr(configID)%eos    
  End Subroutine evalBP_RSS
  
  Subroutine evalBP_EosRSS(refEoS, calcEoS, rss)
! calculate rss between reference EoS and calculated
    Implicit None   ! Force declaration of all variables
! Arg
    Type(eos) :: refEoS, calcEoS
    Real(kind=DoubleReal) :: rss, vol, increment, eCalc, eRef
    Real(kind=DoubleReal), Dimension(1:4) :: calcEoS_Arr, refEoS_Arr
! Private variables
    Integer(kind=StandardInteger) :: i
! Init vars    
    vol = 0.80D0*refEoS%v0
    increment = (0.40D0*refEoS%v0)/99.0D0
    calcEoS_Arr(1) = calcEoS%e0
    calcEoS_Arr(2) = calcEoS%v0
    calcEoS_Arr(3) = calcEoS%b0
    calcEoS_Arr(4) = calcEoS%bp0
    refEoS_Arr(1) = refEoS%e0
    refEoS_Arr(2) = refEoS%v0
    refEoS_Arr(3) = refEoS%b0
    refEoS_Arr(4) = refEoS%bp0
    rss = 0.0D0
! Loop through comparison points
    Do i=1,100    
! increment for next loop
      vol = vol + increment
      eCalc = BirchMurnCalc(vol,calcEoS_Arr)
      eRef = BirchMurnCalc(vol,refEoS_Arr)
      rss = rss + RSSCalc(eCalc,eRef,rssWeighting(8))
      !rss = rss + (eCalc-eRef)**2
    End Do  
  End Subroutine evalBP_EosRSS
  
  
  Subroutine evalBP_PrintRSS()
! calculate rss between reference EoS and calculated
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
! Print
    If(mpiProcessID.eq.0)Then
      print *,"eval rss"
      Do configID=1,configCountBP  
        print *,"RSS Values: BP Config ",configID
        print *,"alat   ",rssBPArr(configID)%aLat
        print *,"v0     ",rssBPArr(configID)%v0
        print *,"e0     ",rssBPArr(configID)%e0
        print *,"b0     ",rssBPArr(configID)%b0
        print *,"bp0    ",rssBPArr(configID)%bp0
        print *,"eos    ",rssBPArr(configID)%eos
        print *,"c11    ",rssBPArr(configID)%c11
        print *,"c12    ",rssBPArr(configID)%c12
        print *,"c44    ",rssBPArr(configID)%c44
        print *,"Total  ",rssBPArr(configID)%total
      End Do
    End If    
  End Subroutine evalBP_PrintRSS

! ------------------------------------------------------------------------!
!                                                                         !
! MODULE FUNCTIONS                                                        !
!                                                                         !
!                                                                         !
! ------------------------------------------------------------------------!
  Function LowerInc(incs) RESULT (lower)
! Get value of function at x
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: incs, lower
    If(Odd(incs))Then
      lower = -1*(incs-1)/2
    Else
      lower = -1*(incs-2)/2
    End If
  End Function LowerInc
! ---------------------------------------------------------------------------------------------------
  Function UpperInc(incs) RESULT (upper)
! Get value of function at x
    Implicit None ! Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: incs, upper
    If(Odd(incs))Then
      upper = (incs-1)/2
    Else
      upper = incs/2
    End If    
  End Function UpperInc
! ---------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------------
End Module evalBP
