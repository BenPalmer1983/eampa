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
  Public :: bpAlat
  Contains
! ---------------------------------------------------------------------------------------------------  
  Subroutine evalBulkProperties()
    Implicit None   ! Force declaration of all variables
    Integer(kind=StandardInteger) :: configID
! Calculate alat (and fit EoS)
    Call bpEoS() !evalBP.f90 
! Loop through configs and calculate RSS
    If(mpiProcessID.eq.0)Then
      Do configID=1,configCountBP
        Call evalBP_RSS(configID)
        totalRSS = totalRSS + rssBPArr(configID)%total
      End Do
    End If
! Distribute value
    Call M_distDouble(totalRSS)


  End Subroutine evalBulkProperties
! ---------------------------------------------------------------------------------------------------
  Subroutine bpEoS()
! Calculates the lattice parameter - also gives the optimum energy, bulk modulus and derivative 
! using Birch Murnaghan fit  
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: dLatInc
    Logical :: makeChart
    Real(kind=DoubleReal) :: bpTimeStart, bpTimeEnd
! Start time    
    Call cpu_time(bpTimeStart)
! One attempt
    aLatBP(1) = 4.04D0
    makeChart = .false.
    dLatInc = 0.02D0
    Call bpAlatProcess(dLatInc, makeChart)
! End time    
    Call cpu_time(bpTimeEnd)
    evalTimeBP = evalTimeBP + bpTimeEnd - bpTimeStart
! Results stored in calcBulkProperties
  End Subroutine bpEoS
! ---------------------------------------------------------------------------------------------------
  Subroutine bpAlat()
! Calculates the lattice parameter - also gives the optimum energy, bulk modulus and derivative 
! using Birch Murnaghan fit  
    Implicit None   ! Force declaration of all variables
! Private variables
    Real(kind=DoubleReal) :: dLatInc
    Logical :: makeChart
    Real(kind=DoubleReal) :: bpTimeStart, bpTimeEnd
! Start time    
    Call cpu_time(bpTimeStart)
! First attempt
    makeChart = .true.
    dLatInc = 0.05D0
    Call bpAlatProcess(dLatInc, makeChart)
! Refine, second attempt
    makeChart = .true.
    If(eosChart)Then
      makeChart = .true.   ! Override if set in user input file
    End If
    dLatInc = 0.01D0
    Call bpAlatProcess(dLatInc, makeChart)  
! End time    
    Call cpu_time(bpTimeEnd)
    evalTimeBP = evalTimeBP + bpTimeEnd - bpTimeStart
! Results stored in calcBulkProperties
  End Subroutine bpAlat
! ---------------------------------------------------------------------------------------------------
  Subroutine bpAlatProcess(dLatInc, makeChart)
! Calculates the lattice parameter - also gives the optimum energy, bulk modulus and derivative 
! using Birch Murnaghan fit  
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID, i, j, k 
    Real(kind=DoubleReal) :: dAlat, aLatCalc, dLatInc
    Real(kind=DoubleReal), Dimension(1:3,1:3) :: dMatrix
    Real(kind=DoubleReal), Dimension(1:aLatSamples,1:2) :: energyVolume
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsBM
    Real(kind=DoubleReal), Dimension(1:100,1:2) :: chartPoints
    Type(chart) :: alatChart
    Logical :: makeChart
    Character(len=32) :: fileName
    Character(len=32) :: tempStr
    Character(len=8) :: tempName
! Init variables
    dMatrix = 0.0D0
    energyVolume = 0.0D0       
! Store neighbour list for each bulk property config 
    Do configID=1,configCountBP    
      Call bpStoreNL(configID)
    End Do
    k = 0
    Do i=LowerInc(aLatSamples),UpperInc(aLatSamples)
      k = k + 1
! Change in aLat
      dAlat = 1.0D0*i*dLatInc
! Volume per atom
      Do configID=1,configCountBP  
        configVolBP(configID,k) = (aLatBP(configID) * bpCopies * (1.0D0+dAlat))**3/&
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
      Call calcEnergiesBP()
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
        End Do  
! Fit points to birch-murn   
        coefficientsBM = BirchMurnFit_R(energyVolume)   ! maths.f90
! Make chart if required
        If(makeChart)Then
! Make points
          chartPoints = BirchMurnPoints(coefficientsBM,&
          energyVolume(1,1),energyVolume(aLatSamples,1),100)   ! maths.f90
! Convert points to ang (alat not volume)        
          Do i=1,100
            chartPoints(i,1) = (bpUnitCellCount(configID)*chartPoints(i,1))**(1.0D0/3.0D0)  
          End Do
! Chart labels
          alatChart%title = "Equation of State"    
          alatChart%xAxis = "Unit Cell Alat/Ang"    
          alatChart%yAxis = "Energy per Atom/eV"    
          Call randFileName(tempName)
          write(tempStr,"(I4)") configID          
          fileName = "EoS_"//trim(adjustl(tempStr))//"_"//tempName
          Call makePlot(outputDirectory, trim(fileName), tempDirectory, &
          chartPoints, 1, 2, 1, 100, alatChart, .false.)
        End If  
        aLatCalc = (bpUnitCellCount(configID)*coefficientsBM(2))**(1.0D0/3.0D0)     
! ---- Store results
        calcBulkProperties(configID)%alat = aLatCalc
        calcBulkProperties(configID)%v0 = coefficientsBM(2)
        calcBulkProperties(configID)%e0 = coefficientsBM(1)
        calcBulkProperties(configID)%b0 = UnitConvert(coefficientsBM(3),"EVAN3","GPA")
        calcBulkProperties(configID)%bp0 = coefficientsBM(4)
! ---- Update lattice parameter
        !aLatBP(configID) = aLatCalc     
      End Do
    End If    
  End Subroutine bpAlatProcess
! ---------------------------------------------------------------------------------------------------
  Subroutine evalBP_RSS(configID)
! calculate rss from bulk property calculated values
    Implicit None   ! Force declaration of all variables
! Private variables
    Integer(kind=StandardInteger) :: configID
    Real(kind=DoubleReal) :: aLatRSS, v0RSS, e0RSS, b0RSS, bp0RSS
    Real(kind=DoubleReal) :: c11RSS, c12RSS, c44RSS
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
    aLatRSS = rssWeighting(4)*RSSCalc(refBulkProperties(configID)%alat,calcBulkProperties(configID)%alat)
    v0RSS = rssWeighting(4)*RSSCalc(refBulkProperties(configID)%v0,calcBulkProperties(configID)%v0)
    e0RSS = rssWeighting(5)*RSSCalc(refBulkProperties(configID)%e0,calcBulkProperties(configID)%e0)
    b0RSS = rssWeighting(6)*RSSCalc(refBulkProperties(configID)%b0,calcBulkProperties(configID)%b0)
    bp0RSS = 1.0D0*RSSCalc(refBulkProperties(configID)%bp0,calcBulkProperties(configID)%bp0)
    c11RSS = rssWeighting(7)*RSSCalc(refBulkProperties(configID)%c11,calcBulkProperties(configID)%c11)
    c12RSS = rssWeighting(7)*RSSCalc(refBulkProperties(configID)%c12,calcBulkProperties(configID)%c12)
    c44RSS = rssWeighting(7)*RSSCalc(refBulkProperties(configID)%c44,calcBulkProperties(configID)%c44)
! Convert RSS of GPA values to eVang3
    b0RSS = UnitConvert(b0RSS, "GPA", "EVAN3")
    bp0RSS = UnitConvert(bp0RSS, "GPA", "EVAN3")
    c11RSS = UnitConvert(c11RSS, "GPA", "EVAN3")
    c12RSS = UnitConvert(c12RSS, "GPA", "EVAN3")
    c44RSS = UnitConvert(c44RSS, "GPA", "EVAN3")
! Store
    rssBPArr(configID)%alat = aLatRSS
    rssBPArr(configID)%v0 = v0RSS
    rssBPArr(configID)%e0 = e0RSS
    rssBPArr(configID)%b0 = b0RSS   
    rssBPArr(configID)%bp0 = bp0RSS
    rssBPArr(configID)%c11 = c11RSS
    rssBPArr(configID)%c12 = c12RSS
    rssBPArr(configID)%c44 = c44RSS
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
    
! Sum
    rssBPArr(configID)%total = aLatRSS+e0RSS+b0RSS+c11RSS+c12RSS+c44RSS  
    
  End Subroutine evalBP_RSS
  

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
  Function RSSCalc(inputA, inputB) RESULT (output)
! Get value of function at x
    Implicit None ! Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: inputA, inputB, output
    output = 0.0D0
    If(inputA.gt.-2.0D20.and.inputB.gt.-2.0D20)Then
      output = (inputA-inputB)**2
    End If    
  End Function RSSCalc
! ---------------------------------------------------------------------------------------------------
End Module evalBP
