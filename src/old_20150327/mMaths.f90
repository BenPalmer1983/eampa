Module mMaths

! --------------------------------------------------------------!
! MPI Maths functions
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!

! ----------------------------------------
! Updated: 13th August 2014
! ----------------------------------------
! Use standard maths module
  Use msubs
  Use maths
! Force declaration of all variables
  Implicit None
! Include MPI header
  Include 'mpif.h'
! Kinds
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer    1 byte
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer    4 bytes -2E31 to 2E31-1
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer 8 bytes -2E63 to 2E63-1
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer

! Make private
  Private
! --Functions--!
  Public :: M_BirchMurnFit

  Contains

! ------------------------------------------------------------------------!
!                                                                        !
! MODULE FUNCTIONS                                                       !
!                                                                        !
!                                                                        !
! ------------------------------------------------------------------------!

! ------------------------------------------------------------------------!
! Fitting, Regression, Interpolation
! ------------------------------------------------------------------------!

  Function M_BirchMurnFit(points, varianceIn, loopsIn, refinementsIn, coeffsIn) RESULT (coefficients)
! Fit Murnaghan EoS to data
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i, n, pointToVary
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS, loopFactor
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:3) :: coefficientsQ
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), optional :: varianceIn
    Integer(kind=StandardInteger), optional :: loopsIn, refinementsIn
    Real(kind=DoubleReal), Dimension(1:4), optional :: coeffsIn
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops, refinements
! MPI Vars
    Integer(kind=StandardInteger) :: processID, processCount, error
! Call mpi subroutines
    Call MPI_Comm_rank(MPI_COMM_WORLD,processID,error)
    Call MPI_Comm_size(MPI_COMM_WORLD,processCount,error)
! Optional argument variables
    variance = 0.01D0
    loops = 1000
    refinements = 5
    If(present(varianceIn))Then
      variance = varianceIn
    End If
    If(present(loopsIn))Then
      loops = loopsIn
    End If
    If(present(refinementsIn))Then
      refinements = refinementsIn
    End If
! Init variables
    loopFactor = 1.0D0
! Starting point
    If(Present(coeffsIn))Then
      coefficients = coeffsIn
    Else
      If(processID.eq.0)Then
! Quadratic fit  (as a starting point)
        coefficientsQ = PolyFit(points,2)
! Random number
        Call RANDOM_NUMBER(randDouble)
! Starting values for fit
        coefficients(2) = (-1.0D0*coefficientsQ(2))/(2.0D0*coefficientsQ(3))    !V0
        coefficients(1) = coefficientsQ(3)*coefficients(2)**2+&                 !E0
        coefficientsQ(2)*coefficients(2)+&
        coefficientsQ(1)
        coefficients(3) = 2.0D0 * coefficientsQ(3) * coefficients(2)            !B0
        coefficients(4) = 2.0D0 + 2.0D0 * randDouble                            !B'0
      End If
! Distribute
      Call M_distDouble1D(coefficients)
    End If
    If(processID.eq.0)Then  ! Just do on root process
! Starting RSS
      optRSS = M_BirchMurnRSS(points,coefficients)
! Adjust points - second pass
      Do n=0,refinements
        loopFactor = exp(-0.5D0*n)
        Do i=1,loops
          If(n.lt.3)Then
            pointToVary = mod(i-1,3)+1
          Else
            pointToVary = mod(i-1,4)+1
          End If
          coefficientsTemp = coefficients
          Call RANDOM_NUMBER(randDouble)
          varyAmount = variance*2.0D0*(-0.5D0+randDouble)*loopFactor
          coefficientsTemp(pointToVary) = &
          (1.0D0 + varyAmount)*coefficientsTemp(pointToVary)
          testRSS = M_BirchMurnRSS(points,coefficientsTemp)
          If(testRSS.lt.optRSS)Then
            optRSS = testRSS
            coefficients = coefficientsTemp
            If(optRSS.lt.1.0D-7)Then
              Exit
            End If
          End If
        End Do
        If(optRSS.lt.1.0D-7)Then
          Exit
        End If
      End Do
    End If
! Distribute
    Call M_distDouble1D(coefficients)
  End Function M_BirchMurnFit

  Function M_BirchMurnCalc(volume,coefficients) RESULT (energy)
! Calculate energy from volume using Murnaghan EoS
    Implicit None  !Force declaration of all variables
! Declare variables
    Real(kind=DoubleReal) :: volume, energy, eta
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients    ! E0, V0, B0, B'0
! Calculate energy
    eta = (volume/coefficients(2))**(1.0D0/3.0D0)
    energy = coefficients(1) + &
    ((9.0D0*coefficients(3)*coefficients(2))/(16.0D0))*&
    ((eta**2-1.0D0)**2)*&
    (6.0D0+coefficients(4)*(eta**2-1.0D0)-4.0D0*eta**2)
  End Function M_BirchMurnCalc

  Function M_BirchMurnRSS(points,coefficients) RESULT (rss)
! Fit Murnaghan EoS to data
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients    ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: volume, energy, energyC, rss
! calculate RSS
    rss = 0.0D0
    Do i=1,size(points,1)
      volume = points(i,1)
      energy = points(i,2)
      energyC = M_BirchMurnCalc(volume,coefficients)
      rss = rss + (energyC-energy)**2
    End Do
  End Function M_BirchMurnRSS

  Function M_BirchMurnFitBP(points, coefficientsIn) RESULT (coefficients)
! Fits B'0 and holds other coefficients
    Implicit None  !Force declaration of all variables
! Declare variables
    Integer(kind=StandardInteger) :: i
! Real(kind=DoubleReal) :: energyMin, volOpt, bm, bmP, randDouble
    Real(kind=DoubleReal) :: randDouble, varyAmount, optRSS, testRSS
    Real(kind=DoubleReal), Dimension(:,:) :: points
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsIn  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficients  ! E0, V0, B0, B'0
    Real(kind=DoubleReal), Dimension(1:4) :: coefficientsTemp  ! E0, V0, B0, B'0
    Real(kind=DoubleReal) :: variance
    Integer(kind=StandardInteger) :: loops
! Optional argument variables
    variance = 0.02D0
    loops = 100
    coefficients = coefficientsIn
! Starting RSS
    optRSS = M_BirchMurnRSS(points,coefficients)
! Adjust points
    Do i=1,loops
      coefficientsTemp = coefficients
      Call RANDOM_NUMBER(randDouble)
      varyAmount = variance*2.0D0*(-0.5D0+randDouble)
      coefficientsTemp(4) = &
      (1.0D0 + varyAmount)*coefficientsTemp(4)
      testRSS = M_BirchMurnRSS(points,coefficientsTemp)
      If(testRSS.lt.optRSS)Then
        optRSS = testRSS
        coefficients = coefficientsTemp
        If(optRSS.lt.1.0D-5)Then
          Exit
        End If
      End If
    End Do
  End Function M_BirchMurnFitBP

End Module mMaths
