MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37
  INTEGER, PARAMETER :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307
  INTEGER, PARAMETER :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  INTEGER, PARAMETER :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer         1 byte    -2^4 to 2^4-1                                           -16 to 16
  INTEGER, PARAMETER :: SmallInteger = Selected_Int_Kind(4)           ! small integer        4 bytes  -2^16 to 2^16-1                                       -65536 to 65535
  INTEGER, PARAMETER :: StandardInteger = Selected_Int_Kind(8)        ! standard integer     8 bytes  -2^31 to 2^31-1                                  -2147483648 to 2147483648
  INTEGER, PARAMETER :: LongInteger = Selected_Int_Kind(12)           ! long integer         16 bytes -2^63 to 2^63-1                         -9223372036854775808 to 9223372036854775807
  INTEGER, PARAMETER :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer    32 bytes -2^63 to 2^63-1     -170141183460469231731687303715884105728 to 170141183460469231731687303715884105727
END MODULE kinds
