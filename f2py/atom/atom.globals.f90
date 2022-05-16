
INTEGER(kind=StandardInteger), PARAMETER ::        mpf = 100 
INTEGER(kind=StandardInteger), PARAMETER ::        lgmax = 20
INTEGER(kind=StandardInteger), PARAMETER ::        tabsize = 1001
CHARACTER*32 , PARAMETER ::                        blank = "                                "

INTEGER(kind=StandardInteger) ::                   fcount = 0

CHARACTER*32 ::                                    pot_fname(1:mpf)
INTEGER(kind=StandardInteger) ::                   pot_is_tab(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                   pot_ftype(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                   pot_a(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                   pot_b(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                   pot_fgroup(1:mpf) = 0
REAL(kind=DoubleReal) ::                           pot_tab(1:mpf, 1:tabsize, 1:3) = 0.0D0

REAL(kind=DoubleReal) ::                           pot_p(1:mpf, 1:100) = 0.0D0
REAL(kind=DoubleReal) ::                           pot_pf(1:mpf, 1:100) = 0.0D0
INTEGER(kind=StandardInteger) ::                   pot_psize(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                   pot_pfsize(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                   group_max = 0

! Pair Index
INTEGER(kind=StandardInteger) ::                   pair_index_a(0:lgmax, 0:lgmax) = -1
INTEGER(kind=StandardInteger) ::                   pair_index_b(0:lgmax * lgmax) = 0

! Dens Index - Any Atom B
INTEGER(kind=StandardInteger) ::                   dens_index_a(0:lgmax, 0:lgmax) = -1 
INTEGER(kind=StandardInteger) ::                   dens_index_b(0:lgmax, 0:lgmax) = 0

! Embe Index
INTEGER(kind=StandardInteger) ::                   embe_index_a(0:lgmax, 0:lgmax) = -1 
INTEGER(kind=StandardInteger) ::                   embe_index_b(0:lgmax, 0:lgmax) = 0

! MD HISTORY
REAL(kind=DoubleReal), ALLOCATABLE ::              chistory(:, :, :)
REAL(kind=DoubleReal), ALLOCATABLE ::              ehistory(:)

