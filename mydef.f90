!----------------------------------------------------------------
! Name : const
! Purpose : Define the constants necessary
!----------------------------------------------------------------
module mydef
  implicit none

  ! precision definition
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)

  ! costant definition
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
  real(dp), parameter :: one = 1.0_dp, zero=0.0_dp, half=0.5_dp, two=2.0_dp
  real(dp), parameter :: e05 = 1.0e-05_dp, e10 = 1.0e-10_dp, e12=1.0E-12_dp

  complex(dpc), parameter :: ic = cmplx(0.0_dp, 1.0_dp)
  complex(dpc), parameter :: itwopi = cmplx(zero, twopi)

end module mydef

