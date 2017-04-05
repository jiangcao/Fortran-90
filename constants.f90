module constants
use types, only: dp
implicit none
! Constants contain more digits than double precision, so that
! they are rounded correctly:
real(dp), parameter :: pi   = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e    = 2.7182818284590452353602874713527_dp
real(dp), parameter :: BOLTZ = 8.61734e-05_dp              !eV K-1
REAL(dp), PARAMETER :: hbar  = 6.58211899E-16_dp           !eV s
REAL(dp), PARAMETER :: m0    = 5.6856E-16_dp               !eV s^2 / cm^2  rest mass of electron
real(dp), parameter :: elch  = 1.60217653e-19_dp           ![C]
real(dp), parameter :: eps0  = 8.854187817e-14_dp/elch     ![F/cm] --> [eV/V2/cm] electrical constant of void
complex(dp), parameter :: I = (0, 1)
complex(dp), parameter :: J = (1, 0)
complex(dp), parameter :: O = (0, 0)
end module
